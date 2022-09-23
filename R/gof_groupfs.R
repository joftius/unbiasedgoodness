#' Goodness-of-fit test using RPtest package.
#'
#' @description
#' `gof_groupfs` uses `selectiveInference::groupfs` to compute a goodness-of-fit test.
#'
#' @details
#' First uses `groupfs` to select a model via forward stepwise with AIC,
#' then tests the selected model against one with a larger sparsity.
#'
#' @param x A matrix of predictors.
#' @param y A vector of outcomes.
#' @param m Multiplier of null model sparsity.
#' @param b Additive increase to alternative sparsity.
#' @return A p-value computed using `groupfs`.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
#' sim_data <- rXb(n, p, s0)
#' x <- sim_data$x
#' beta <- sim_data$beta
#' y <- x %*% beta + rnorm(n)
#' gof_groupfs(x, y)
gof_groupfs <- function(obj, against, sigma = NULL, verbose = TRUE) {

  n <- nrow(obj$x)
  p <- ncol(obj$x)
  maxsteps <- attr(obj, "maxsteps")
  k <- attr(obj, "k")
  index <- obj$index
  x <- obj$x
  y <- obj$y
  Ep <- sum(index %in% obj$action)

  pvals = dfs = dfs2 = Tstats = numeric(maxsteps)
  supports <- list()
  fullinds <- union(obj$action, against)

  if (!is.null(sigma)) {
    type <- "TC"
    if (!is.null(obj$sigma)) {
      cat(paste("Using specified value", sigma, "for sigma in place of the value", obj$sigma, "used by groupfs()\n"))
    }
  } else {
    if (is.null(obj$sigma)) {
      type <- "TF"
      Pf <- svdu_thresh(obj$x[,which(obj$index %in% fullinds), drop = FALSE])
      dffull <- ncol(Pf)
      df2 <- n - dffull - obj$intercept - 1
      Pfull <- Pf %*% t(Pf)
    } else {
      type <- "TC"
      sigma <- obj$sigma
    }
  }

  i <- against

  if (type == "TC") {
    # Form projection onto active set minus i
    # and project x_i orthogonally
    x_i <- obj$x[,which(obj$index %in% i), drop = FALSE]
    if (length(obj$action) > 1) {
      minus_i <- setdiff(obj$action, i)
      x_minus_i <- svdu_thresh(obj$x[,which(obj$index %in% minus_i), drop = FALSE])
      x_i <- x_i - x_minus_i %*% t(x_minus_i) %*% x_i
    }

    # Project y onto what remains of x_i
    Ugtilde <- svdu_thresh(x_i)
    R <- t(Ugtilde) %*% obj$y
    TC <- sqrt(sum(R^2))
    eta <- Ugtilde %*% R / TC
    Z <- obj$y - eta * TC
    dfi <- ncol(Ugtilde)
    Tstats <- TC
    dfs <- dfi

    ydecomp <- list(Z=Z, eta=eta)

  } else {

    if (length(obj$action) > 1) {
      minus_i <- setdiff(obj$action, i)
      Psub <- svdu_thresh(obj$x[,which(obj$index %in% minus_i), drop = FALSE])
      Z <- Psub %*% t(Psub) %*% obj$y
      df1 <- dffull - ncol(Psub)
    } else {
      Z <- rep(0, n)
      df1 <- dffull + obj$intercept + 1
    }

    C <- df1/df2
    R1 <- obj$y - Z
    R2 <- obj$y - Pfull %*% obj$y
    R1sq <- sum(R1^2)
    R2sq <- sum(R2^2)
    R <- sqrt(R1sq)
    delta <- R1-R2
    Vdelta <- delta/sqrt(sum(delta^2))
    V2 <- R2/sqrt(R2sq)
    TF <- (R1sq-R2sq)/(C*R2sq)
    Tstats <- TF
    dfs <- df1

    ydecomp <- list(R=R, Z=Z, Vd=Vdelta, V2=V2, C=C)

  }

  intervallist <- truncationRegion(obj, ydecomp, type)

  # Additional constraints from AIC stopping
  if (attr(obj, "stopped")) {
    aicintervals <- vector("list", maxsteps)
    aicstop <- attr(obj, "aicstop")
    if (type == "TC") {
      pen0 <- k * obj$intercept
      aic.begin <- aic.last <- sum(obj$y^2)/sigma^2 - n + k * obj$intercept
    } else {
      pen0 <- exp(k * (1+obj$intercept)/n)
      aic.begin <- n*(log(2*pi) + log(mean(obj$y^2))) + k * (1 + n + obj$intercept)
    }
    AICs <- c(aic.begin, obj$AIC)

    ulist <- c(list(matrix(0, n, 1)), obj$maxprojs)
    penlist <- c(pen0, obj$maxpens)
    zlist <- vector("list", maxsteps+1)
    zlist[[1]] <- zlist[[2]] <- Z
    if (type == "TC") {
      etalist <- vector("list", maxsteps+1)
      etalist[[1]] <- etalist[[2]] <- eta
    } else {
      vdlist <- v2list <- vector("list", maxsteps+1)
      vdlist[[1]] <- vdlist[[2]] <- Vdelta
      v2list[[1]] <- v2list[[2]] <- V2
    }
    if (maxsteps > 1) {
      for (step in 1:(maxsteps-1)) {
        cproj <- obj$cumprojs[[step]]
        zlist[[step+2]] <- cproj %*% Z
        if (type == "TC") {
          etalist[[step+2]] <- cproj %*% eta
        } else {
          vdlist[[step+2]] <- cproj %*% Vdelta
          v2list[[step+2]] <- cproj %*% V2
        }
      }
    }

    for (step in 1:maxsteps) {
      # Compare AIC at s+1 to AIC at s
      # roots() functions assume g indexes smaller AIC
      # this is step+1 until the last step
      peng <- penlist[[step+1]]
      Ug <- ulist[[step+1]]
      Uh <- ulist[[step]]
      Zg <- zlist[[step+1]]
      Zh <- zlist[[step]]

      if (type == "TC") {
        penh <- 0
        etag <- etalist[[step+1]]
        etah <- etalist[[step]]
        coeffs <- quadratic_coefficients(obj$sigma, Ug, Uh, peng, penh, etag, etah, Zg, Zh)

        if (AICs[step] < AICs[step+1]) {
          coeffs <- lapply(coeffs, function(coeff) -coeff)
        }

        intstep <- quadratic_roots(coeffs$A, coeffs$B, coeffs$C, tol = 1e-15)

      } else {
        penh <- 1
        Vdg <- vdlist[[step+1]]
        Vdh <- vdlist[[step]]
        V2g <- v2list[[step+1]]
        V2h <- v2list[[step]]
        coeffs <- TF_coefficients(R, Ug, Uh, peng, penh, Zg, Zh, Vdg, Vdh, V2g, V2h)

        if (AICs[step] < AICs[step+1]) {
          coeffs <- lapply(coeffs, function(coeff) -coeff)
        }

        intstep <- TF_roots(R, C, coeffs)
      }

      aicintervals[[step]] <- intstep
    }
    intervallist <- c(intervallist, aicintervals)
  }

  # Compute intersection:
  region <- do.call(interval_union, intervallist)
  region <- interval_union(region, Intervals(c(-Inf,0)))
  E <- interval_complement(region, check_valid = FALSE)

  if (length(E) == 0) {
    stop("Empty support")
  }
  supports <- E

  # E is now potentially a union of intervals
  if (type == "TC") {
    pvals <- TC_surv(TC, sigma, dfi, E)
  } else {
    # write TF_surv function first
    pvals <- TF_surv(TF, df1, df2, E)
  }

  if (any(is.nan(pvals))) {
    nanp <- which(is.nan(pvals))
    pvals[nanp] <- 0
    warning(paste0("P-value NaNs of the form 0/0 converted to 0 for group(s) ", paste(obj$action[nanp], collapse=","), ". This typically occurs for numerical reasons in the presence of a large signal-to-noise ratio."))
  }

  names(pvals) <- "Goodness of fit"
  out <- list(vars = obj$action, pv=pvals)
  if (type == "TC") {
    out$TC <- Tstats
    out$sigma <- sigma
  } else {
    out$TF <- Tstats
    out$df2 <- df2
  }
  out$df <- dfs
  out$support <- supports
  class(out) <- "gof_groupfs"
  if (!is.null(attr(obj, "varnames"))) {
    attr(out, "varnames") <- attr(obj, "varnames")
  }
  return(out)
}

print.gof_groupfs <- function(obj) {
  out <- data.frame(df = obj$df, TF = obj$TF, Min = min(obj$support), Max = max(obj$support), pv = obj$pv)
  rownames(out) <- NULL
  colnames(out)[5] <- "Pr(>F)"
  print(out)
  invisible()
}


another_gof_RPtest <- function(x, y, m = 1.2, b = 5) {

  cv_fit <- glmnet::cv.glmnet(x, y)
  lasso_fit <- glmnet::glmnet(x, y)

  beta_min <- coef(lasso_fit, s = cv_fit$lambda.1se)
  sparsity_min <- sum(beta_min != 0) - 1

  # RPtest cannot handle sparsity == 1
  if (sparsity_min < 2) {
    warning("Sparsity of 1 chosen by cv.glmnet, forcing to 2 for RPtest compatibility\n")
    lambda_2 <- lasso_fit$lambda[which(lasso_fit$df >= 2)[1]]
    beta_min <- coef(lasso_fit, s = lambda_2)
    sparsity_min <- sum(beta_min != 0) - 1
  }

  # Try to get an alternative model with desired larger sparsity
  larger_sparsity <- ceiling(m * sparsity_min + b)
  larger_lambda_indices <- which(lasso_fit$df >= larger_sparsity)
  if (length(larger_lambda_indices) >= 1) {
    larger_fit_lambda <- lasso_fit$lambda[larger_lambda_indices[1]]
    #print(larger_fit_lambda)
  } else {
    warning("Desired larger support size not attained, trying m = 1\n")
    second_larger_sparsity <- sparsity_min + b
    second_larger_lambda_indices <- which(lasso_fit$df >= second_larger_sparsity)
    if (length(second_larger_lambda_indices) >= 1) {
      larger_fit_lambda <- lasso_fit$lambda[second_larger_lambda_indices[1]]
      #print(larger_fit_lambda)
    } else {
      warning("Larger support using m = 1 also failed, using lambda.min/10\n")
      larger_fit_lambda <- cv_fit$lambda.min / 10
    }
    #print(c(cv_fit$lambda.1se, larger_fit_lambda))
  }

  beta_larger <- coef(lasso_fit, s = larger_fit_lambda)
  sparsity_larger <- sum(beta_larger != 0) - 1

  support_min <- which(beta_min[-1] != 0)
  support_larger <- which(beta_larger[-1] != 0)
  support_new <- setdiff(support_larger, support_min)

  if (length(support_new) == 0) {
    #print(c(sparsity_min, larger_sparsity))
    stop("Failed to include additional variables in larger model\n")
  }

  pval_full <- RPtests::RPtest(x[, support_min, drop = FALSE],
                              y,
                              resid_type = "OLS",
                              test = "group",
                              x_alt = x[, -support_min, drop = FALSE])

  pval_full_lasso <- RPtests::RPtest(x[, support_min, drop = FALSE],
                               y,
                               resid_type = "Lasso",
                               test = "group",
                               x_alt = x[, -support_min, drop = FALSE])

  # pass estimated beta
  pval_OLS_beta <- RPtests::RPtest(x[, support_min, drop = FALSE],
                              y,
                              resid_type = "OLS",
                              test = "group",
                              x_alt = x[, support_new, drop = FALSE],
                              beta_est = beta_min)

  pval_OLS <- RPtests::RPtest(x[, support_min, drop = FALSE],
                          y,
                          resid_type = "OLS",
                          test = "group",
                          x_alt = x[, support_new, drop = FALSE])

  pval_lasso <- RPtests::RPtest(x[, support_min, drop = FALSE],
                 y,
                 resid_type = "Lasso",
                 test = "group",
                 x_alt = x[, support_new, drop = FALSE])

  # these _rev ones are nested, should probably be removed
  suppressWarnings(
  pval_OLS_rev <- RPtests::RPtest(x[, support_larger, drop = FALSE],
                                  y,
                                  resid_type = "OLS",
                                  test = "group",
                                  x_alt = x[, support_min, drop = FALSE])
  )
  suppressWarnings(
  pval_lasso_rev <- RPtests::RPtest(x[, support_larger, drop = FALSE],
                                y,
                                resid_type = "Lasso",
                                test = "group",
                                x_alt = x[, support_min, drop = FALSE])
  )
  list(selected_support = support_min,
       alt_support = support_new,
       pval_full = pval_full,
       pval_full_lasso = pval_full_lasso,
       pval_lasso = pval_lasso,
       pval_OLS = pval_OLS,
       pval_OLS_beta = pval_OLS_beta,
       pval_lasso_rev = pval_lasso_rev,
       pval_OLS_rev = pval_OLS_rev)
}
