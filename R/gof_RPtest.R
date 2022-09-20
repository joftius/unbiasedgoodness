#' Goodness-of-fit test using RPtest package.
#'
#' @description
#' `gof_RPtest` uses the `RPtest` package to compute a goodness-of-fit test.
#'
#' @details
#' First uses `cv.glmnet` to select a model with `lambda.1se`, then
#' tests the selected model against one with a larger sparsity.
#'
#' @param x A matrix of predictors.
#' @param y A vector of outcomes.
#' @param m Multiplier of null model sparsity.
#' @param b Additive increase to alternative sparsity.
#' @return A p-value computed using `RPtest`.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
#' sim_data <- rXb(n, p, s0)
#' x <- sim_data$x
#' beta <- sim_data$beta
#' y <- x %*% beta + rnorm(n)
#' gof_RPtest(x, y)
#'
gof_RPtest <- function(x, y, m = 1.2, b = 5) {

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
