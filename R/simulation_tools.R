#' Simulation for high-dimensional linear regression model variable selection
#'
#' @description
#' `this` repeats `that` multiple times.
#'
#' @details
#' After generating simulated data this calls `something` for inference.
#'
#' @param niters Number of simulation iterations.
#' @param n Sample size, number of rows in design matrix.
#' @param p Number of random predictor variables or columns in design matrix.
#' @param s0 Sparsity or number of nonzero coefficients in the true linear model.
#' @param DGPargs Other data-generation process arguments passed to `hdi::rXb`.
#' @param fit_fun List of functions to apply to simulated data
#' @param fit_args List of arguments to pass to `fit_fun`
#' @return Variable selection information and p-values for goodness-of-fit tests.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
simulate_hdr <- function(
    niters,
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
                   "toeplitz"  = 1/3,
                   "equi.corr" = 1/20,
                   "exp.decay" = c(0.4, 5)
    ),
    y_fun = y_standard_linear,
    y_args = NULL,
    fit_fun = fit_glmnet_cvmin,
    fit_args = NULL,
    post_fun = post_glmnet_cvmin,
    post_args = NULL,
    cores = NULL
)


{

  xtype <- tolower(match.arg(xtype))
  time_start <- Sys.time()
  local_files <- list.files(path = paste0(getwd(), "/R"), pattern = ".R")
  if (is.null(cores)) cores <- floor(detectCores()/2)
  cl <- makeCluster(cores)
  registerDoParallel(cl)
  clusterCall(cl, function() library(glmnet))
  clusterCall(cl, function() library(hdi))
  clusterCall(cl, function() library(selectiveInference))
  clusterCall(cl, function() library(RPtests))
  #clusterCall(cl, function() library(unbiasedgoodness))
  clusterCall(cl, function() source(paste0(getwd(), "/R/simulation_tools.R")))

  # test_run <- instance_hdr(
  #   n,
  #   p,
  #   s0,
  #   sigma,
  #   xtype = c("toeplitz", "exp.decay", "equi.corr"),
  #   btype = "U[-2,2]",
  #   permuted = TRUE,
  #   x.par = switch(xtype,
  #                  "toeplitz"  = 1/3,
  #                  "equi.corr" = 1/20,
  #                  "exp.decay" = c(0.4, 5)
  #   ),
  #   fit_fun,
  #   fit_args
  # )

  output <- foreach(icount(niters)) %dopar% {
            instance_hdr(
              n,
              p,
              s0,
              xtype = c("toeplitz", "exp.decay", "equi.corr"),
              btype = "U[-2,2]",
              permuted = TRUE,
              x.par = switch(xtype,
                             "toeplitz"  = 1/3,
                             "equi.corr" = 1/20,
                             "exp.decay" = c(0.4, 5)
              ),
              y_fun,
              y_args,
              fit_fun,
              fit_args,
              post_fun,
              post_args
            )
  }
  stopImplicitCluster()
  time_end <- Sys.time()

  print(time_end - time_start)
  return(output)
}


#' Simulation for high-dimensional linear regression model variable selection
#'
#' @description
#' `this` repeats `that` multiple times.
#'
#' @details
#' After generating simulated data this calls `something` for inference.
#'
#' @param niters Number of simulation iterations.
#' @param n Sample size, number of rows in design matrix.
#' @param p Number of random predictor variables or columns in design matrix.
#' @param s0 Sparsity or number of nonzero coefficients in the true linear model.
#' @param DGPargs Other data-generation process arguments passed to `hdi::rXb`.
#' @param fit_fun List of functions to apply to simulated data
#' @param fit_args List of arguments to pass to `fit_fun`
#' @return Variable selection information and p-values for goodness-of-fit tests.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
instance_hdr <- function(
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
                   "toeplitz"  = 1/3,
                   "equi.corr" = 1/20,
                   "exp.decay" = c(0.4, 5)
    ),
    y_fun,
    y_args = NULL,
    fit_fun,
    fit_args = NULL,
    post_fun,
    post_args = NULL)
{
  xtype <- tolower(match.arg(xtype))

  sim_data <- rXb(n = n, p = p, s0 = s0, xtype = xtype,
                  btype = btype, permuted = permuted, x.par = x.par)
  x <- sim_data$x
  beta <- sim_data$beta
  if (is.null(y_args)) {
    y <- y_fun(x, beta)
  } else {
    y <- y_fun(x, beta, y_args)
  }
  # min_beta <- min(abs(beta[1:s0]))
  # max_corr <- max(sapply(1:s0,
  #                        function(j) {
  #                          max(abs(cor(x = x[, j], y = x[, (s0+1):p])))
  #                        }))

  if (is.null(fit_args)) {
    fit_obj <- fit_fun(x, y, beta)
  } else {
    fit_obj <- fit_fun(x, y, beta, fit_args)
  }
  if (is.null(post_args)) {
    return(post_fun(fit_obj, x, y, beta))
  } else {
    return(post_fun(fit_obj, x, y, beta, post_args))
  }
}

y_standard_linear <- function(x, beta, y_args = NULL) {
  sigma <- 1
  if (!is.null(y_args)) sigma <- y_args$sigma
  x %*% beta + sigma * rnorm(nrow(x))
}

fit_cvglmnet <- function(x, y, beta, fit_args = NULL) {
  if (is.null(fit_args)) {
    return(cv.glmnet(x, y))
  } else {
    return(cv.glmnet(x, y, fit_args))
  }
}

post_cvglmnet <- function(fit_obj, x, y, beta, post_args = NULL) {
  if (class(fit_obj) != "cv.glmnet") {
    warning("fit_obj should be of class cv.glmnet")
  }
}

fit_glmnet_cvmin <- function(x, y, beta, fit_args = NULL) {
  if (is.null(fit_args)) {
    cv_fit <- cv.glmnet(x, y)
    glmnet_fit <- glmnet(x, y)
  } else {
    cv_fit <- cv.glmnet(x, y, fit_args)
    glmnet_fit <- glmnet(x, y, fit_args)
  }
  list(cv_fit = cv_fit, glmnet_fit = glmnet_fit)
}

post_glmnet_cvmin <- function(fit_obj, x, y, beta, post_args = NULL) {
  lambda_type <- "lambda.1se"
  if (!is.null(post_args)) lambda_type <- post_args$lambda
  lambda_min <- fit_obj$cv_fit[[lambda_type]]
  list(true_beta = beta,
       beta_hat = coef(fit_obj$glmnet_fit, s = lambda_min))
}

# test run
#source("R/simulation_tools.R")
# one_thing  <- instance_hdr(100, 200, 5, y_fun = y_standard_linear, y_args = NULL, fit_fun = fit_glmnet_cvmin,fit_args = NULL, post_fun = post_glmnet_cvmin,post_args = NULL)
# things  <- simulate_hdr(10, 100, 200, 5)

#####################

#' Simulation for high-dimensional linear model variable selection
#'
#' @description
#' `simulate_gof` repeats `simulate_gof_instance` multiple times.
#'
#' @details
#' After generating simulated data this calls `gof_RPtest` for inference.
#'
#' @param niters Number of simulation iterations.
#' @param n Sample size, number of rows in design matrix.
#' @param p Number of random predictor variables or columns in design matrix.
#' @param s0 Sparsity or number of nonzero coefficients in the true linear model.
#' @param m Multiplier of null model sparsity.
#' @param b Additive increase to alternative sparsity.
#' @param xtype Design matrix correlation structure passed to `hdi::rXb`.
#' @param ... Optional additional parameters passed to `hdi::rXb`.
#' @return Variable selection information and p-values for goodness-of-fit tests.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
simulate_gof <- function(
    niters,
    m,
    b,
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
                   "toeplitz"  = 1/3,
                   "equi.corr" = 1/20,
                   "exp.decay" = c(0.4, 5)
    ))
{
  xtype <- tolower(match.arg(xtype))
  niters |>
    rerun(simulate_gof_instance(m = m,
                                b = b,
                                n = n,
                                p = p,
                                s0 = s0,
                                xtype = xtype,
                                btype = btype,
                                permuted = permuted,
                                x.par = x.par)) |>
    map_df(~data.frame(.x))
}

#' Simulation for high-dimensional linear model variable selection
#'
#' @description
#' `simulate_gof_instance` generates simulated data , applies lasso for variable selection, and computes goodness of fit tests.
#'
#' @details
#' After generating simulated data this calls `gof_RPtest` for inference.
#'
#' @param n Sample size, number of rows in design matrix.
#' @param p Number of random predictor variables or columns in design matrix.
#' @param s0 Sparsity or number of nonzero coefficients in the true linear model.
#' @param m Multiplier of null model sparsity.
#' @param b Additive increase to alternative sparsity.
#' @param ... Optional additional parameters passed to `hdi::rXb`.
#' @return Variable selection information and p-values for goodness-of-fit tests.
#' @export
#' @examples
#' n <- 100
#' p <- 200
#' s0 <- 5
simulate_gof_instance <- function(
    m,
    b,
    n,
    p,
    s0,
    xtype = c("toeplitz", "exp.decay", "equi.corr"),
    btype = "U[-2,2]",
    permuted = TRUE,
    x.par = switch(xtype,
      "toeplitz"  = 1/3,
      "equi.corr" = 1/20,
      "exp.decay" = c(0.4, 5)
      ))
{
  xtype <- tolower(match.arg(xtype))
  sim_data <- rXb(n = n, p = p, s0 = s0, xtype = xtype,
                  btype = btype, permuted = permuted, x.par = x.par)
  x <- sim_data$x
  beta <- sim_data$beta
  min_beta <- min(abs(beta[1:s0]))
  max_corr <- max(sapply(1:s0,
                         function(j) {
                           max(abs(cor(x = x[, j], y = x[, (s0+1):p])))
                  }))

  y <- x %*% beta + rnorm(n)
  gof_fit <- gof_RPtest(x, y)

  larger_support <- union(gof_fit$selected_support, gof_fit$alt_support)
  null_intersection <- length(intersect(gof_fit$selected_support, 1:s0))
  alt_intersection <- length(intersect(larger_support, 1:s0))

  data.frame(
    max_corr = max_corr,
    min_beta = min_beta,
    s_null = length(gof_fit$selected_support),
    s_alt = length(gof_fit$alt_support),
    null_intersection = null_intersection,
    alt_intersection = alt_intersection,
    null = null_intersection == s0,
    alt = alt_intersection == s0,
    pval_lasso = gof_fit$pval_lasso,
    pval_OLS = gof_fit$pval_OLS,
    pval_lasso_rev = gof_fit$pval_lasso_rev,
    pval_OLS_rev = gof_fit$pval_OLS_rev,
    pval_full = gof_fit$pval_full,
    pval_full_lasso = gof_fit$pval_full_lasso,
    pval_OLS_beta = gof_fit$pval_OLS_beta
  )
}
