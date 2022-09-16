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
    pval_full = gof_fit$pval_full
  )
}
