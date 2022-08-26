#' Goodness-of-fit test using RPtest package.
#'
#' @description
#' `gof_RPtest` uses the `RPtest` package to compute a goodness-of-fit test.
#'
#' @details
#' First uses `cv.glmnet` to select a model with `lambda.min`, then
#' tests the selected model against one with a value of `lambda`
#' chosen to increase the model sparsity.
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
gof_RPtest <- function(x, y, m = 2, b = 1) {
  cv_fit <- glmnet::cv.glmnet(x, y)
  lasso_fit <- glmnet::glmnet(x, y)
  beta_min <- coef(lasso_fit, s = cv_fit$lambda.min)
  sparsity_min <- sum(beta_min != 0) - 1

  # RPtest cannot handle sparsity = 1
  if (sparsity_min < 2) {
    lambda_2 <- lasso_fit$lambda[which(lasso_fit$df >= 2)[1]]
    beta_min <- coef(lasso_fit, s = lambda_2)
    sparsity_min <- sum(beta_min != 0) - 1
  }

  larger_sparsity <- m * sparsity_min + b
  larger_fit_lambda <- lasso_fit$lambda[which(lasso_fit$df >= larger_sparsity)[1]]
  beta_larger <- coef(lasso_fit, s = larger_fit_lambda)
  sparsity_larger <- sum(beta_larger != 0) - 1
  support_min <- which(beta_min[-1] != 0)
  support_larger <- which(beta_larger[-1] != 0)
  support_new <- setdiff(support_larger, support_min)
  pval <- RPtests::RPtest(x[, support_min, drop = FALSE],
                 y,
                 resid_type = "Lasso",
                 test = "group",
                 x_alt = x[, support_new, drop = FALSE])
  list(selected_support = support_min, alt_support = support_new, pval = pval)
}
