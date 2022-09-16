## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
set.seed(1)
x <- runif(50)
y <- x + rnorm(50)
summary(lm(y ~ x))

## ----setup--------------------------------------------------------------------
library(RPtests)
library(hdi) # riboflavin dataset
library(purrr)
library(tidyverse)
devtools::load_all() # this shouldn't be necessary?
library(unbiasedgoodness)

## -----------------------------------------------------------------------------
set.seed(1)
n <- 100
p <- 200
s0 <- 5
sim_data <- rXb(n, p, s0, xtype = "equi.corr")
x <- sim_data$x
beta <- sim_data$beta
y <- x %*% beta + rnorm(n)
# True support
which(beta != 0)

## -----------------------------------------------------------------------------
gof_RPtest(x, y)

## -----------------------------------------------------------------------------
devtools::load_all() # this shouldn't be necessary?
library(unbiasedgoodness)
set.seed(1)
n <- 100
p <- 200
s0 <- 5
m <- 1.2
b <- 5

simulate_gof_instance(m, b = 10, n, p, s0)

## -----------------------------------------------------------------------------
devtools::load_all() # this shouldn't be necessary?
library(unbiasedgoodness)
set.seed(1)
niters <- 20
output <- simulate_gof(niters, m, b = 10, n, p, s0)

## -----------------------------------------------------------------------------
output |> 
  summarize(null = mean(null),
            alt = mean(alt),
            s_null = mean(s_null),
            s_alt = mean(s_alt))

## -----------------------------------------------------------------------------
output |>
  ggplot(aes(pval_OLS)) +
  geom_histogram() +
  facet_grid(~null)

## -----------------------------------------------------------------------------
output |>
  ggplot(aes(pval_full)) +
  geom_histogram() +
  facet_grid(~null)

## -----------------------------------------------------------------------------
output |>
  ggplot(aes(pval_lasso)) +
  geom_histogram() +
  facet_grid(~null)

## -----------------------------------------------------------------------------
library(tidyverse)

## -----------------------------------------------------------------------------
output %>% pull(null_intersection) %>% qplot()
output %>% pull(alt_intersection) %>% qplot()

## -----------------------------------------------------------------------------
output |> pull(s_alt) |> qplot()

## -----------------------------------------------------------------------------
output |> ggplot(aes(max_corr, null_intersection)) + 
  geom_point() + geom_smooth(method = "lm")

## -----------------------------------------------------------------------------
output |> ggplot(aes(min_beta, null_intersection)) + 
  geom_point() + geom_smooth(method = "lm")

## -----------------------------------------------------------------------------
output %>% pull(null_intersection) |> qplot()

