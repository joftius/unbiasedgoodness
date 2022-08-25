
library(tidyverse)
set.seed(1)
N <- 100
ar1 <- arima.sim(n = N, list(ar = c(0.7)))
ar2 <- arima.sim(n = N, list(ar = c(0.7, 0.2)))
df <- data.frame(t = rep(1:N, 2), Xt = c(ar1, ar2), p = factor(c(rep(1, N), rep(2, N))))

arplot <- ggplot(df, aes(t, Xt)) + geom_line() + facet_grid(p~.) + theme_minimal()

ggsave("ar_paths.pdf", arplot)
