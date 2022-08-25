

library(tidyverse)

set.seed(1)
C <- 2

Z <- sort(rnorm(1000000))

cutoff <- qnorm(.975)
tcutoff <- quantile(abs(Z[abs(Z) < C]), .95)

#mean(abs(Z[abs(Z) < C]) > tcutoff)

mu <- seq(from = -2, to = 2, length.out = 100)

Zs <- outer(Z, mu, FUN = "+")

tpow <- apply(Zs, 2, function(col) mean(abs(col[abs(col) < C]) > tcutoff))
pow <- apply(Zs, 2, function(col) mean(abs(col[abs(col) < C]) > cutoff))

df <- data.frame(mu = rep(mu, 2), power = c(pow, tpow), cutoff = c(rep("unconditional", 100), rep("conditional", 100)))


powerm1 <- ggplot(df, aes(mu, power)) + 
  geom_line(aes(linetype = cutoff)) +
  ylab("conditional power") + theme_minimal() +
  ggtitle("Hard thresholding with C = 2")

ggsave(filename = "slides/powerm1.pdf", powerm1)