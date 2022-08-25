
tdnorm <- function(x, a, b, mean = 0, sd = 1) {
  if ((x < a) | (x > b)) return(0)
  dnorm(x, mean, sd)/(pnorm(b, mean, sd) - pnorm(a, mean, sd))
}

tpnorm <- function(x, a, b, mean = 0, sd = 1) {
  if (x < a) return(0)
  if (x > b) return(1)
  (pnorm(x, mean, sd) - pnorm(a, mean, sd))/(pnorm(b, mean, sd) - pnorm(a, mean, sd))
}

tnmean <- function(a, b, mean = 0, sd = 1) {
  A <- (a - mean)/sd
  B <- (b - mean)/sd
  mean + sd*(dnorm(A) - dnorm(B))/(pnorm(B) - pnorm(A))
}

tnvar <- function(a, b, mean = 0, sd = 1) {
  A <- (a - mean)/sd
  B <- (b - mean)/sd
  BA <- pnorm(B) - pnorm(A)
  sd^2*(1 + (A*dnorm(A) - B*dnorm(B))/BA - ((dnorm(A) - dnorm(B))/BA)^2)
}

tnsqmean <- function(a, b, mean = 0, sd = 1) {
  A <- (a - mean)/sd
  B <- (b - mean)/sd
  BA <- pnorm(B) - pnorm(A)
  mean^2 + 2*mean*sd*(dnorm(A) - dnorm(B))/BA + sd^2*(1 + (A*dnorm(A) - B*dnorm(B))/BA)
}


sd <- 3
mean <- -1

Z <- sd*rnorm(10000) + mean

c(mean(Z[abs(Z) < 1]), tnmean(-1,1,mean,sd))
c(var(Z[abs(Z) < 1]), tnvar(-1,1,mean,sd))
c(mean(Z[abs(Z) < 1]^2), tnsqmean(-1,1,mean,sd))





Z1 <- rnorm(900000)
Z2 <- rnorm(900000)
Z3 <- rnorm(900000)
Z4 <- rnorm(900000)
Z5 <- rnorm(900000)
Z6 <- rnorm(900000)
Z7 <- rnorm(900000)

Z1s <- Z1[abs(Z1) < 2]
Z2s <- Z2[abs(Z2) < 2]
Z3s <- Z3[abs(Z3) < 2]
Z4s <- Z4[abs(Z4) < 2]
Z5s <- Z5[abs(Z5) < 2]
Z6s <- Z6[abs(Z6) < 2]
Z7s <- Z7[abs(Z7) < 2]

ml <- min(length(Z1s), length(Z2s)) #, length(Z3s)) , length(Z4s), length(Z5s), length(Z6s), length(Z7s))


Zs <- Z1s[1:ml]^2 + Z2s[1:ml]^2 #+ Z3s[1:ml]^2 + Z4s[1:ml]^2 + Z5s[1:ml]^2 + Z6s[1:ml]^2 + Z7s[1:ml]^2

library(tidyverse)

tdchisq <- function(x, df, C) {
  
  output <- dchisq(x, df)/pchisq(C^2, 1)^df
  
  #output <- (pchisq(C^2, 1)^df - ((pnorm(C) - pnorm(-C))^df))/2
  
  output[x > C^2] <- pchisq(C^2, 1)^df - ((pnorm(C) - pnorm(-C))^df - output[x > C^2])
  
  #output <- dchisq(x, df)/pchisq(df*C^2, df)
  output[x > df*C^2] <- 0
  return(output)
}

ggplot(data.frame(Zs), aes(Zs)) + stat_density() +
  stat_function(fun = tdchisq, args = list(df = 2, C = 2)) +
  stat_function(fun = dchisq, args = list(df = 2), color = "orange") 



Zcrit <- 7*tnsqmean(-2,2) + 1.96*sqrt(7*tnvar(-2,2))

mean(Zs > Zcrit)

quantile((Zs/7 - tnsqmean(-1,1))/sqrt(tnvar(-2,2)/7), .95)


