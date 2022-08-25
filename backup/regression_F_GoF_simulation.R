
library(tidyverse)
library(knitr)
library(selectiveInference)
source("../R-software/selectiveInference/R/funs.common.R")
source("../R-software/selectiveInference/R/funs.groupfs.R")
source("../R-software/selectiveInference/R/funs.quadratic.R")

# New groupfsGoF code
# ---- 

groupfsGoF <- function(obj, against, sigma = NULL, verbose = TRUE) {
  
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
  class(out) <- "groupfsGoF"
  if (!is.null(attr(obj, "varnames"))) {
    attr(out, "varnames") <- attr(obj, "varnames")
  }
  return(out)
}

print.groupfsGoF <- function(obj) {
  out <- data.frame(df = obj$df, TF = obj$TF, Min = min(obj$support), Max = max(obj$support), pv = obj$pv)
  rownames(out) <- NULL
  colnames(out)[5] <- "Pr(>F)"
  print(out)
  invisible()
}

# ---- 

# Simulation code

# ----
# Low-dimension simulation code
# ----
set.seed(1)
n <- 100
p <- 10
active <- 1:5
nz <- 1/5

x <- matrix(rnorm(n*p), nrow = n)
beta <- rep(0, p)
beta[active] <- nz
beta[1:2] <- 1
mu <- x %*% beta

instance <- function() {
  y <- mu + rnorm(n)
  ybar <- mean(y)
  y <- y - ybar
  data <- data.frame(y = y, x = x)
  fit <- groupfs(x, y, index = 1:ncol(x), intercept = FALSE, center = FALSE, normalize = FALSE, aicstop = 1, k = log(n))
  against <- setdiff(fit$index, fit$action)
  adjpv <- groupfsGoF(fit, against)
  mselected <- lm(as.formula(paste("y~-1+", paste(paste0("x.", fit$action), collapse="+"))), data = data)
  mfull <- lm(y ~ .-1, data = data)
  unadjpv <- anova(mselected, mfull)[2,c(3,5,6)]
  vars <- sort(fit$action)
  return(list(tests = c(unadjpv, c(adjpv[c("df", "TF", "pv")],
                                   min(adjpv$support),
                                   max(adjpv$support))),
              vars = vars))
}


# ----
# High-dimension simulation code
# ----
set.seed(1)
n <- 100
p <- 200
active <- 1:5
nz <- 1/5

x <- matrix(rnorm(n*p), nrow = n)
beta <- rep(0, p)
beta[active] <- nz
beta[1:3] <- 1
mu <- x %*% beta

instance <- function() {
  y <- mu + rnorm(n)
  ybar <- mean(y)
  y <- y - ybar
  data <- data.frame(y = y, x = x)
  fit <- groupfs(x, y, index = 1:ncol(x), intercept = FALSE, center = FALSE, normalize = FALSE, aicstop = 1, k = log(n))
  against <- setdiff(fit$index, fit$action)
  adjpv <- groupfsGoF(fit, against)
  mselected <- lm(as.formula(paste("y~-1+", paste(paste0("x.", fit$action), collapse="+"))), data = data)
  mfull <- lm(y ~ .-1, data = data)
  unadjpv <- anova(mselected, mfull)[2,c(3,5,6)]
  vars <- sort(fit$action)
  return(list(tests = c(unadjpv, c(adjpv[c("df", "TF", "pv")],
                                   min(adjpv$support),
                                   max(adjpv$support))),
              vars = vars))
}


#### This actually runs the simulation
# When code is working correctly change from 100 below to ~10000 or something large
# for the final simulation to include in the slides/paper
starttime <- Sys.time()
simdata <- replicate(100, instance())
print(Sys.time() - starttime)
###

# ----
# Storing/plotting results
# ----
results <- cbind(
  t(sapply(simdata[2,], function(vars) c(length(vars), length(intersect(vars, active))))),
  do.call(rbind, simdata[1,]))

colnames(results) <- c("khat", "overlap", "Df.1", "F.1", "Pr.1", "Df.2", "F.2", "Pr.2", "min", "max")

results <- data.frame(results)

summary(unlist(results$min))
summary(unlist(results$max))

df <- results %>% 
  mutate(id = 1:n()) %>%
  gather(key, value, -khat, -overlap, -id, -min, -max) %>%
  separate(key, into = c("var", "pvalue")) %>%
  spread(key = var, value) %>% 
  unnest() %>%
  mutate(pvalue = recode(as.factor(pvalue), `1` = "Naive", `2` = "Adjusted"))


models <- table(df$khat, df$overlap)/2
tb <- expand.grid(rownames(models), colnames(models))
tb$value <- sapply(1:nrow(tb), function(i) models[tb[i,1], tb[i,2]])

fstepbic <- ggplot(tb, aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  geom_text(aes(label = value), color="white") +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_discrete(expand = c(0,0)) +
  scale_fill_gradient("Instances", low = "white", high = "black") +
  xlab("Selected sparsity") + ylab("True positives") +
  theme_minimal()


ggsave(filename = "slides/fstepbic2.pdf", fstepbic)

fdists1 <- ggplot(df, aes(Pr)) +
  geom_density() + #aes(color = as.ordered(overlap))) +
  facet_wrap(~ as.factor(pvalue) + as.factor(overlap), nrow = 2, scales = "free") +
  xlab("p-value") + ggtitle("Overlap of BIC selected model") +
  theme_minimal()

ggsave(filename = "slides/fdists_both2.pdf", fdists1)

df_subset <- df %>% filter(pvalue == "Adjusted")

fdists2 <- ggplot(df_subset, aes(Pr)) +
  geom_density() + #aes(color = as.ordered(overlap))) +
  facet_wrap(~ as.factor(overlap), nrow = 1) +
  xlab("p-value") + ggtitle("Selected overlap") + theme_minimal()

ggsave(filename = "slides/fdists_selected2.pdf", fdists2)

df %>%
  group_by(pvalue, overlap) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on size of overlap")

df %>%
  filter(overlap == khat) %>%
  group_by(pvalue, overlap) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on size of overlap")


# df %>% 
#   filter(khat < 5, overlap == khat) %>%
#   group_by(alternative, overlap) %>%
#   summarize(rejected = mean(Pr < 0.1))

df %>% 
  filter(overlap < 5) %>%
  group_by(pvalue) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on overlap less than 5")

df %>% 
  filter(overlap <= 4, overlap == khat) %>%
  group_by(pvalue) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on khat = overlap, and overlap less than 5")

# ----