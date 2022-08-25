
library(tidyverse)
library(knitr)

set.seed(1)
n <- 100
p <- 10
active <- 1:5
nz <- 1/8

x <- matrix(rnorm(n*p), nrow = n)
beta <- rep(0, p)
beta[active] <- nz
beta[1:2] <- 1
mu <- x %*% beta

instance <- function() {
  y <- mu + rnorm(n)
  data <- data.frame(y = y, x = x)
  mfull <- lm(y ~ ., data = data)
  moracle <- lm(as.formula(paste("y~", paste(paste0("x.", active[1:2]), collapse="+"))), data = data)
  mselected <- step(mfull, k = log(n), trace = 0)
  testfull <- anova(mselected, mfull)[2,c(3,5,6)]
  testoracle <- anova(moracle, mfull)[2,c(3,5,6)]
  vars <- sort(as.integer(gsub("x.", "", names(mselected$coeff[-1]))))
  return(list(tests = c(testfull, testoracle), vars = vars))
}

starttime <- Sys.time()
simdata <- replicate(10000, instance())
print(Sys.time() - starttime)

results <- cbind(
      t(sapply(simdata[2,], function(vars) c(length(vars), length(intersect(vars, active))))),
      do.call(rbind, simdata[1,]))

colnames(results) <- c("khat", "overlap", "Df.1", "F.1", "Pr.1", "Df.2", "F.2", "Pr.2")

results <- data.frame(results)

df <- results %>% 
  mutate(id = 1:n()) %>%
  gather(key, value, -khat, -overlap, -id) %>%
  separate(key, into = c("var", "alternative")) %>%
  spread(key = var, value) %>% 
  unnest() %>%
  mutate(alternative = recode(as.factor(alternative), `1` = "selected", `2` = "2-sparse"))

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

ggsave(filename = "slides/fstepbic.pdf", fstepbic)

fdists1 <- ggplot(df, aes(Pr)) +
  geom_density() + #aes(color = as.ordered(overlap))) +
  facet_wrap(~ as.factor(alternative) + as.factor(overlap), nrow = 2, scales = "free") +
  xlab("p-value") + ggtitle("Overlap of BIC selected model") +
  theme_minimal()

ggsave(filename = "slides/fdists_both.pdf", fdists1)
 
# nonscreen <- df %>% 
#   filter(overlap < 5)
# 
# ggplot(nonscreen, aes(Pr)) +
#   geom_density() + #aes(color = as.ordered(overlap))) +
#   facet_wrap(~ as.factor(alternative) + as.factor(overlap), nrow = 2)

df_subset <- df %>% filter(alternative == "selected")

fdists2 <- ggplot(df_subset, aes(Pr)) +
  geom_density() + #aes(color = as.ordered(overlap))) +
  facet_wrap(~ as.factor(overlap), nrow = 1) +
  xlab("p-value") + ggtitle("Selected overlap") + theme_minimal()

ggsave(filename = "slides/fdists_selected.pdf", fdists2)

df %>%
  group_by(alternative, overlap) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on size of overlap")

# df %>% 
#   filter(khat < 5, overlap == khat) %>%
#   group_by(alternative, overlap) %>%
#   summarize(rejected = mean(Pr < 0.1))

df %>% 
  filter(overlap < 5) %>%
  group_by(alternative) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on overlap less than 5")

df %>% 
  filter(overlap < 5, overlap == khat) %>%
  group_by(alternative) %>%
  summarize(`Pr(reject)` = mean(Pr < 0.1)) %>%
  kable(digits = 3, format = "latex", caption = "Probability of rejection at level 0.1, conditional on khat = overlap, and overlap less than 5")
