library(MASS)
library(foreach)
library(doParallel)
library(tidyverse)
library(latex2exp)

source("estimators.R")

set.seed(2021)
mu <- 0.1
N <- 200000
rho <- 0.9

generate.data <- function(n, delta) {
  data <-
    mvrnorm(n, mu = c(mu, 0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  Ys <- data[, 1]
  Ys <- pmin(Ys, 50)
  Ys <- pmax(Ys,-50)
  
  Xs <- data[, 2]
  ps <- 1 / (1 + exp(-2 * Xs))
  ps <- pmax(ps, delta)
  ps <- pmin(ps, 1 - 1e-6)
  
  indicators <- rbinom(n, 1, ps)
  return(list(
    "Ys" = Ys,
    "ps" = ps,
    "indicators" = indicators
  ))
}

lambda.grid <- c(0.25, 0.5, 0.75)

sample.estimators <- function(...) {
  data <- generate.data(...)
  Ys <- data$Ys
  ps <- data$ps
  indicators <- data$indicators
  
  qs <- (1 - ps) / ps
  S.hat <- sum(Ys * indicators / ps)
  n.hat <- sum(indicators / ps)
  pi.hat <- sum(qs * indicators / ps) / n.hat
  T.hat <- sum(Ys * qs * indicators / ps) / n.hat
  
  mu.hats <-
    (S.hat) / ((1 - lambda.grid) * length(Ys) + lambda.grid * n.hat)
  
  mu.hat.an <- fp.estimator(Ys, ps, indicators)
  return(c(mu.hats, mu.hat.an))
}

delta.grid = c(0.01, 0.05, 0.1)
n.grid = seq(10, 200, 5)

registerDoParallel(cores = 40)
results <- foreach (n = n.grid, .combine = "rbind") %dopar% {
  bias_estimates = foreach (delta = delta.grid, .combine = "cbind") %do% {
    estimators <- replicate(N, sample.estimators(n, delta))
    rowMeans(estimators - mu)
  }
  
  c(bias_estimates, n)
} %>% data.frame()

colnames(results) = c(apply(expand.grid(c(lambda.grid, "AN"), delta.grid), 1, function(x)
  paste(x, collapse = "_")), "n")
  

results %>% 
  pivot_longer(cols = -n, names_sep = "_", names_to = c("lambda", "delta")) %>%
  filter(n>20) %>% 
  mutate(value = abs(value)) %>% 
  ggplot(aes(x = n, y = value, color = delta)) + 
  geom_line() + 
  scale_y_continuous(trans = "log10") +
  labs(x = "Sample size, n", y = "|bias|", color = "Delta") + 
  facet_grid(lambda~., scale = "free_y") + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = ggplot2::margin(-5, 0, 0, 0),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(vjust = 0.5),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = grid::unit(c(0.1, 0.1, 0, 0), "cm"),
    panel.spacing.x = unit(8, "mm"),
    strip.text.y = element_text(size = 14)            
  )