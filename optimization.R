library(MASS)
library(foreach)
library(doParallel)
library(tidyverse)
library(latex2exp)

source("estimators.R")

set.seed(2021)
n <- 250
mu <- 1
N <- 5000
rho <- -0.5

generate.data <- function(){
  data <- mvrnorm(n, mu = c(mu, 0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  Ys <- data[,1]
  Ys <- pmin(Ys, 50)
  Ys <- pmax(Ys, -50)
  
  ps <- pnorm(data[,2], 0, 1)
  ps <- pmax(ps, 1e-2)
  ps <- pmin(ps, 1-1e-6)
  
  indicators <- rbinom(n, 1, ps)
  return(list("Ys" = Ys, "ps" = ps, "indicators" = indicators))
}

lambda.grid <- seq(-0.5, 2.5, 0.01)

sample.estimators <- function(){
  data <- generate.data()
  Ys <- data$Ys
  ps <- data$ps
  indicators <- data$indicators
  
  qs <- (1-ps)/ps
  S.hat <- sum(Ys*indicators/ps)
  n.hat <- sum(indicators/ps)
  pi.hat <- sum(qs*indicators/ps)/n.hat
  T.hat <- sum(Ys*qs*indicators/ps)/n.hat
  
  mu.hats <- (S.hat)/((1-lambda.grid)*n + lambda.grid*n.hat)
  
  lambda.star.hat <- T.hat/((S.hat/n)*pi.hat)
  lambda.fp.hat <- (T.hat)/(pi.hat*fp.estimator(Ys, ps, indicators))
  
  return(mu.hats)
}

mu.hats <- replicate(N, sample.estimators())

mses <- rowMeans((mu.hats - mu)^2)
vars <- apply(mu.hats, 1, var)
biases <- rowMeans(mu.hats - mu)

plot.df <- data.frame(lambda.grid, mses, vars, biases)

colnames(plot.df) <- c("Lambda", "MSE", "Variance", "Bias")

plot.df %>% 
  pivot_longer(cols = c("MSE", "Variance", "Bias")) %>% 
  mutate(name = factor(name, levels = c("MSE", "Bias", "Variance"))) %>% 
  ggplot(aes(x = Lambda, y = value)) +
  geom_line() + 
  facet_grid(.~name) + 
  theme_bw() + 
  labs(x = "Lambda", y = "Value")

ggsave("optimization.pdf", device = "pdf", width = 6, height = 4, units = "in")