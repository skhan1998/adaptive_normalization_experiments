source("estimators.R")
library(MASS)
library(foreach)
library(doParallel)
library(tidyverse)

set.seed(2021)

mu <- 1

generate.data <- function(n, rho = 0.25){
  data <- mvrnorm(n, mu = c(mu, 0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  Ys <- data[,1]
  ps <- pnorm(data[,2], 0, 1)
  indicators <- rbinom(n, 1, ps)
  return(list("Ys" = Ys, "ps" = ps, "indicators" = indicators, "mu" = mu))
}

sample.estimators <- function(n, ...){
  data <- generate.data(n, ...)
  mu.hat.fp <- fp.estimator(data$Ys, data$ps, data$indicators)

  iterations <- lambda.iterative.estimator(data$Ys, data$ps, data$indicators, 50)
  return(c(mu.hat.fp, iterations))
}

N <- 10000
n <- 100

estimators <- replicate(N, sample.estimators(n))

iters <- (1:(nrow(estimators)-1))-1
vars <- apply(estimators, 1, var)[2:nrow(estimators)]

ggplot(data = data.frame(x = iters, y = vars), aes(x=x, y=y)) + 
  geom_line() + scale_y_continuous(trans = "log10") + 
  labs(x = "Number of iterations, t", y = "Variance") + theme_bw()

ggsave("iterative_variance.pdf", device = "pdf",
    width = 6, height = 4, units = "in")

