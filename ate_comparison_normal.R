library(MASS)
library(foreach)
library(doParallel)
library(tidyverse)
library(mgcv)
library(latex2exp)

source("estimators.R")

set.seed(2022)
n <- 500
mu <- 1
tau <- 0.5

generate.data <- function(rho, tau) {
  data <-
    mvrnorm(n, mu = c(mu, 0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  Ys1 <- data[, 1]
  Ys1 <- pmin(Ys1, 50)
  Ys1 <- pmax(Ys1, -50)
  
  Xs <- data[, 2]
  ps <- 1 / (1 + exp(-2 * Xs))
  ps <- pmax(ps, 1e-2)
  ps <- pmin(ps, 1 - 1e-6)
  
  indicators <- rbinom(n, 1, ps)
  
  ps.hat <- glm(indicators ~ Xs, family = "binomial")$fitted.values
  
  Ys0 <- Ys1 - tau
  return(
    list(
      "Ys1" = Ys1,
      "Ys0" = Ys0,
      "ps" = ps,
      "ps_hat" = ps.hat,
      "Xs" = Xs,
      "indicators" = indicators
    )
  )
}

sample.estimators <- function(estimators, ...) {
  data <- generate.data(...)
  estimates1 <- rep(0, length(estimators))
  estimates0 <- rep(0, length(estimators))
  
  for (i in seq_along(estimators)) {
    if (substr(estimators[i], 1, 4) == "aipw") {
      estimates1[i] <- do.call(
        estimators[i],
        list(
          Ys = data$Ys1,
          ps = data$ps_hat,
          Xs = data$Xs,
          indicators = data$indicators
        )
      )
      
      estimates0[i] <- do.call(
        estimators[i],
        list(
          Ys = data$Ys0,
          ps = 1 - data$ps_hat,
          Xs = data$Xs,
          indicators = 1 - data$indicators
        )
      )
      
    } else {
      estimates1[i] <- do.call(estimators[i],
                               list(
                                 Ys = data$Ys1,
                                 ps = data$ps,
                                 indicators = data$indicators
                               ))
      estimates0[i] <- do.call(estimators[i],
                               list(
                                 Ys = data$Ys0,
                                 ps = 1 - data$ps,
                                 indicators = 1 - data$indicators
                               ))
    }
  }
  
  #double.fp = double.fp.estimator(data$Ys1, data$Ys0, data$ps, data$indicators)
  return(c(estimates1 - estimates0))
}

estimators <-
  c("hajek.estimator",
    "fp.estimator",
    "aipw.estimator",
    "aipw.fp.estimator")
names <-
  sapply(estimators, function(x)
    str_split(x, ".estim")[[1]][1])
names <- c(names)

N <- 10000
rho.grid <- seq(-0.95, 0.95, 0.05)

registerDoParallel(cores = 4)

start <- Sys.time()
results <- foreach (rho = rho.grid, .combine = "rbind") %dopar% {
  print(rho)

  estimators <- replicate(N, sample.estimators(estimators, rho, tau))

  bias <- rowMeans(estimators-tau)
  vars <- apply(estimators, 1, var)
  mses <- rowMeans((estimators - tau)^2)
  c(bias, vars, mses, rho)
} %>% data.frame()
finish <- Sys.time()

bias_names <- sapply(names, function(x)
  sprintf("%s_bias", x))
var_names <- sapply(names, function(x)
  sprintf("%s_var", x))
mse_names <- sapply(names, function(x)
  sprintf("%s_mse", x))

colnames(results) <-
  unlist(unname(c(bias_names, var_names, mse_names, "rho")))

color_map <- c("aipw" = rgb(214, 39, 40, maxColorValue = 256),
               "aipw.fp" = rgb(148, 103, 189, maxColorValue = 256),
               "hajek" = rgb(31, 119, 180, maxColorValue = 256),
               "fp" = rgb(44, 160, 44, maxColorValue = 256))

name_map_tau <- c("aipw" = unname(TeX('$\\hat{\\tau}_{AIPW}$')),
                  "aipw.fp" = unname(TeX('$\\hat{\\tau}_{AIPW, AN}$')),
                  "hajek" = unname(TeX('$\\hat{\\tau}_{Hajek}$')),
                  "fp" = unname(TeX('$\\hat{\\tau}_{AN}$')))

lty_map <- c("aipw" = "solid",
             "aipw.fp" = "twodash",
             "hajek" = "dashed",
             "fp" = "dotdash")

results %>%
  mutate(hajek_bias = abs(hajek_bias),
         fp_bias = abs(fp_bias),
         aipw_bias = abs(aipw_bias),
         aipw.fp_bias = abs(aipw.fp_bias)) %>% 
  
  pivot_longer(
    cols = unlist(unname(c(
      bias_names, var_names, mse_names
    ))),
    names_sep = "_",
    names_to = c("estimator", "metric")
  ) %>%
  mutate(metric_label = case_when(metric == "mse" ~ "MSE",
                                  metric == "bias" ~ "|bias|",
                                  metric == "var" ~ "Var."),
         metric_label = factor(metric_label, levels = c("MSE", "|bias|", "Var."))) %>% 
  ggplot(aes(
    y = value,
    x = rho,
    linetype = estimator,
    col = estimator
  )) +
  geom_line(lwd = 1) +
  # scale_y_continuous(trans = "log10") +
  labs(x = TeX('Correlation, $\\theta$'),
       y = "Value",
       col = "Estimator",
       linetype = "Estimator") +
  scale_colour_manual(values = color_map, labels = name_map_tau) +
  scale_linetype_manual(values = lty_map, labels = name_map_tau) + 
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = ggplot2::margin(-15, 0, 0, 0),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(vjust = 0),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = grid::unit(c(0, 0, 0, 0), "cm"),
    panel.spacing.x = unit(8, "mm"),
    strip.text.x = element_text(size = 14)            
  ) +
  facet_grid(. ~ metric_label)