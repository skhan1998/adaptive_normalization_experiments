library(MASS)
library(foreach)
library(doParallel)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(mgcv)
source("estimators.R")
library(logger)

args = commandArgs(trailingOnly=TRUE) 

set.seed(args[2])

mu <- 1
rho <- 0.1

generate.data <- function(n) {
  data <-
    mvrnorm(n, mu = c(mu, 0), Sigma = matrix(c(1, rho, rho, 1), nrow = 2))
  Ys <- data[, 1]
  Ys <- pmin(Ys, 50)
  Ys <- pmax(Ys, -50)
  
  Xs <- data[, 2]
  ps <- 1 / (1 + exp(-2 * Xs))
  ps <- pmax(ps, 1e-2)
  ps <- pmin(ps, 1 - 1e-6)
  

  indicators <- rbinom(n, 1, ps)
  
  return(list(
    "Ys" = Ys,
    "ps" = ps,
    "Xs" = Xs,
    "indicators" = indicators,
    "mu" = mu
  ))
}

sample.estimators <- function(estimators, n) {
  data <- generate.data(n)
  estimates <- rep(0, length(estimators))
  
  for (i in seq_along(estimators)) {
    if (substr(estimators[i], 1, 4) == "aipw") {
      estimates[i] <- do.call(
        estimators[i],
        list(
          Ys = data$Ys,
          ps = data$ps,
          Xs = data$Xs,
          indicators = data$indicators
        )
      )
    } else {
      estimates[i] <- do.call(estimators[i],
                              list(
                                Ys = data$Ys,
                                ps = data$ps,
                                indicators = data$indicators
                              ))
    }
  }
  return(estimates)
}

estimators <-
  c(
    "fp.estimator",
    "hajek.estimator",
    "ht.estimator",
    "aipw.estimator",
    "aipw.fp.estimator",
    "fp.var.estimator",
    "hajek.var.estimator",
    "ht.var.estimator",
    "aipw.var.estimator"
  )

N <- 20000
n.grid <- seq(100, 750, 10)

registerDoParallel(cores = 10)
results <- foreach (n = n.grid, .combine = "rbind") %dopar% {
  log_info(sprintf("Sampling for n={n}"))
  estimators <- replicate(N, sample.estimators(estimators, n))


  fp.cover = mu < estimators[1,] + 1.96 * sqrt(estimators[6,]) / sqrt(n)
  fp.cover = fp.cover &
    (estimators[1,] - 1.96 * sqrt(estimators[6,]) / sqrt(n) < mu)

  hajek.cover = mu < estimators[2,] + 1.96 * sqrt(estimators[7,]) / sqrt(n)
  hajek.cover = hajek.cover &
    (estimators[2,] - 1.96 * sqrt(estimators[7,]) / sqrt(n) < mu)

  ht.cover = mu < estimators[3,] + 1.96 * sqrt(estimators[8,]) / sqrt(n)
  ht.cover = hajek.cover &
    (estimators[3,] - 1.96 * sqrt(estimators[8,]) / sqrt(n) < mu)

  aipw.cover = mu < estimators[4,] + 1.96 * sqrt(estimators[9,]) / sqrt(n)
  aipw.cover = aipw.cover &
    (estimators[4,] - 1.96 * sqrt(estimators[9,]) / sqrt(n) < mu)

  aipw.fp.cover = mu < estimators[5,] + 1.96 * sqrt(estimators[9,]) /
    sqrt(n)
  aipw.fp.cover = aipw.fp.cover &
    (estimators[5,] - 1.96 * sqrt(estimators[9,]) / sqrt(n) < mu)

  c(
    mean(fp.cover),
    mean(hajek.cover),
    mean(ht.cover),
    mean(aipw.cover),
    mean(aipw.fp.cover),
    var(estimators[1,]),
    var(estimators[2,]),
    var(estimators[3,]),
    var(estimators[4,]),
    var(estimators[5,]),
    mean(estimators[6, ])/n,
    mean(estimators[7, ])/n,
    mean(estimators[8, ])/n,
    mean(estimators[9, ])/n,
    n
  )
} %>% data.frame()

colnames(results) <- c("fp_cover", "hajek_cover", "ht_cover", "aipw_cover", "aipw.fp_cover",
                       "fp_var", "hajek_var", "ht_var", "aipw_var", "aipw.fp_var",
                       "fp_asy", "hajek_asy", "ht_asy", "aipw_asy",
                       "n")

color_map <- c("aipw" = rgb(214, 39, 40, maxColorValue = 256),
               "aipw.fp" = rgb(148, 103, 189, maxColorValue = 256),
               "hajek" = rgb(31, 119, 180, maxColorValue = 256),
               "fp" = rgb(44, 160, 44, maxColorValue = 256))

name_map_tau <- c("aipw" = unname(TeX('$\\hat{\\mu}_{AIPW}$')),
                  "aipw.fp" = unname(TeX('$\\hat{\\mu}_{AIPW, AN}$')),
                  "hajek" = unname(TeX('$\\hat{\\mu}_{Hajek}$')),
                  "fp" = unname(TeX('$\\hat{\\mu}_{AN}$')))

lty_map <- c("aipw" = "solid",
             "aipw.fp" = "twodash",
             "hajek" = "dashed",
             "fp" = "dotdash")

plot_df = results %>% pivot_longer(cols = c(-n), names_sep = "_", names_to = c("estimator", "quantity"))  %>%
  filter(n > 200)

plot_df %>%
  ggplot(aes(y = value, x = n, color = estimator, lty = estimator)) +
  geom_line(lwd=1) +
  labs(x = TeX('Sample size, $n$'),
       y = "Value",
       col = "Estimator",
       lty = "Estimator") +
  scale_colour_manual(values = color_map, labels = name_map_tau) +
  scale_linetype_manual(values = lty_map, labels = name_map_tau) + 
  geom_hline(data = plot_df %>% filter(quantity == "cover"), aes(yintercept = 0.95), linetype = "dashed") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.margin = ggplot2::margin(-10, 0, 0, 0),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    plot.margin = grid::unit(c(0.2, 0, 0, 0.1), "cm"),
    strip.text.y = element_text(size = 14)            
  ) +
  facet_grid(case_when(quantity == "asy" ~ "Estimated var.",
                                  quantity == "cover" ~ "Coverage",
                                  quantity == "var" ~ "Empirical var.") ~ ., 
             scales = "free_y")
