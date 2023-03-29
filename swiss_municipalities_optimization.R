library(sampling)
library(tidyverse)
source("estimators.R")

N <- 10000
data("swissmunicipalities")

df <- swissmunicipalities %>%
  select(HApoly, Surfacesbois, Surfacescult, Alp, Airbat, 
         Airind, POPTOT) %>%
  mutate(ps50 = inclusionprobabilities(HApoly, 50),
         ps250 = inclusionprobabilities(HApoly, 250))

n <- nrow(df)
estimators <- c("hajek.estimator", "ht.estimator", "fp.estimator")
response <- c("Surfacesbois", "Airind")
ps <- c("ps50", "ps250")

n <- nrow(df)
lambda.grid <- seq(-0.25, 1.25, 0.01)

registerDoParallel(cores = 40)
results <- foreach (trials = 1:N, .combine = "rbind") %dopar% {
  Ys <- df$Surfacesbois
  ps <- df$ps250
  indicators <- rbinom(n, 1, ps)
  
  S.hat <- sum(Ys*indicators/ps)
  n.hat <- sum(indicators/ps)
  
  S.hat/((1-lambda.grid)*n + lambda.grid*n.hat)
} %>% data.frame()

mses <- colMeans((results-mean(df$Surfacesbois))^2)
plot.df <- data.frame(lambda.grid, mses)
colnames(plot.df) <- c("lambda", "mse")

ggplot(data = plot.df, aes(x = lambda, y = mse)) + 
  geom_line() + 
  scale_y_continuous(trans = "log10") +
  labs(title = unname(TeX('$\\sum p_k=250, Y_1$')), 
       x = "Lambda", y = "MSE") + 
  theme_bw() + 
  theme(axis.title = element_text(size=14),
        axis.text = element_text(size=12),
        plot.margin=grid::unit(c(0.2,0,0,0.1),"cm"))

