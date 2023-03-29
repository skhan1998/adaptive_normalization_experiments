library(sampling)
library(tidyverse)
source("estimators.R")

N <- 100000
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

params <- expand.grid(list(estimators = estimators, response = response, ps = ps), stringsAsFactors = FALSE)

estimate.mses <- function(){
  registerDoParallel(cores = 40)
  results <- foreach (trials = 1:N, .combine = "rbind") %dopar% {
    foreach (i = 1:nrow(params), .combine = "cbind") %do% {
      param <- params[i,]
      Ys <- with(df, get(param$response))
      ps <- with(df, get(param$ps))
      indicators <- rbinom(n, 1, ps)
      do.call(param$estimators, list(Ys = Ys, ps = ps, indicators = indicators))
    }
  } %>% data.frame()
  
  
  colnames(results) <- apply(params, 1, function(x) paste(x, collapse = "_"))
  
  mses <- results %>%
    pivot_longer(cols = everything(), names_to = c("Estimator", "Y", "ps"), names_sep = "_") %>%
    mutate(mean = case_when(Y == "Surfacesbois" ~ mean(df$Surfacesbois),
                            Y == "Airind" ~ mean(df$Airind)),
           error = (value-mean)^2) %>%
    group_by(Estimator, Y, ps) %>%
    summarise(mse = sqrt(mean(error))) %>%
    arrange(Y, ps)
  
  return(mses)
}

estimate.standard.error <- function(){
  results <- foreach (i = 1:10, .combine = "cbind") %do% {
    mses <- estimate.mses()
    mses$mse
  }
  errors <- sqrt(apply(results, 1, var))
  return(errors)
}

error <- estimate.standard.error()
