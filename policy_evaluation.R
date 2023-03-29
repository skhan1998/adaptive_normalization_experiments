library(foreach)
library(doParallel)
library(MASS)

generate.data <- function(n){
  Xs <- mvrnorm(n, mu = c(0, 0, 0), Sigma = matrix(c(1, 0, 0,
                                                     0, 1, 0,
                                                     0, 0, 1), nrow = 3))
  ps <- 1/(1+exp(-Xs[,1]))
  
  Ys1 <- Xs[,1]
  
  Ys0 <- Ys1 - sign(Xs[,2]+Xs[,3])
  
  indicators <- rbinom(n, 1, ps)
  return(list(Ys1=Ys1, Ys0=Ys0, Xs = Xs, ps = ps, indicators = indicators))
}

learn.policy <- function(Ys1, Ys0, Xs, ps, indicators){
  
  #estimated value of policy that thresholds X[,2] at t
  value <- function(t){
    treatments = rep(1, length(Ys1))
    treatments[Xs[,2] < t] = 0 
    Ys <- Ys1*indicators + Ys0*(1-indicators)
    
    obs = indicators == treatments
    obs_ps = treatments*ps + (1-treatments)*(1-ps)
    obs_qs = (1-obs_ps)/(obs_ps)
    
    V.ipw <- mean((obs*Ys)/(obs_ps))
    
    V.ipw.an <- V.ipw + (sum(Ys*obs_qs*obs/obs_ps)/sum(obs_qs*obs/obs_ps))*(1-mean(obs/obs_ps))
    
    return(c(V.ipw, V.ipw.an))
  }
  
  t.grid <- seq(-1, 1, 0.005)
  values <- sapply(t.grid, value)
  
  t.ipw <- t.grid[which.max(values[1,])]
  t.ipw.an <- t.grid[which.max(values[2,])]  
  return(c(t.ipw, t.ipw.an))
}

learn.policy.aipw <- function(Ys1, Ys0, Xs, ps, indicators){
  Ys <- Ys1*indicators + Ys0*(1-indicators)
  
  df <- data.frame(Ys1, Ys0, Xs)
  
  model1 <- lm(Ys1 ~ Xs, subset = (indicators==1))
  model0 <- lm(Ys0 ~ Xs, subset = (indicators==0))
  
  Ys1.hat <- predict(model1, df)
  Ys0.hat <- predict(model0, df)
  

  
  #estimated value of policy that thresholds X[,2] at t
  value <- function(t){
    treatments = rep(1, length(Ys1))
    treatments[Xs[,2] < t] = 0 
    
    weights <- Ys1.hat - Ys0.hat 
    weights <- weights + indicators*(Ys1-Ys1.hat)/ps
    weights <- weights + (1-indicators)*(Ys0-Ys0.hat)/(1-ps)
    
    A.hat.aipw <- sum((2*treatments - 1)*weights)

    return(A.hat.aipw)
  }
  
  t.grid <- seq(-1, 1, 0.005)
  values <- sapply(t.grid, value)
  
  t.aipw <- t.grid[which.max(values)]  
  return(t.aipw)
}

N <- 100000

estimate.mean.thresholds <- function(n){
  registerDoParallel(cores = 30)
  
  results <- foreach (i = 1:N, .combine = "rbind") %dopar% {
    data <- generate.data(n)
    learn.policy(data$Ys1, data$Ys0, data$Xs, data$ps, data$indicators)
  } %>% data.frame()
  
  return(colMeans(results))
}

estimate.standard.error <- function(n){
  results <- foreach (i = 1:10, .combine = "rbind") %do% {
    estimate.mean.thresholds(n)
  }
  return(results)
}

for (n in c(250, 500, 750, 1000)){
  print(estimate.mean.thresholds(n))
}

for (n in c(250, 500, 750, 1000)){
  print(n)
  results <- estimate.standard.error(n) 
  errors <- sqrt(apply(results, 2, var))
  cat(paste(n, " ", errors[[1]], ",", errors[[2]], "\n"),
      file = "policy_results.txt",
      append = T)
}



