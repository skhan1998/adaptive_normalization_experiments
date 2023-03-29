hajek.estimator <- function(Ys, ps, indicators){
  n <- length(Ys)
  S.hat <- sum(indicators*Ys/ps)
  n.hat <- sum(indicators/ps)
  
  return(S.hat/n.hat)
}

hajek.var.estimator <- function(Ys, ps, indicators){
  mu.hat <- hajek.estimator(Ys, ps, indicators)
  qs = (1-ps)/ps
  n.hat <- sum(indicators/ps)
  return(sum((Ys - mu.hat)^2*indicators/ps)/n.hat + sum(qs*(Ys-mu.hat)^2*indicators/ps)/(n.hat))
}

ht.estimator <- function(Ys, ps, indicators){
  n <- length(Ys)
  S.hat <- sum(indicators*Ys/ps)
  return(S.hat/n)
}

ht.var.estimator <- function(Ys, ps, indicators){
  qs = (1-ps)/ps
  n.hat <- sum(indicators/ps)
  mu.hat <- ht.estimator(Ys, ps, indicators)
  
  return(sum((Ys - mu.hat)^2*indicators/ps)/n.hat + sum(qs*(Ys)^2*indicators/ps)/(n.hat))
}

lambda.oracle.estimator <- function(Ys, ps, indicators){
  n <- length(Ys)
  S.hat <- sum(indicators*Ys/ps)
  n.hat <- sum(indicators/ps)
  
  lambda.star <- -cov(Ys, qs)/(mean(Ys)*mean(qs))
  mu.hat.lambda.star <- S.hat/(lambda.star*n + (1-lambda.star)*n.hat)
  
  return(mu.hat.lambda.star)
}

lambda.iterative.estimator <- function(Ys, ps, indicators, iter){
  n <- length(Ys)
  S.hat <- sum(indicators*Ys/ps)
  n.hat <- sum(indicators/ps)
  qs = (1-ps)/ps
  pi.hat <- sum(qs*indicators/ps)/n.hat
  T.hat <- sum(Ys*qs*indicators/ps)/n.hat
  mu.hat.0 <- S.hat/n
  rho.hat <- sum((Ys-mu.hat.0)*(qs)*indicators/ps)/n.hat
  
  mu.hats <- rep(0, iter)
  mu.hats[1] <- mu.hat.0
  
  idx = 2
  while (idx <= iter){
    mu.hat <- mu.hats[idx-1]
    
    A <- S.hat/n
    B <- (T.hat/pi.hat)*(1-(n.hat/n))
    
    mu.hat.next <- A/(1-(B/mu.hat))
    
    # rho.hat <- sum((Ys-mu.hat)*(qs)*indicators/ps)/n.hat
    # lambda.star.hat <- -rho.hat/(mu.hat*pi.hat)
    # 
    # mu.hat.next <- S.hat/(lambda.star.hat*n + (1-lambda.star.hat)*n.hat)
    mu.hats[idx] <- mu.hat.next
    idx <- idx + 1
  }
  
  return(mu.hats)
}

fp.estimator <- function(Ys, ps, indicators){
  n <- length(Ys)
  S.hat <- sum(indicators*Ys/ps)
  n.hat <- sum(indicators/ps)
  qs = (1-ps)/ps

  T.hat <- sum(Ys*qs*indicators/ps)/n.hat
  pi.hat <- sum(qs*indicators/ps)/n.hat

  mu.hat.fixed <- (S.hat/n) + (T.hat/pi.hat)*(1-n.hat/n)
  return(mu.hat.fixed)
}

fp.var.estimator <- function(Ys, ps, indicators){
  n <- length(Ys)
  
  S.hat <- sum(indicators*Ys/ps)
  n.hat <- sum(indicators/ps)
  qs = (1-ps)/ps

  T.hat <- sum(Ys*qs*indicators/ps)/n.hat
  pi.hat <- sum(qs*indicators/ps)/n.hat

  mu.hat <- fp.estimator(Ys, ps, indicators)
  
  
  return(sum((Ys - mu.hat)^2*indicators/ps)/n.hat + sum(qs*(Ys-T.hat/pi.hat)^2*indicators/ps)/(n.hat))
}

aipw.estimator <- function(Ys, ps, Xs, indicators){
  n = length(Ys)
  fold = sample(1:n, n/2, replace = F)
  train = rep(0, n)
  train[fold] = 1
  
  fold.1 <- aipw.estimator.fold(Ys, ps, Xs, train, indicators)
  fold.2 <- aipw.estimator.fold(Ys, ps, Xs, 1-train, indicators)
  
  mu.hat.1 <- mean(fold.1$Ys.hat) + sum(fold.1$resid.hat)/fold.1$n.hat
  mu.hat.2 <- mean(fold.2$Ys.hat) + sum(fold.2$resid.hat)/fold.2$n.hat
  
  return(0.5*(mu.hat.1 + mu.hat.2))
}

aipw.var.estimator <- function(Ys, ps, Xs, indicators){
  n = length(Ys)
  
  mu.hat <- fp.estimator(Ys, ps, indicators)
  n.hat <- sum(indicators/ps)
  
  sigma.y.hat <- sum((Ys - mu.hat)^2*indicators/ps)/n.hat
  rho <- sum(Ys*Xs*indicators/ps)/n.hat - mu.hat*mean(Xs)
  
  if (rho > 1){
    rho <- 1
  }
  
  if (rho < -1){
    rho <- -1
  }
  
  #ps.hat <- glm(indicators ~ Xs, family = "binomial")$fitted.values
  
  var.hat <- rho^2*sigma.y.hat + sigma.y.hat*(1-rho)*mean(1/ps)

  return(var.hat)
}


aipw.estimator.fold <- function(Ys, ps, Xs, train, indicators){
  n <- length(Ys)
  Ys.observed <- Ys[indicators*train == 1]
  ps.observed <- ps[indicators*train == 1]
  Xs.observed <- Xs[indicators*train == 1]
  indicators.observed <- indicators[train == 1]

  outcome.model.df <- data.frame(Xs.observed, Ys.observed)
  colnames(outcome.model.df) <- c("Xs", "Ys")
  
  #outcome.model <- lm(Ys ~ Xs, data = outcome.model.df)
  outcome.model <- gam(Ys ~ s(Xs, k=5), data = outcome.model.df)
  #outcome.model <- randomForest(Ys ~ Xs, data = model.df, surface = "direct")
  
  propensity.model.df <- data.frame(Xs[train==1], indicators[train==1])
  colnames(propensity.model.df) <- c("Xs", "indicators")
  
  propensity.model <- glm(indicators ~ Xs, family = "binomial", data = propensity.model.df)
  
  #print(summary(propensity.model))
  
  predict.df <- data.frame(Xs[train==0])
  colnames(predict.df) = c("Xs")
  Ys.hat <- predict(outcome.model, newdata = predict.df)
  ps.hat <- predict(propensity.model, type = "response", newdata = predict.df)
  
  n.hat <- sum(indicators[train==0]/ps.hat)
  #n.hat <- 0.5*length(Ys)
  
  # mu.hat.reg <- mean(Ys.hat) + sum((Ys[train==0] - Ys.hat)*indicators[train==0]/ps.hat)/n.hat
  
  scores <- Ys.hat + (Ys[train==0]-Ys.hat)*indicators[train==0]/ps.hat

  return(list("scores"=scores, 
              "Ys.hat" = Ys.hat, 
              "resid.hat" = (Ys[train==0]-Ys.hat)*indicators[train==0]/ps.hat,
              "n.hat" = n.hat))
}

aipw.fp.estimator <- function(Ys, ps, Xs, indicators){
  n = length(Ys)
  fold = sample(1:n, n/2, replace = F)
  train = rep(0, n)
  train[fold] = 1

  mu.hat.1 <- aipw.fp.estimator.fold(Ys, ps, Xs, train, indicators)
  mu.hat.2 <- aipw.fp.estimator.fold(Ys, ps, Xs, 1-train, indicators)

  return(0.5*(mu.hat.1 + mu.hat.2))
}

aipw.fp.estimator.fold <- function(Ys, ps, Xs, train, indicators){
  n = length(Ys)
  Ys.observed <- Ys[indicators*train == 1]
  ps.observed <- ps[indicators*train == 1]
  Xs.observed <- Xs[indicators*train == 1]
  indicators.observed <- indicators[train == 1]
  
  outcome.model.df <- data.frame(Xs.observed, Ys.observed)
  colnames(outcome.model.df) <- c("Xs", "Ys")
  
  #outcome.model <- lm(Ys ~ Xs, data = outcome.model.df)
  outcome.model <- gam(Ys ~ s(Xs, k=5), data = outcome.model.df)
  #model <- loess(Ys ~ ps, data = model.df, surface = "direct")
  
  propensity.model.df <- data.frame(Xs[train==1], indicators[train==1])
  colnames(propensity.model.df) <- c("Xs", "indicators")
  
  propensity.model <- glm(indicators ~ Xs, family = "binomial", data = propensity.model.df)
  
  #print(summary(propensity.model))
  
  predict.df <- data.frame(Xs[train==0])
  colnames(predict.df) = c("Xs")
  Ys.hat <- predict(outcome.model, newdata = predict.df)
  ps.hat <- predict(propensity.model, type = "response", newdata = predict.df)
  
  #mu.hat.reg <- mean(Ys.hat) + sum((Ys[train==0] - Ys.hat)*indicators[train==0]/ps[train==0])/n.hat
  
  n.hat <- sum(indicators[train==0]/ps.hat)
  
  qs = (1-ps.hat)/ps.hat
  pi.hat <- sum(qs*indicators[train==0]/ps.hat)

  Zs <- Ys[train==0] - Ys.hat
  TZ.hat <- sum(Zs*qs*indicators[train==0]/ps.hat)
  
  correction <- (TZ.hat/pi.hat)*(1-2*n.hat/n)
  
  mu.hat <- mean(Ys.hat) + sum(Zs*indicators[train==0]/ps.hat)/(0.5*length(Ys))
  mu.hat <- mu.hat + correction
  return(mu.hat)
}



double.fp.estimator <- function(Ys1, Ys0, ps, indicators){
  n <- length(Ys1)
  
  n1 <- sum(indicators/ps)
  n0 <- sum((1-indicators)/(1-ps))
  
  qs1 <- (1-ps)/ps
  qs0 <- (ps)/(1-ps)
  
  pi1 <- mean(qs1*indicators/ps)
  pi0 <- mean(qs0*(1-indicators)/(1-ps))
  
  S1 <- sum(Ys1*indicators/ps)
  S0 <- sum(Ys0*(1-indicators)/(1-ps))
  
  T1 <- mean(Ys1*qs1*indicators/ps)
  T0 <- mean(Ys0*qs0*(1-indicators)/(1-ps))
  
  C <- pi0*pi1-1
  
  A <- matrix(c(1 + (1-n1/n)/(C),
                -(pi1*(1-n0/n))/(C),
                -(pi0*(1-n1/n))/(C),
                1 + (1-n0/n)/(C)), nrow = 2)
  
  b <- c(S1/n + (T1*pi0 - T0)*(1-n1/n)/(C),
         S0/n + (T0*pi1 - T1)*(1-n0/n)/(C))
  
  mu.hat <- solve(A) %*% b
  return(mu.hat[1,1] - mu.hat[2,1])
}


color_map <- c("hajek" = rgb(31, 119, 180, maxColorValue = 256),
               "ht" = rgb(255, 127, 14, maxColorValue = 256),
               "fp" = rgb(44, 160, 44, maxColorValue = 256),
               "aipw" = rgb(214, 39, 40, maxColorValue = 256),
               "aipw.fp" = rgb(148, 103, 189, maxColorValue = 256),
               "double.fp" = rgb(140, 86, 75, maxColorValue = 256),               
               "super.fp" = rgb(227, 119, 194, maxColorValue = 256))

name_map <- c("hajek" = unname(TeX('$\\hat{\\mu}_{Hajek}$')),
              "ht"= unname(TeX('$\\hat{\\mu}_{HT}$')),
              "fp" = unname(TeX('$\\hat{\\mu}_{AN}$')),
              "aipw" = unname(TeX('$\\hat{\\mu}_{AIPW}$')),
              "aipw.fp" = unname(TeX('$\\hat{\\mu}_{AIPW, AN}$')),
              "super.fp" = unname(TeX('$\\hat{\\mu}_{Super, AN}$')))

name_map_tau <- c("hajek" = unname(TeX('$\\hat{\\tau}_{Hajek}$')),
              "ht"= unname(TeX('$\\hat{\\tau}_{HT}$')),
              "fp" = unname(TeX('$\\hat{\\tau}_{AN}$')),
              "aipw" = unname(TeX('$\\hat{\\tau}_{AIPW}$')),
              "aipw.fp" = unname(TeX('$\\hat{\\tau}_{AIPW, AN}$')),
              "double.fp" = unname(TeX('$\\hat{\\tau}_{AN^2}$')),
              "lin" = unname(TeX('$\\hat{\\tau}_{Interact}$')))



