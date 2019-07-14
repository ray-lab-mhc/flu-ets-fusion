## A collection of updating functions that are required for the estimation of parameters.

## 1) Local level + Cycle(52)

## model_lc <- SSModel(flu_ts~SSMtrend(1, Q = list(matrix(NA))) + SSMcycle(period = 52, Q = matrix(NA), type = "common"), data = flu_ts)

uplik_lc <- function(pars, model, estimate = TRUE){
  Q <- diag(exp(pars[1:5]))
  Q[upper.tri(Q)] <- pars[6:15]
  model$Q[1:5, 1:5, 1] <- crossprod(Q)
  diag(model$Q[6:7, 6:7, 1]) <- exp(pars[16])
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
  
}

## 2) Local level + Cycle(52 and 26)

##model_lcc <- SSModel(flu_ts~SSMtrend(1, Q = list(matrix(NA))) + SSMcycle(period = 52, 
                                                                         Q = matrix(NA), type = "common")+SSMcycle(period = 26, Q = matrix(NA),type = "common"), data = flu_ts)
##model_lcc$P1inf[] <- 0  ##removes the diffuse initialization
##diag(model_lcc$P1) <- 100  ##prior variance of the coefficient
uplik_model_lcc <- function(pars, model, estimate = TRUE){
  Q <- diag(exp(pars[1:5]))
  Q[upper.tri(Q)] <- pars[6:15]
  model$Q[1:5, 1:5, 1] <- crossprod(Q)##update Q in SSMtrend
  diag(model$Q[6:7, 6:7, 1]) <- exp(pars[16])##update Q in SSMcycle period 52
  diag(model$Q[8:9, 8:9, 1]) <- exp(pars[17])##update Q in SSMcycle period 26
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
  
}

## 3) Local level + Cycle(52 and 26 and 13)


##model_lc3 <- SSModel(flu_ts~SSMtrend(1, Q = list(matrix(NA))) + SSMcycle(period = 52, 
                                                                         Q = matrix(NA), type = "common")+SSMcycle(period = 26, Q = matrix(NA),type = "common") + SSMcycle(period = 13, Q = matrix(NA), type = "common") , data = flu_ts)
##model_lc3$P1inf[] <- 0  ##removes the diffuse initialization
##diag(model_lc3$P1) <- 100  ##prior variance of the coefficient

uplik_model_lc3 <- function(pars, model, estimate = TRUE){
  Q <- diag(exp(pars[1:5]))
  Q[upper.tri(Q)] <- pars[6:15]
  model$Q[1:5, 1:5, 1] <- crossprod(Q)##update Q in SSMtrend
  diag(model$Q[6:7, 6:7, 1]) <- exp(pars[16])##update Q in SSMcycle period 52
  diag(model$Q[8:9, 8:9, 1]) <- exp(pars[17])##update Q in SSMcycle period 26
  diag(model$Q[10:11, 10:11, 1]) <- exp(pars[18])##update Q in SSMcycle period 26
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
  
}

## 4) Model with dummy variables(christmas only): SScycle(period = 52)+ Local level model

##d_model_5 <- SSModel(flu_ts3[,1:5]~ SSMregression(~christmas_week+ postchristmas_week+thanksgiving_week + postthanksgiving_week, Q = matrix(NA), data = flu_ts3, index = 1:5) + SSMtrend(1, Q = list(matrix(NA)))
## + SSMcycle(period = 52, Q = matrix(NA), type = "common"), data = flu_ts3)

uplik_lcc_d_mod <- function(pars, model, estimate = TRUE){
  model$Q[1,1,1] <- exp(pars[1])
  Q <- diag(exp(pars[2:6]))
  Q[upper.tri(Q)] <- pars[7:16]
  model$Q[2:6, 2:6, 1] <- crossprod(Q)
  diag(model$Q[7:8, 7:8, 1]) <- exp(pars[17])
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
}
## 5) Model with dummy variables(christmas only): SScycle(period = 52 and period = 26)+ Local level model
##m1 <- SSModel(flu_ts2[, 1:5] ~ SSMregression(~flu_ts2[,6:7], Q = matrix(NA), index = 1:5) + SSMtrend(1, Q = list(matrix(NA))) + 
                SSMcycle(period = 52, Q = matrix(NA), type = "common")+ SSMcycle(period = 52, Q = matrix(NA), type = "common"), data = flu_ts2)
##m1$P1inf[] <- 0  ##removes the diffuse initialization
##diag(m1$P1) <- 100

uplik_m1 <- function(pars, model, estimate = TRUE){
  model$Q[1,1,1] <- exp(pars[1])
  Q <- diag(exp(pars[2:6]))
  Q[upper.tri(Q)] <- pars[7:16]
  model$Q[2:6, 2:6, 1] <- crossprod(Q)
  diag(model$Q[7:8, 7:8, 1]) <- exp(pars[17])
  diag(model$Q[9:10, 9:10,1])<- exp(pars[18])
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
}

## 6) Model with dummy variables(christmas only): SScycle(period = 52 and period = 26 and period = 13)+ Local level model

##m2 <- SSModel(flu_ts2[, 1:5] ~ SSMregression(~flu_ts2[,6:7], Q = matrix(NA), index = 1:5) + SSMtrend(1, Q = list(matrix(NA))) + 
                SSMcycle(period = 52, Q = matrix(NA), type = "common")+ SSMcycle(period = 26, Q = matrix(NA), type = "common")+
                SSMcycle(period = 13, Q = matrix(NA), type = "common"), data = flu_ts2)
##m1$P1inf[] <- 0  ##removes the diffuse initialization
##diag(m1$P1) <- 100

uplik_m2 <- function(pars, model, estimate = TRUE){
  model$Q[1,1,1] <- exp(pars[1])
  Q <- diag(exp(pars[2:6]))
  Q[upper.tri(Q)] <- pars[7:16]
  model$Q[2:6, 2:6, 1] <- crossprod(Q)
  diag(model$Q[7:8, 7:8, 1]) <- exp(pars[17])
  diag(model$Q[9:10, 9:10,1])<- exp(pars[18])
  diag(model$Q[11:12,11:12,1]) <- exp(pars[19])
  if(estimate){
    -logLik(model)
  }
  else{
    model
  }
}
