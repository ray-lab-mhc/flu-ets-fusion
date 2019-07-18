
#' simulate_incidence function 
#' @param model_fit
#' @param horizon
#' @param n number of paths to simulate n
#' @param dummy a data frame(or ts) of only regression variables
#' @return A (number of time series) by h by n array of sampled incidence.


simulate_incidence <- function(model_fit,horizon, n, num_ts, dummy){
  sampled_inc <- array(NA, dim = c(n, horizon, num_ts)) #create an array to contain simulated series
  sampled_inc2 <- array(NA, dim = c(n, horizon, num_ts))
  sampled_inc3 <- array(NA, dim = c(n, horizon, num_ts))
  xmas_week <- when_is_xmas(dummy)
  out <- KFS(model_fit)
  for ( i in 1: n){
    alpha_t <- unname(out$a[nrow(out$a),]) 
    
    for(j in 1:horizon) {
      update_results <- update_function(model_fit, alpha_t, xmas_week)
      alpha_t <- update_results[[1]]
      sampled_inc[i, j, ] <- update_results[[2]]
      sampled_inc2[i,j, ]  <- update_results[[3]]
      sampled_inc3[i,j, ]  <- update_results[[4]]
    }
  }
  return(list(sampled_inc,sampled_inc2, sampled_inc3))
}

#' update_function updates alpha_t and y_t (by randomizing eta_t and epsilon_t). 
#' @param model_fit
#' @param alpha_t
#' @param xmas_week 
#' @return a list of updated alpha_t and y_t for h-step forwards

update_function <- function(model_fit, alpha_t,xmas_week){
  eta_t <- t(rmvnorm(n = 1, mean = rep(0, nrow(model_fit$Q)), sigma = model_fit$Q[, , 1], method = "chol"))
  
  epsilon_t <- t(rmvnorm(n = 1, mean = rep(0, nrow(model_fit$H)), sigma = model_fit$H[, , 1], method = "chol"))
  
  alpha_t <- (model_fit$T[, , 1])%*%alpha_t + (model_fit$R[, , 1]) %*% eta_t
  
  y_t <- model_fit$Z[,,1] %*% alpha_t + epsilon_t
  
  y_t_prime <- model_fit$Z[,,unlist(sample(xmas_week,1))] %*% alpha_t + epsilon_t #calculate y_t_prime from Z matrix which has xmas effect
  
  y_t_prime2 <- model_fit$Z[,,65]  %*% alpha_t + epsilon_t# obtain Z matrix at t = 65 (christmas_week)
  
  
  return (list(alpha_t, y_t, y_t_prime, y_t_prime2))
}


#' when_is_xmas function returns a list of christmas week from a given data frame 
#' @param xmas_cov data frame(or ts) of regression models(e.g xmas)
#' @return xmas_week 
when_is_xmas<- function(xmas_cov){
  xmas_week <- list()
  for ( i in 1: length(xmas_cov)){
    if(xmas_cov[[i]] == TRUE){
      
      xmas_week[[i]] <- i
    }
  }
  
  xmas_week <- xmas_week[-which(sapply(xmas_week, is.null))]
  return(xmas_week)
}


# Example: simulate_incidence(model_optim_d5, 4,4,5,flu_data2$christmas_week)
# You can use the model below to test simulate_incidence
# model for testing:
#d_model_5 <- SSModel(flu_ts2[,1:5]~ SSMregression(~christmas_week+ postchristmas_week, Q = matrix(NA), data = flu_ts2, index = 1:5) + SSMtrend(1, Q = list(matrix(NA)))
# + SSMcycle(period = 52, Q = matrix(NA), type = "common"), data = flu_ts2)

#uplik_model5 <- function(pars, model, estimate = TRUE){
# model$Q[1,1,1] <- exp(pars[1])##updating Q in regression
#  Q <- diag(exp(pars[2:6]))
#  Q[upper.tri(Q)] <- pars[7:16]
#  model$Q[2:6, 2:6, 1] <- crossprod(Q)##updating Q in trend
#  diag(model$Q[7:8,7:8,1]) <- exp(pars[17])##updating Q in cycle
#  if(estimate){
#    -logLik(model)
#  }
#  else{
#    model
#  }
#}

#fit_opt_d5 <- optim(par = c(inits_d2[1:17]), fn = uplik_model5, method = "BFGS", model = d_model_5)##use optim to optimize initial values for the params over an update function( which calculate the likelihoo)
#model_optim_d5 <- uplik_model5(fit_opt_d5$par,d_model_5, estimate = FALSE)##update the model given the optimized params.



