#' Cross-valid function uses a series of training and test sets to calculate the forecast accuracy
#' measures(RMSE) on each test set (of each time series) and the results are averaged across all test sets (of each series).
#' This CV procedure is available for multivariate time series using exponential smoothing 
#' models.

#' @param: orig_ts original multivariate time series (must be ts object)
#' @param: missing- TRUE if there is one observation of wiki data than of non-wiki data
#' @param: is_dummy
#' @param: dummy - dummy variables(e.g: holiday effects) in the form of ts object. 0 if no dummny variables needed
#' @param: window the minimum size of training set
#' @param: horizon the size of each test set
#' @param: level. 1 for local level model and 0 for only seasonal component
#' @param: seasonal.  
#' @param: freq. the frequency of time series
#' @param: lamb. If the original time series needs transforming, enter the value of lambda. If no, enter 37.
#' @inits: inits initial values for covariance structures
#' @update_func: updating function to update the model given the parameters



CV_es <- function(orig_ts,missing, is_dummy, dummy, window, horizon, seasonal,cycle, freq, lamb,inits, update_func){
  #create the indices for time series data splitting
  time_slice <- createTimeSlices(orig_ts[,1], initialWindow = window, horizon = horizon, fixedWindow = FALSE)
  train <- time_slice$train
  test <- time_slice$test
  train_lst <- list()
  test_lst  <- list()
  pred_list <- list()
  avg_RMSE <- list()
  dummy_lst <- list()
  #create a series of training sets based on the indices above
  for (i in 1: length(train)){
    train_lst[[i]] <- orig_ts[range(train[i])[1]: range(train[i])[2], 1:ncol(orig_ts)]
    if(is_dummy == "yes"){
      
      dummy_lst[[i]] <- dummy[range(train[i])[1]: range(train[i])[2],]
      
    }
    if(missing == TRUE){
      train_lst[[i]][nrow(train_lst[[i]]), 1:2] <- NA #omit the last observations of ili data
    }
  }
  
  
  #fit a desired model on each training set and predict h-step forecasts  
  for ( i in 1: length(train_lst)){
    
    pred_list[[i]] <- model(y = train_lst[[i]], h = horizon, seasonal= seasonal,cycle = cycle,is_dummy = is_dummy, dummy = dummy_lst[[i]], inits = inits,freq = freq, update_func = update_func, lambda = lamb)
  }
  
  
  ##convert the forecasts into the original scale
  
  if(lamb != 37){  
    fc_trans <- undo_BoxCox(horizon = horizon, pred = pred_list, lambda = lamb)
  }
  
  
  
  tmp <- list() #create a temporary list
  #rmse <- matrix(nrow = length(test),ncol = ncol(orig_ts))#create a matrix to hold RMSE 
  result <- matrix(0, nrow = horizon, ncol = ncol(orig_ts))
  for (i in 1: length(test)){
    test_lst[[i]] <- orig_ts[range(test[i])[1]: range(test[i])[2], 1:ncol(orig_ts)] # list of test sets
    if(lamb != 37){ #for transformed ts
      test_lst[[i]] - fc_trans[[i]] -> error # compute the errors(residuals)
      tmp[[i]] <- error
    } else{
      if( horizon == 1){
        test_lst[[i]] - as.data.frame(pred_list[[i]]) -> error # compute the errors(residuals)
        tmp[[i]] <- error  
      } else{
        as.data.frame(test_lst[i]) - as.data.frame(pred_list[[i]]) -> error # compute the errors for the forecast observations(residuals)
        tmp[[i]] <- error  
      }
    }  
  } 
  
  for(j in 1:length(tmp)){
    result <- result + (tmp[[j]])^2  #calculate MSE for each horizon
      
    }
    
  
  RMSE <- sqrt(result/length(test))
  
  return(RMSE)
}







# pred_mod function fits a desired model for multivariate time series and return h-step forecasts.
#' @param: y multivariate time series
#' @param: h forecast horizon
#' @param: level. 1 for local level model and 0 for only seasonal component
#' @param: seasonal.  
#' @param: freq. the frequency of time series
#' @param: lamb. If the original time series needs transforming, enter the value of lambda. If no, enter 37.
#' @param: inits initial values for covariance structures
#' @param: update_func updating function to update the model given the parameters
#' @return: return h-step forecasts

pred_mod <- function(y, horizon, seasonal,cycle, inits,freq, update_func) {
  
  # local level model
  if(seasonal == "trig"){
    fit_model <- SSModel(y ~ SSMtrend(1, Q = list(matrix(NA)))+ SSMseasonal(period = freq, Q = NA, sea.type = "trigonometric", type = "common"), data = as.data.frame(y))
    print(fit_model)
  }
  else if (seasonal == "cycle"){
    if(cycle == 1){
      fit_model <- SSModel(y ~ SSMtrend(1, Q = list(matrix(NA)))+ SSMcycle(period = freq, Q = NA,type = "common"), data = as.data.frame(y))
      
    }
    else if( cycle == 2){
      fit_model <- SSModel(y~SSMtrend(1, Q = list(matrix(NA))) + SSMcycle(period = freq, 
                                                                          Q = matrix(NA), type = "common")+SSMcycle(period = freq/2, Q = matrix(NA),type = "common"), data = as.data.frame(y))
      fit_model$P1inf[] <- 0  ##removes the diffuse initialization
      diag(fit_model$P1) <- 100
    }
    else if(cycle == 3){
      fit_model <- SSModel(y~SSMtrend(1, Q = list(matrix(NA))) + SSMcycle(period = freq, 
                                                                          Q = matrix(NA), type = "common")+SSMcycle(period = freq/2, Q = matrix(NA),type = "common")+ SSMcycle(period = freq/4, Q  = matrix(NA), type = "common") , data = as.data.frame(y))
      fit_model$P1inf[] <- 0  ##removes the diffuse initialization
      diag(fit_model$P1) <- 100
    }
  }
  
  
  
  
  fit_opt <- optim(p = inits, f = update_func,
                   method = "L-BFGS-B", model = fit_model)
  model_optim <- update_func(fit_opt$par,fit_model, estimate = FALSE)
  pred_test <- predict(model_optim, n.ahead = horizon)
  return (pred_test)
}


#' model function performs transfomation if valid value of lambda is passed and  calls
#' pred_mod function inside to fit an appropriate model for the series and return forecasts.
#' @param:  y multivariate ts
#' @param   h forecast horizon
#' @param: level. 1 for local level model and 0 for only seasonal component
#' @param: seasonal  
#' @param: cycle
#' @param: is_dummy
#' @param: dummy
#' @param: freq. the frequency of time series
#' @param: lamb. If the original time series needs transforming, enter the value of lambda. If no, enter 37.
#' @param: inits initial values for covariance structures
#' @param: update_func updating function to update the model given the parameters
model <- function(y, h, seasonal,cycle,is_dummy, dummy, inits,freq, update_func, lambda){
  if(lambda != 37){
    y <- BoxCox(y, lambda = lambda)
  }
  if(is_dummy == "no"){
    pred_mod(y, h, seasonal = seasonal,cycle = cycle, inits = inits,freq = freq, update_func = update_func) 
  }else{
    dummy_mod(y,h,seasonal = seasonal, cycle = cycle, dummy = dummy, inits = inits, freq = freq, update_func = update_func)
  }
}


#' undo BoxCox function converts forecasting values into the original scale.
#' @param: horizon forecast horizon
#' @param: pred forecasting values
#' @param: lambda 
#' 
#' @return forecasts at original scales
undo_BoxCox <- function(horizon, pred, lambda){
  fc <- data.frame(matrix(nrow = horizon, length(pred)))
  
  for (i in 1:length(pred)){
    fc[[i]] <- InvBoxCox(as.data.frame(pred[[i]]), lambda = lambda) 
  }
  return(fc)
}

#' dummy_mod function fits a desired model with dummy variables for multivariate time series . 
#' Note: only available for one dummy variable. You may need to write your own dummy_mod function if more than one dummy variable are used.
#' @param: y multivariate time series
#' @param: dummy - dummy variables
#' @param: h forecast horizon
#' @param: level. 1 for local level model and 0 for only seasonal component
#' @param: seasonal.  
#' @param: freq. the frequency of time series
#' @param: lamb. If the original time series needs transforming, enter the value of lambda. If no, enter 37.
#' @param: inits initial values for covariance structures
#' @param: update_func updating function to update the model given the parameters
#' @return: return h-step forecasts




dummy_mod <- function(y,dummy,horizon, seasonal,cycle, inits,freq, update_func) {
  d <- cbind.data.frame(y ,dummy)
  fc <- ts(matrix(NA, nrow = horizon, ncol = ncol(y)), frequency = freq)
  if(seasonal == "trig"){
    fit_model <- SSModel(y ~ SSMregression(~ dummy, Q = matrix(NA), index = 1:ncol(y)) + SSMtrend(1, Q = list(matrix(NA)))+ 
                           SSMseasonal(period = freq, Q = NA, sea.type = "trigonometric", type = "common"))
    fit_model$P1inf[] <- 0  ##removes the diffuse initialization
    diag(fit_model$P1) <- 100
    
    fit_opt <- optim(p = inits, f = update_func,
                     method = "L-BFGS-B", model = fit_model)
    model_optim <- update_func(fit_opt$par,fit_model, estimate = FALSE)
    
    new_data <- SSModel(fc ~ SSMregression(~ dummy, Q = model_optim$Q[1,1,1], data = as.data.frame(d[1:horizon,])) + SSMtrend(1, Q = list(model_optim$Q[2:6, 2:6,1]))+ 
                          SSMseasonal(period = freq, Q = model_optim$Q[7, 7, 1], sea.type = "trigonometric", type = "common"))
    
  }
  else if(seasonal == "cycle"){
    if (cycle == 1){
      fit_model <- SSModel(y ~ SSMregression(~ dummy, Q = matrix(NA), index = 1:ncol(y))+ 
                             SSMtrend(1, Q = list(matrix(NA)))+ SSMcycle(period = freq, Q = NA,type = "common"))
      fit_model$P1inf[] <- 0  ##removes the diffuse initialization
      diag(fit_model$P1) <- 100
      fit_opt <- optim(p = inits, f = update_func,
                       method = "L-BFGS-B", model = fit_model)
      model_optim <- update_func(fit_opt$par,fit_model, estimate = FALSE)
      
      new_data <- SSModel(fc ~ SSMregression(~ dummy, Q = model_optim$Q[1,1,1], data = as.data.frame(d[1:horizon,]))
                          + SSMtrend(1, Q = list(model_optim$Q[2:6, 2:6,1]))
                          + SSMcycle(period = freq, Q = model_optim$Q[7,7,1], type = "common"))
      print(new_data)
      
    }
    else if (cycle == 2){
      fit_model <- SSModel(y ~ SSMregression(~ dummy, Q = matrix(NA), index = 1:ncol(y))+ 
                             SSMtrend(1, Q = list(matrix(NA)))+ SSMcycle(period = freq, Q = NA,type = "common") + SSMcycle(period = freq/2, Q = NA, type = "common"))
      
      fit_model$P1inf[] <- 0  ##removes the diffuse initialization
      diag(fit_model$P1) <- 100
      
      fit_opt <- optim(p = inits, f = update_func,
                       method = "L-BFGS-B", model = fit_model)
      model_optim <- update_func(fit_opt$par,fit_model, estimate = FALSE)
      
      new_data <- SSModel(fc ~ SSMregression(~ dummy, Q = model_optim$Q[1,1,1], data = as.data.frame(d[1:horizon,]))
                          + SSMtrend(1, Q = list(model_optim$Q[2:6, 2:6,1]))
                          + SSMcycle(period = freq, Q = model_optim$Q[7,7,1], type = "common") + SSMcycle(period = freq/2, Q = model_optim$Q[9,9,1], type = "common"))
      
    }
    else if (cycle == 3){
      fit_model <- SSModel(y ~ SSMregression( ~ dummy, Q = matrix(NA), index = 1:ncol(y))+ 
                             SSMtrend(1, Q = list(matrix(NA)))+ SSMcycle(period = freq, Q = NA,type = "common") + SSMcycle(period = freq/2, Q = NA, type = "common")+
                             SSMcycle(period = freq/4, Q =NA, type = "common"))
      
      fit_opt <- optim(p = inits, f = update_func,
                       method = "L-BFGS-B", model = fit_model)
      model_optim <- update_func(fit_opt$par,fit_model, estimate = FALSE)
      
      new_data <- SSModel(fc ~ SSMregression(~dummy, Q = model_optim$Q[1,1,1], data = as.data.frame(d[1:horizon,]))
                          + SSMtrend(1, Q = list(model_optim$Q[2:6, 2:6,1]))
                          + SSMcycle(period = freq, Q = model_optim$Q[7,7,1], type = "common") + SSMcycle(period = freq/2, Q = model_optim$Q[9,9,1], type = "common")+
                            SSMcycle(period = freq/4, Q = model_optim$Q[11,11,1], type = "common"))
      
    }
  }
  
  pred_test <- predict(model_optim, newdata = new_data)
  return (pred_test)
}


