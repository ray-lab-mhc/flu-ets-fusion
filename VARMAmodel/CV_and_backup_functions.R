
#' Cross-valid function uses a series of training and test sets to calculate the forecast accuracy
#' measures(RMSE) on each test set (of each time series) and the results are averaged across all test sets (of each series).
#' This function is only available for multivariate time series using VAR,VARMA, and VMA model. Only first-order difference is available 
#' in this function.
#' @param: orig_ts original multivariate time series (must be ts object)
#' @param: window the minimum size of training set
#' @param: val_p value of order p
#' @param: val_q: value of order q
#' @param: freq. If the original time series needs differencing, enter the lag of interest. If no, enter 0
#' @param: lamb. If the original time series needs transforming, enter the value of lambda. If no, enter 37.

cross_valid <- function(orig_ts, window, h, val_p, val_q,freq, lamb){
  #create the indices for time series data splitting
  time_slice <- createTimeSlices(orig_ts[,1], initialWindow = window, horizon = h, fixedWindow = FALSE)
  train <- time_slice$train
  test <- time_slice$test
  train_lst <- list()
  test_lst  <- list()
  pred_list <- list()
  avg_RMSE <- list()
  #create a series of training sets based on the indices above
  for (i in 1: length(train)){
    train_lst[[i]] <- orig_ts[range(train[i])[1]: range(train[i])[2], 1:ncol(orig_ts)]
  }
  
  #transform training sets with a given lambda
  if (lamb != 37){
    transformed_ts <- list()
    for (i in 1: length(train_lst)){
      transformed_ts[[i]] <- BoxCox(train_lst[[i]], lamb)
    }
  }
  #fit a desired model on each training set and predict h-step forecasts  
  for ( i in 1: length(train_lst)){
    
    pred_list[[i]] <- model(train_lst[[i]], h = h, val_p = val_p,val_q = val_q, freq = freq, lambda = lamb)
  }
  
  
  ##convert the forecasts into the original scale
  if(freq > 0){
    temp_list <- list()
    if(lamb == 37){  #invert differencing only 
      temp_list <- invert_forecast(pred_list, train_lst, freq = freq, h)
    }
    else{ ## invert differencing first, then Box-Cox transformation
      temp_list <- invert_forecast(pred_list, transformed_ts, freq = freq, h) ## invert forecast back to the transformed scale
      
      for( i in 1: length(temp_list)){
        temp_list[[i]] <- InvBoxCox(temp_list[[i]], lambda = lamb)##invert to the original scale
      }
    }
  }
  if (freq == 0 && lamb != 37){ #inverse Box-Cox transformation only
    pred_list[[i]] <- InvBoxCox(pred_list[[i]], lambda = lamb)
  } 
  
  
  tmp <- list()#####create a temporary list
  rmse <- matrix(nrow = length(test),ncol = ncol(orig_ts))#create a matrix to hold RMSE 
  result <- matrix(0, nrow = h, ncol = ncol(orig_ts))
  for (i in 1: length(test)){
    test_lst[[i]] <- orig_ts[range(test[i])[1]: range(test[i])[2], 1:ncol(orig_ts)]##### list of test sets
    if(freq  > 0){
      if(h ==1 ){
        test_lst[[i]] - temp_list[[i]] -> error #compute the errors(residuals)
        tmp[[i]] <- error
      }
      else{
        as.data.frame(test_lst[[i]]) - as.data.frame(temp_list[[i]]) -> error #####compute the errors(residuals)
        tmp[[i]] <- error
      }
    }
    else{
      if( h == 1){
        test_lst[[i]] - pred_list[[i]] -> error #####compute the errors(residuals)
        tmp[[i]] <- error  #
      }
      else{
        as.data.frame(test_lst[i]) - as.data.frame(pred_list[i]) -> error #####compute the errors for the forecast observations(residuals)
        tmp[[i]] <- error  #
      }
    }  
  } 
  for(j in 1:length(tmp)){
    result <- result + (tmp[[j]])^2  #calculate MSE for each horizon
    
  }
  
  
  RMSE <- sqrt(result/length(test))
  
  return(RMSE)
}




#' invert_forecast function invert first-order differenced forecasts of a series of training  sets
#' @param pred_lst a list containing first-order differenced forecasts of a series of training sets
#' @param train_lst a list containg a series of training data set
#' @param freq  ts_frequency
#' @param h  forecast horizon
#' @return time series at original scales
invert_forecast <- function(pred_lst, train_lst, freq, h){
  lst <- list()
  for ( i in 1: length(pred_lst)){
    init_series <- train_lst[[i]]
    diff <- pred_lst[[1]]
    if (h > 1){
      lst[[i]] <- diff + init_series[(nrow(init_series) - freq + 1): (nrow(init_series) - freq + nrow(diff)), 1: ncol(init_series)]
    }
    else{
      lst[[i]] <- diff + init_series[(nrow(init_series) - freq + 1), 1: ncol(init_series)]
    }
  }
  return(lst)
}

#' invert_fc function invert first-ordered differenced forecasts
#' @param dy first-ordered differenced forecasts
#' @param y original time series
#' @param freq ts- frequency
#' @param h forecast horizon

invert_fc <- function(dy,y,freq,h){
  if (h > 1){
    orig_fc <- dy + y[(nrow(y) - freq + 1): (nrow(y) - freq + nrow(dy)), 1: ncol(y)]
  }
  else{
    orig_fc <- dy + y[(nrow(y) - freq + 1), 1: ncol(y)]
  }
}


#' pred_var function fits an apppropriate model for multivariate time series and return h-step forecasts.
#' @param: y multivariate time series
#' @param: h forecast horizon
#' @param: val_p value of p
#' @param: val_q value of q
#' @return: h-step forecast
pred_var <- function(y, h, val_p, val_q) {
  if (val_q == 0){
    fitVAR <- VAR(y, p = val_p)
    return(VARpred(fitVAR, h = h)$pred)
  }
  else if (val_p ==0){
    fitVMA <- VMA(y, q = val_q)
    return(VMApred(fitVMA, h = h)$pred)
  }
  else{
    fitVARMA <- VARMA(y, p = val_p, q = val_q)
    return(VARMApred(fitVARMA, h = h)$pred)
  }
  
  
}
#' model function performs differencing on time series if a valid value of frequency is passed and calls
#' pred_var function inside to fit an appropriate model for the series and return forecasts.
#' @paramL y multivariate ts
#' @param h forecast horizon
#' @param val_p: value of order p
#' @param val_q: value of order q
#' @param freq ts frequence
#' @param lambda tranformation param
model <- function(y,h, val_p, val_q, freq, lambda){
  y <- BoxCox_trans(y, lambda = lambda)
  
  if(freq != 0){
    y %>% diff(lag = freq) -> y_diff
    pred_var(y_diff,h, val_p, val_q)
  }
  else{
    pred_var(y,h, val_p, val_q)
  }
}

#' BoxCox_trans returns a Box-Cox transformed time series if a valid lambda is passed.
#' @param y multivariate time series.
#' @param lambda transformation parameter.
#' @return transformed y

BoxCox_trans <- function(y, lambda){
  if (lambda != 37){
    y <- BoxCox(y,lambda)
    return (y)
  }
  return(y)
}

