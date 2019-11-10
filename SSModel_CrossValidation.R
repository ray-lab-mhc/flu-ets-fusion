#'Compute cross-validation statistics (RMSE) using a rolling time origin procedure.
#'
#' @param ts_data matrix of time series values, one column per time series and
#'   one row per time point.  NA's are allowed.
#' @param extra_wiki whether there are extra wiki data in the data set
#' @param initial_train_size the initial number of values in each training set sample
#' @param test_size the number of values in test set sample
#' @param regression_var character vector of variables to use in a
#'   regression component of the model.  If missing, no regression component is
#'   included.
#' @param regression_data data frame with variables for regression
#' @param trend_degree degree of trend component; if missing, no trend
#'   component is included.
#' @param seasonal_period period of seasonality term in mdoel; if missing, no
#'   seasonal component is included.
#' @param seasonal_type whether seasonal states are shared across time series
#'   ("common") or distinct for each time series ("distinct")
#' @param cycle_num number of cycle components to include.  The \code{k}'th
#'   cycle component will have period \code{cycle_base_period / k}.
#' @param cycle_base_period base period for cycle components.
#' @param cycle_type whether cycle states are shared across time series
#'   ("common") or distinct for each time series ("distinct")
#' @return a matrix containing RMSE for h-steps ahead
#'  
#' @export

SSmodel_CV <- function( ts_data,
                        extra_wiki = TRUE, 
                        initial_train_size,
                        test_size,
                        lambda_val,
                        init_value_H,
                        init_value_Q,
                        regression_var,
                        regression_data,
                        trend_degree,
                        seasonal_period,
                        seasonal_type = "common",
                        cycle_num,
                        cycle_base_period,
                        cycle_type = "common"){
                        
                    
  #create the indices for time series data splitting
  time_slice <- caret::createTimeSlices(ts_data[,1], initialWindow = initial_train_size, horizon = test_size, fixedWindow = FALSE)
  train <- time_slice$train
  test <- time_slice$test
  train_set <- list()
  test_set  <- list()
  forecast_lst <- list()
  avg_RMSE <- list()
  regression_df <- list()
  if(!missing(lambda_val)){
    ts_data_orig <- ts_data
    ts_data <- forecast::BoxCox(ts_data, lambda = lambda_val)
    
  }
  
  #create a series of training sets based on the indices above
  for (i in 1: length(train)){
    
    train_set[[i]] <- ts_data[range(train[i])[1]: range(train[i])[2], 1:ncol(ts_data)]
    
    if(!missing(lambda_val)){
      
      test_set[[i]] <-  ts_data_orig[range(test[i])[1]: range(test[i])[2], 1:ncol(ts_data_orig)] # list of test sets
      
    }else{
      test_set[[i]] <-  ts_data[range(test[i])[1]: range(test[i])[2], 1:ncol(ts_data)] # list of test sets
    }
    
    if(!missing(regression_var)){
      
      regression_df[[i]] <- regression_data[range(train[i])[1]: range(train[i])[2],1:ncol(regression_data)]
      
    }
    
    
    
    if(extra_wiki){
      train_set[[i]][nrow(train_set[[i]]), 1:2] <- NA #omit the last observations of ili data
      
    }
  }
  
  
  #fit a desired model on each training set and predict h-step forecasts  
  for ( i in 1: length(train_set)){
    
    unfit_mod <- create_SSMmodel(train_set[[i]],
                                 regression_variables = regression_var,
                                 regression_data = regression_df[[i]],
                                 init_value_H = init_value_H,
                                 init_value_Q = init_value_Q,
                                 trend_degree = trend_degree,
                                 seasonal_period = seasonal_period,
                                 seasonal_type = "common",
                                 cycle_num = cycle_num,
                                 cycle_base_period = cycle_base_period,
                                 cycle_type = "common",
                                 return_optim_args = TRUE)
    
    
    forecast_val <- SSModel_prediction(ts_data = train_set[[i]],
                                       unfit_mod, 
                                       horizon = test_size, 
                                       regression_variables = regression_var,
                                       regression_data = regression_df[[i]],
                                       init_value_H = init_value_H,
                                       trend_degree = trend_degree,
                                       seasonal_period = seasonal_period,
                                       seasonal_type = "common",
                                       cycle_num = cycle_num,
                                       cycle_base_period = cycle_base_period,
                                       cycle_type = "common") 
    
    if(!missing(lambda_val)){
      forecast_val <- forecast::InvBoxCox(as.data.frame(forecast_val), lambda = lambda_val)
    }
    
    forecast_lst[[i]] <- forecast_val
  
  }
  
  error_lst <- list() 
  
  result <- matrix(0, nrow = test_size, ncol = ncol(ts_data))
  
  for (i in 1: length(test)){
    
    if( test_size == 1){
      
      error <- test_set[[i]] - as.data.frame(forecast_lst[[i]])  #compute (y - y_hat)
      
      error_lst[[i]] <- error
    }else{
      error <- as.data.frame(test_set[i]) - forecast_lst[[i]]  #compute (y - y_hat)
      
      error_lst[[i]] <- error 
    }
    
  } 
  
  
  for(j in 1:length(error_lst)){
    
    result <- result + (error_lst[[j]])^2  #calculate MSE for each horizon
    
  }
  
  
  RMSE <- sqrt(result/length(test))
  
  return(RMSE)
  
}


