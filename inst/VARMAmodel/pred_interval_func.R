
## CI_function calculated 95% prediction interval.
#@param: pred h-step forecasts
#@param: standard error of each values
#@param: h forecast horizon
#@param: ts_frequency for differenced time series (enter 0 if the series do not need experiencing)
#@param: orig_series original series
#@param: lambda transformation param (enter 37 if no transformation is needed.)
CI_function <- function(pred, se, h, ts_frequency, orig_series, lambda){
  
  if (h >1){
    upper_bound <- matrix(nrow = nrow(pred), ncol = ncol(pred))
    lower_bound <- matrix(nrow = nrow(pred), ncol = ncol(pred))
    for(i in 1: ncol(pred)){
      for (j in 1: nrow(pred)){
        upper_bound[j,i] = pred[j,i] + 1.96*se[j,i] 
      }
    }
    for(i in 1: ncol(pred)){
      for (j in 1: nrow(pred)){
        lower_bound[j,i] = pred[j,i] - 1.96*se[j,i] 
      }
    }
  }
  else{
    upper_bound <- matrix(nrow = nrow(se), ncol = ncol(orig_series))
    lower_bound <- matrix(nrow = nrow(se), ncol = ncol(orig_series))
    for ( i in 1:ncol(orig_series)){
      
      upper_bound[i] = pred[i] + 1.96*se[i] 
      lower_bound[i] = pred[i] - 1.96*se[i]
      
    }
  }
  
  if(ts_frequency > 0){
    invert_fc(upper_bound, orig_series, freq = ts_frequency, h = h)-> upper_bound
    invert_fc(lower_bound, orig_series, freq = ts_frequency, h = h)-> lower_bound
  }
  if(lambda != 37){
    InvBoxCox(upper_bound, lambda = lambda) -> upper_bound
    InvBoxCox(lower_bound, lambda = lambda) -> lower_bound
  }
  
  
  
  print("Here are upper limits of 95% prediction interval")
  print(upper_bound)
  print("Here are lower limits of 95% prediction interval")
  print(lower_bound)
  
}





## function un-doing differencing in multivariate ts
undo_diff_m <- function(diff, orig, ts_frequency){
  fc_lst <- list()
  for (i in 1:ncol(orig)){
    
    fc_each <- invert_seasonal_difference(diff[,i], orig[,i], ts_frequency = ts_frequency)
    fc_lst[[i]] <- as.data.frame(fc_each)
    
  }
  return(fc_lst)
}