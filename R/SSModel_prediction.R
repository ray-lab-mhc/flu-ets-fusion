
#'Predict the future observations of a state space model of class SSMmodel.
#'
#' @param ts_data matrix of time series values, one column per time series and
#'   one row per time point.  NA's are allowed.
#'   @param unfit_model a state space object of class SSModel
#' @param regression_variables character vector of variables to use in a
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
#' @return a matrix containing the predictions
#'  
#' @export


SSModel_prediction <- function(ts_data,
                               unfit_model, 
                               horizon, 
                               regression_variables,
                               regression_data,
                               init_value_H,
                               trend_degree,
                               seasonal_period,
                               seasonal_type = "common",
                               cycle_num,
                               cycle_base_period,
                               cycle_type = "common",
                               return_disturbance_matrix = FALSE){
  
    
  
  
  unfit_model$model$P1inf[] <- 0  ##removes the diffuse initialization
  diag(unfit_model$model$P1) <- 100  ##prior variance of the coefficient
  
  model_fit_pars <- optim(
    fn = unfit_model$calc_loglik,
    par = unfit_model$init_pars,
    model = unfit_model$model,
    control = list(fnscale = -1)
  )
  if(!missing(init_value_H)){
    model_fit <- unfit_model$update_pars(pars = model_fit_pars$par, model = unfit_model$model,init_value_H = init_value_H)
  }else{
    model_fit <- unfit_model$update_pars(pars = model_fit_pars$par, model = unfit_model$model)
  }
  print(model_fit$H)
  print(model_fit$Q)
  
  if(!missing(regression_variables) && !missing(regression_data)){  
    combined_data <- cbind.data.frame(ts_data, regression_data)
    reg_dat <- as.data.frame(combined_data[1:horizon,]) 
  }else{
    reg_dat <- NULL 
  }
  forecast_val <- ts(matrix(NA, nrow = horizon, ncol = ncol(ts_data)), frequency = frequency(ts_data))
    
  new_data <- create_SSMmodel_pred(ts_data, 
                                  forecast_val,
                                  model_fit, 
                                  horizon, 
                                  regression_variables, 
                                  reg_dat, 
                                  trend_degree,
                                  seasonal_period,
                                  seasonal_type = "common",
                                  cycle_num,
                                  cycle_base_period,
                                  cycle_type = "common")
  
  predicted_val <- predict(model_fit, new_data, n.ahead = horizon)
  
  if(return_disturbance_matrix){
    return(predicted_val)
  }else{
    return(list(
      predicted_val = predicted_val,
      H = model_fit$H,
      Q = model_fit$Q,
      model_fit = model_fit
    ))
  }
}


#'Create a compatible SSModel object to be added in the end of the old object for which the predictions are required.
#'
#' @param ts_data matrix of time series values, one column per time series and
#'   one row per time point.  NA's are allowed.
#' @param old_model a state space object of class SSModel for which the predictions are required
#' @param horizon number of steps ahead at which to predict
#' @param regression_variables character vector of variables to use in a
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
#' @return a SSModel object
#'  
#' @export

create_SSMmodel_pred <- function(ts_data,
                                 pred_matrix,
                                 old_model,
                                 horizon,
                                 regression_variables,
                                 regression_data,
                                 trend_degree,
                                 seasonal_period,
                                 seasonal_type = "common",
                                 cycle_num,
                                 cycle_base_period,
                                 cycle_type = "common"){
  p <- ncol(ts_data)
  cur_par_ind <- 1

  
  if(!missing(regression_variables) && !is.null(regression_data)) {
    reg_formula <- paste0(
      "SSMregression(~ ", paste0(regression_variables, collapse = " + "), ", Q = old_model$Q[1,1,1], ",
      "data = regression_data", ")" )
    cur_par_ind_l <- cur_par_ind + 1
    cur_par_ind <- cur_par_ind + p
    print(reg_formula)
  } 
  else {
    reg_formula <- NULL
    cur_par_ind_l <- 1 
    cur_par_ind <- p
  }
  
 
  if(!missing(trend_degree)) {
    # par_inds  <- seq_along(from = cur_par_ind, length = p)
    # cur_par_ind <- cur_par_ind + p
    
    trend_formula <- paste0(
      "SSMtrend(", trend_degree,
      ", Q = old_model$Q[",cur_par_ind_l, ": ",cur_par_ind, ",", cur_par_ind_l, ":", cur_par_ind, ",1])")
    print(trend_formula)
    cur_par_ind <- cur_par_ind + 1
  } else {
    trend_formula <- NULL
    cur_par_ind <- cur_par_ind_l
  }
  
  if(!missing(seasonal_period)) {
    seasonal_formula <- paste0(
      "SSMseasonal(period = ", seasonal_period,
      ", Q = old_model$Q[", cur_par_ind, ",", cur_par_ind, ",1],", "sea.type = \"trigonometric\", type = \"", seasonal_type, "\")"
    )
    
    
  } else {
    seasonal_formula <- NULL
  }
  
  if(!missing(cycle_num)) {
   
    cycle_formula <- ""
    for(i in seq_len(cycle_num)) {
      cycle_formula <- paste0(
        cycle_formula,
        if(i > 1) {" + "} else {NULL},
        "SSMcycle(period = ", cycle_base_period / i,
        ", Q =  old_model$Q[", cur_par_ind,",", cur_par_ind, ",1],", "type = \"", cycle_type, "\")"
      )
      
      cur_par_ind <- cur_par_ind + 2
    
  }} else {
    cycle_formula <- NULL
  }
  
  ssmodel_formula <- paste0(
    "pred_matrix ~ ",
    paste(c(reg_formula, trend_formula, seasonal_formula, cycle_formula),
          collapse = " + ")
  )
  ssmodel_formula <- as.formula(ssmodel_formula)
  print(ssmodel_formula)
  model <- KFAS::SSModel(ssmodel_formula, H = old_model$H[,,1])
  return(model)
  
  
}
   

  
 


