#' Create an SSMmodel object and auxiliary information for a call to optim
#' 
#' @param ts_data matrix of time series values, one column per time series and
#'   one row per time point.  NA's are allowed.
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
#' @param return_optim_args boolean; if TRUE, also return functions and initial
#'   parameter values that could be used in a call to optim for parameter
#'   estimation.
#' 
#' @return named list.  If \code{return_optim_args} is FALSE, contains only
#'   \code{model}, an SSMmodel object.  If \code{return_optim_args} is TRUE,
#'   also contains: \code{update_pars}, a function to update the \code{model}
#'   parameters; \code{calc_loglik}, a function to calculate the \code{model}
#'   log likelihood; and \code{init_pars}, a vector of initial parameter values.
#' 
#' @export
create_SSMmodel <- function(
  ts_data,
  regression_variables,
  regression_data,
  init_value_H,
  init_value_Q,
  trend_degree,
  seasonal_period,
  seasonal_type = "common",
  cycle_num,
  cycle_base_period,
  cycle_type = "common",
  return_optim_args = FALSE) {
  p <- ncol(ts_data)
  
  if (sum(is.na(ts_data)) > 0){
    print("true")
    start_date <- start(ts_data)
    ts_frequency <- frequency((ts_data))
    ts_data <- as.data.frame(ts_data)
    ts_data <- ts_data %>%
      filter_all(all_vars(!is.na(.)))
    ts_data <- ts(ts_data, start = start_date, frequency = ts_frequency)
    print(ts_data)
  } 
  
  
  
  if(missing(init_value_Q)){
    init_pars <- NULL
  }else{
    init_pars <- init_value_Q
  }
  
  
  if(!missing(regression_variables)) {
    if(missing(regression_data)) {
      stop("If regression_variables are provided, regression_data must be too.")
    }
    eval(regression_data)
    reg_formula <- paste0(
      "SSMregression(~ ", paste0(regression_variables, collapse = " + "), ", Q =matrix(NA), ",
      "data = regression_data, ",
      "index = 1:", p, ")")
    if(return_optim_args) {
      if(missing(init_value_Q)){
        for (i in 1: length(regression_variables)){
          regression_init_pars <- 0.001
          init_pars <- c(init_pars, regression_init_pars)
        }
      }else{
        init_pars <- c(init_pars) 
      }
      
    }
  }else {
    reg_formula <- NULL
  }
  
  if(!missing(trend_degree)) {
    trend_formula <- paste0(
      "SSMtrend(", trend_degree,
      ", Q = rep(list(matrix(NA, nrow = ", p,
      ", ncol = ", p, ")), ", trend_degree, "))")
    
    if(return_optim_args) {
      if(missing(init_value_Q)){
        ts_chol <- chol(cov(ts_data) / 10)
        trend_init_pars <- c(log(diag(ts_chol)), ts_chol[upper.tri(ts_chol)])
        init_pars <- c(init_pars, trend_init_pars)
      } else{
        init_pars <- c(init_pars)
      }
    }
  } else {
    trend_formula <- NULL
  }
  
  if(!missing(seasonal_period)) {
    seasonal_formula <- paste0(
      "SSMseasonal(period = ", seasonal_period,
      ", Q = NA, sea.type = \"trigonometric\", type = \"", seasonal_type, "\")"
    )
    
    if(return_optim_args) {
      if(missing(init_value_Q)){
        seasonal_init_pars <- log(0.1)
        init_pars <- c(init_pars, seasonal_init_pars)
      }else{
        init_pars <- c(init_pars)
      }
    }
  } else {
    seasonal_formula <- NULL
  }
  
  if(!missing(cycle_num)) {
    if(missing(cycle_base_period)) {
      stop("If cycle_num is provided, cycle_base_period must be too.")
    }
    cycle_formula <- ""
    for(i in seq_len(cycle_num)) {
      cycle_formula <- paste0(
        cycle_formula,
        if(i > 1) {" + "} else {NULL},
        "SSMcycle(period = ", cycle_base_period / i,
        ", Q = matrix(NA), type = \"", cycle_type, "\")"
      )
      
      if(return_optim_args) {
        if(missing(init_value_Q)){
          cycle_init_pars <- 0.01
          init_pars <- c(init_pars, cycle_init_pars)
        }else{
          init_pars <- c(init_pars)
        }
      }
      
    }
  } else {
    cycle_formula <- NULL
  }
  
  ssmodel_formula <- paste0(
    "ts_data ~ ",
    paste(c(reg_formula, trend_formula, seasonal_formula, cycle_formula),
          collapse = " + ")
  )
  ssmodel_formula <- as.formula(ssmodel_formula)
  #init_pars <- c(init_pars, log(diag(cov(ts_data))/100)) 
  if(!missing(init_value_H)){
    init_pars <- c(init_pars, log(init_value_H)) 
  }else{
    init_pars <- c(init_pars)
  }
  
  model <- KFAS::SSModel(ssmodel_formula, H = diag(NA_real_, nrow = ncol(ts_data), ncol = ncol(ts_data)))
  
  if(return_optim_args) {
    update_pars <- function(pars, model, init_value_H) {
      cur_par_ind <- 1
      cur_q_ind <- 1
      
      
      if(!is.null(reg_formula)){
        
        model$Q[1,1,1] <- exp(pars[1]) ##updating Q in regression
        
        cur_par_ind <- cur_par_ind + 1
        cur_q_ind <- cur_q_ind + 1
      }
      
      
      if(!is.null(trend_formula)) {
        Q <- diag(exp(pars[seq(from = cur_par_ind, length = p)]), nrow = p, ncol = p)
        
        cur_par_ind <- cur_par_ind + p
        
        Q[upper.tri(Q)] <- pars[seq(from = cur_par_ind, length = choose(p, 2))]
        
        cur_par_ind <- cur_par_ind + choose(p, 2)
        
        q_inds <- seq(from = cur_q_ind, len = p)
        model$Q[q_inds, q_inds, 1] <- crossprod(Q)
        cur_q_ind <- cur_q_ind + p
        
      }
      
      if(!is.null(seasonal_formula)) {
        q_inds <- seq(from = cur_q_ind,
                      length = sum(attr(model, "eta_types") == "seasonal"))
        #model$Q[q_inds, q_inds, 1] <- diag(exp(pars[cur_par_ind]), nrow = p, ncol = p)
        diag(model$Q[q_inds, q_inds, 1]) <- exp(pars[cur_par_ind])
        cur_par_ind <- cur_par_ind + 1
        cur_q_ind <- cur_q_ind + length(q_inds)
      }
      
      if(!is.null(cycle_formula)) {
        for(i in seq_len(cycle_num)) {
          q_inds <- seq(from = cur_q_ind, length = 2)
          
          diag(model$Q[q_inds, q_inds, 1]) <- exp(pars[cur_par_ind])
          
          
          
          cur_par_ind <- cur_par_ind + 1
          cur_q_ind <- cur_q_ind + 2
        }
      }
      if(!missing(init_value_H)){
        diag(model$H[,,1]) <- exp(tail(pars, ncol(ts_data)))
      }else{
        diag(model$H[,,1]) <- 1
      }
      
      return(model)
    }
    calc_loglik <- function(pars, model) {
      model <- update_pars(pars, model)
      return(logLik(model))
    }
  }
  
  if(return_optim_args) {
    return(list(
      model = model,
      update_pars = update_pars,
      calc_loglik = calc_loglik,
      init_pars = init_pars))
  } else {
    return(list(model = model))
  }
}

