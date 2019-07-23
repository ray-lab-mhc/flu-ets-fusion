#' simulate_incidence function 
#' 
#' @param model_fit
#' @param backcast_steps integer number of steps to "backcast"; simulations of
#'   incidence will be generated for the last backcast_steps time points 
#' @param forecast_steps integer number of steps to forecast
#' @param n_sim integer number of paths to simulate
#' @param new_regression_data a data frame of regression variables only.
#' 
#' @return A horizon x num_ts x n_sim array of sampled incidence.
#' 
#' @export
simulate_incidence <- function(model_fit, backcast_steps, forecast_steps, n_sim, new_regression_data) {
  results <- array(NA, dim = c(backcast_steps + forecast_steps, attr(model_fit, "p"), n_sim))
  
  filter_results <- KFAS::KFS(model_fit)
  
  for(step in rev(seq_len(backcast_steps))) {
    alpha_t <- matrix(filter_results$alphahat[nrow(filter_results$alphahat) - step + 1, ])
    for(sim_ind in seq_len(n_sim)) {
      results[backcast_steps - step + 1, , sim_ind] <- simulate_y_t(model_fit, alpha_t)
    }
  }
  
  if(forecast_steps > 0) {
    # newdata object contains system matrices for future times
    # in particular, newdata_SSMmodel$Z contains covariates which may differ at future times.
    if(!missing(new_regression_data)) {
      regression_data <- new_regression_data
    }
    newdata_matrix <- matrix(NA, nrow = forecast_steps, ncol = attr(model_fit, "p"))
    newdata_SSMmodel <- KFAS::SSModel(as.formula(paste0("newdata_matrix ~ ", as.character(model_fit$terms)[3])))
    
    for(sim_ind in seq_len(n_sim)) {
      alpha_t <- matrix(filter_results$alphahat[nrow(filter_results$alphahat), ])
      for(step in seq_len(forecast_steps)) {
        alpha_t <- simulate_alpha_t(model_fit, alpha_t)
        results[backcast_steps + step, , sim_ind] <- simulate_y_t(model_fit, alpha_t)
      }
    }
  }
  
  return(results)
}

#' generate a single simulated value for y_t based on a mode fit and state vector alpha_t
#' 
#' @param model_fit a SSMmodel fit object
#' @param alpha_t a vector (or matrix with one column) with states at time t
#' 
#' @return simulated values for the observed quantities y_t
simulate_y_t <- function(model_fit, alpha_t) {
  epsilon_t <- t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nrow(model_fit$H)), sigma = as.matrix(model_fit$H[, , 1]), method = "chol"))
  Z <- model_fit$Z[,,1]
  if(is.null(dim(Z))) {
    Z <- matrix(Z, ncol = nrow(alpha_t))
  }
  y_t <- Z %*% alpha_t + epsilon_t
  return(y_t)
}

#' generate a single simulated value for alpha_t based on a mode fit and state vector alpha_{t-1}
#' 
#' @param model_fit a SSMmodel fit object
#' @param alpha_tm1 a vector (or matrix with one column) with states at time t - 1
#' 
#' @return simulated values for the state at time t
simulate_alpha_t <- function(model_fit, alpha_tm1) {
  eta_t <- t(mvtnorm::rmvnorm(n = 1, mean = rep(0, nrow(model_fit$Q)), sigma = as.matrix(model_fit$Q[, , 1]), method = "chol"))
  T <- model_fit$T[, , 1]
  if(is.null(dim(T))) {
    T <- matrix(T, ncol = nrow(alpha_tm1))
  }
  alpha_t <- T %*% alpha_tm1 + model_fit$R[, , 1] %*% eta_t
}
