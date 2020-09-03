

sim_forecast <- function(ts_data,
                         unfit_model,
                         covariate = TRUE,
                         horizon, 
                         lambda,
                         regression_variables,
                         regression_data,
                         init_value_H,
                         trend_degree,
                         seasonal_period,
                         seasonal_type = "common",
                         cycle_num,
                         cycle_base_period,
                         cycle_type = "common",
                         return_disturbance_matrix = FALSE,
                         nsim){
  #obtain forecast values and convert back to original scales 
  forecast_f <- SSModel_prediction(ts_data,
                                       unfit_model,
                                       covariate = TRUE,
                                       horizon, 
                                       lambda,
                                       regression_variables,
                                       regression_data,
                                       init_value_H,
                                       trend_degree,
                                       seasonal_period,
                                       seasonal_type = "common",
                                       cycle_num,
                                       cycle_base_period,
                                       cycle_type = "common",
                                       return_disturbance_matrix = FALSE)
  forecast_value <- forecast_f$predicted_val
  traj_data <- simulate_incidence(forecast_value$model_fit,backcast_steps = 0, forecast_steps = horizon, n_sim = nsim) 
  traj_data <- InvBoxCox(traj_data,lambda = lambda)
  return(list(forecast_value = forecast_value,traj_data = traj_data))
  
}


#horizon <- rep(i, length(bin_counts))
#location <- rep(location,length(bin_counts))
#dist_df <- cbind.data.frame(binlwr_df,horizon,location)
#print(dist_df[1:131,])


log_dist <- function(traj_val,location,epi_week){
  # bin endpoint definitions
  lwr <- seq(from = 0.0, to = 13.0, by = 0.1)
  bins <- c(lwr, 100)
  
  for(i in 1:nrow(traj_val)){
    
    x <- traj_val[i,1,]
    file_name <- sprintf("logdist_data/%s_%d_%s.csv",location,i,epi_week)
    print(file_name[1])
    # get probability of each bin based on sample values
    
    # get count and probability for each bin
    bin_counts <- hist(x, bins, right = FALSE, plot = FALSE)$counts 
    binlwr_df <- data.frame(bin_start_incl = lwr,
                            bin_end_notincl = bins[2:length(bins)],
                            value = bin_counts / sum(bin_counts),
                            target = paste0(i, " wk ahead"),
                            type = "Bin",
                            unit = "percent",
                            location = location)
    write.csv(binlwr_df, file = file_name[1],append = FALSE, row.names = FALSE)
    
  }
}






log_dist(sc_sim_inv_o,"SC")
log_dist(ca_sim_inv_o, "CA")











































    
write.csv(final_df,'logdist_data/SC_1.csv')
# example -- x would be your forecast samples at a particular horizon (like h = 1)
x <- c(0.05, 12.87, 13.2)

log_dist(sc_sim_inv_o,"SC")

x <- sc_sim_inv_o[2,1,]
# get probability of each bin based on sample values
# bin endpoint definitions
lwr <- seq(from = 0.0, to = 13.0, by = 0.1)
bins <- c(lwr, 100)

# get count and probability for each bin
bin_counts <- hist(x, bins, right = FALSE, plot = FALSE)$counts 
binlwr_df <- data.frame(lwr = lwr,
                        prob = bin_counts / sum(bin_counts) )
h <- rep("1", length(bin_counts))
l <- rep("SC", length(bin_counts))
final_df <- cbind.data.frame(binlwr_df,h,l)





sct <- SSModel_prediction(ts_data = sc_quidel_ts_o,
                               unfit_model = sc_t_o,
                               covariate = TRUE,
                               horizon = 4, 
                               lambda = 0,
                               init_value_H = c(mse1_o,mse2_o),
                               trend_degree = 1,
                               
                               cycle_num = 2,
                               cycle_base_period = 52,
                               cycle_type = "common",
                               return_disturbance_matrix = FALSE)

