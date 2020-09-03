# CA

# July 2018 (until week 26 of 2018)
ca <- c(states[5], quidel_avail[5])

ca_quidel <- training_data[1:404, ca]
ca_quidel[1][404,] <- NA
ca_test <- training_data[405:408,ca]
ca_quidel_transform <- BoxCox(ca_quidel + 0.05, 0)
ca_quidel_ts <- ts(ca_quidel_transform, start = c(2018,01), frequency = 52)
# model with identity H matrix
ca_nonH <- create_SSMmodel(ca_quidel_ts,
                           covariate = TRUE,
                           trend_degree = 1,
                           cycle_num = 2,
                           cycle_base_period = 52,
                           cycle_type = "common",
                           return_optim_args = TRUE)

ca_nonH$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_nonH$model$P1) <- 100  ##prior variance of the coefficient

ca_mod <- optim(
  fn = ca_nonH$calc_loglik,
  par = ca_nonH$init_pars,
  model =ca_nonH$model,
  control = list(fnscale = -1)
)
ca_fit <- ca_nonH$update_pars(pars = ca_mod$par, model = ca_nonH$model)
ca_pred_4 <- predict(ca_fit, n.ahead = 4)
ca_pred_inv <- InvBoxCox(as.data.frame(ca_pred_4), lambda = 0) - 0.05

#get MSE
ca_sse <- (ca_pred_inv - ca_test)^2
mse1 <- mean(unlist(ca_sse[1])) 
mse2 <- mean(unlist(ca_sse[2]))
#update H matrix
ca_t <- create_SSMmodel(ca_quidel_ts,
                        covariate = TRUE,
                        init_value_H = c(mse1,mse2),
                        
                        trend_degree = 1,
                        
                        cycle_num = 2,
                        cycle_base_period = 52,
                        cycle_type = "common",
                        return_optim_args = TRUE)

ca_t$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_t$model$P1) <- 100

ca_opt <- optim(
  fn = ca_t$calc_loglik,
  par = ca_t$init_pars,
  model = ca_t$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



ca_fit_2 <- ca_t$update_pars(pars = ca_opt$par, model = ca_t$model,init_value_H = c(mse1,mse2))
ca_pred_2 <- predict(ca_fit_2, n.ahead = 4)
ca_pred_inv_2 <- InvBoxCox(as.data.frame(ca_pred_2), lambda = 0) - 0.05
ca_pred_ts_updated_h <- ts(ca_pred_inv_2[1], start = c(2018,27), frequency=52)
ca_pred_ts_nonH <- ts(ca_pred_inv[1], start = c(2018,27), frequency=52)

autoplot(ts(ca_test[1], start = c(2018,27), frequency = 52), ylab = "forecasted ILI")+
  autolayer(ca_pred_ts_updated_h)+
  autolayer(ca_pred_ts_nonH)+
  xlab("horizon")
  ylab("forecasted ILI")
 
ca_simulate_incidence <- simulate_incidence(ca_fit_2,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
ca_sim_inv <- InvBoxCox(ca_simulate_incidence,lambda = 0)
ca_getCI <- quantile(ca_sim_inv[1,1,],.025)
ca_getCI_u <-quantile(ca_sim_inv[1,1,],.975) 

ca_getCI_2 <- quantile(ca_sim_inv[2,1,],.025)
ca_getCI_u2 <-quantile(ca_sim_inv[2,1,],.975) 

ca_getCI_3 <- quantile(ca_sim_inv[3,1,],.025)
ca_getCI_u3 <-quantile(ca_sim_inv[3,1,],.975)

ca_getCI_4 <- quantile(ca_sim_inv[4,1,],.025)
ca_getCI_u4 <-quantile(ca_sim_inv[4,1,],.975) 

ca_pred_df <- as.data.frame(ca_pred_ts_updated_h)
horizon <- c(1,2,3,4)
ca_pred_df <- cbind.data.frame(ca_pred_df, horizon,ca_test[1])
ggplot()+
  geom_line(mapping = aes( x = ca_pred_df$horizon, y = ca_pred_df$fit, color = "Forecasted ILI"))+
  geom_line(mapping = aes( x = ca_pred_df$horizon, y = ca_pred_df$wili_ca, color = "Observed ILI"))+
  geom_ribbon(aes(ymin = c(ca_getCI,ca_getCI_2,ca_getCI_3,ca_getCI_4), 
                  ymax=c(ca_getCI_u, ca_getCI_u2,ca_getCI_u3, ca_getCI_u4),
                  x=ca_pred_df$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("horizon")+
  ylab("ILI")


# Octorber 2018

ca <- c(states[5], quidel_avail[5])

ca_quidel_o <- full_data[1:417, ca]
ca_quidel_o[1][418,] <- NA
ca_test_o <- full_data[418:421,ca]
ca_quidel_transform_o <- BoxCox(ca_quidel_o + 0.05, 0)
ca_quidel_ts_o <- ts(ca_quidel_transform_o, start = c(2018,01), frequency = 52)
# model with identity H matrix
ca_nonH_o <- create_SSMmodel(ca_quidel_ts_o,
                           covariate = TRUE,
                           trend_degree = 1,
                           cycle_num = 2,
                           cycle_base_period = 52,
                           cycle_type = "common",
                           return_optim_args = TRUE)

ca_nonH_o$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_nonH_o$model$P1) <- 100  ##prior variance of the coefficient

ca_mod_o <- optim(
  fn = ca_nonH_o$calc_loglik,
  par = ca_nonH_o$init_pars,
  model =ca_nonH_o$model,
  control = list(fnscale = -1)
)
ca_fit_o <- ca_nonH_o$update_pars(pars = ca_mod_o$par, model = ca_nonH_o$model)
ca_pred_4_o <- predict(ca_fit_o, n.ahead = 4)
ca_pred_inv_o <- InvBoxCox(as.data.frame(ca_pred_4_o), lambda = 0) - 0.05

#get MSE
ca_sse_o <- (ca_pred_inv_o - ca_test_o)^2
mse1_o <- mean(unlist(ca_sse_o[1])) 
mse2_o <- mean(unlist(ca_sse_o[2]))
#update H matrix
ca_t_o <- create_SSMmodel(ca_quidel_ts_o,
                        covariate = TRUE,
                        init_value_H = c(mse1_o,mse2_o),
                        
                        trend_degree = 1,
                        
                        cycle_num = 2,
                        cycle_base_period = 52,
                        cycle_type = "common",
                        return_optim_args = TRUE)

ca_t_o$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_t_o$model$P1) <- 100

ca_opt_o <- optim(
  fn = ca_t_o$calc_loglik,
  par = ca_t_o$init_pars,
  model = ca_t_o$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



ca_fit_2_o <- ca_t_o$update_pars(pars = ca_opt_o$par, model = ca_t_o$model,init_value_H = c(mse1_o,mse2_o))
ca_pred_2_o <- predict(ca_fit_2_o, n.ahead = 4)
ca_pred_inv_2_o <- InvBoxCox(as.data.frame(ca_pred_2_o), lambda = 0) - 0.05
ca_pred_ts_updated_h_o <- ts(ca_pred_inv_2_o[1], start = c(2018,27), frequency=52)
ca_pred_ts_nonH_o <- ts(ca_pred_inv_o[1], start = c(2018,27), frequency=52)

autoplot(ts(ca_test[1], start = c(2018,27), frequency = 52))+
  autolayer(ca_pred_ts_updated_h)+
  autolayer(ca_pred_ts_nonH)

ca_simulate_incidence_o <- simulate_incidence(ca_fit_2_o,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
ca_sim_inv_o <- InvBoxCox(ca_simulate_incidence_o,lambda = 0)
ca_getCI_o <- quantile(ca_sim_inv_o[1,1,],.025)
ca_getCI_u_o <-quantile(ca_sim_inv_o[1,1,],.975) 

ca_getCI_2_o <- quantile(ca_sim_inv_o[2,1,],.025)
ca_getCI_u2_o <-quantile(ca_sim_inv_o[2,1,],.975) 

ca_getCI_3_o <- quantile(ca_sim_inv_o[3,1,],.025)
ca_getCI_u3_o <-quantile(ca_sim_inv_o[3,1,],.975)

ca_getCI_4_o <- quantile(ca_sim_inv_o[4,1,],.025)
ca_getCI_u4_o <-quantile(ca_sim_inv_o[4,1,],.975) 

ca_pred_df_o <- as.data.frame(ca_pred_ts_updated_h_o)
horizon <- c(1,2,3,4)
ca_pred_df_o <- cbind.data.frame(ca_pred_df_o, horizon,ca_test_o[1])
ggplot()+
  geom_line(mapping = aes( x = ca_pred_df_o$horizon, y = ca_pred_df_o$fit,color = "Observed data"))+
  geom_line(mapping = aes( x = ca_pred_df_o$horizon, y = ca_pred_df_o$wili_ca, color = "Forecasted data"))+
  geom_ribbon(aes(ymin = c(ca_getCI_o,ca_getCI_2_o,ca_getCI_3_o,ca_getCI_4_o), 
                  ymax=c(ca_getCI_u_o, ca_getCI_u2_o,ca_getCI_u3_o, ca_getCI_u4_o),
                  x=ca_pred_df_o$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")

# January 2019 (until final week of 2018)

ca <- c(states[5], quidel_avail[5])

ca_quidel_j <- full_data[1:430, ca]
ca_quidel_j[1][430,] <- NA
ca_test_j <- full_data[431:434,ca]
ca_quidel_transform_j <- BoxCox(ca_quidel_j + 0.05, 0)
ca_quidel_ts_j <- ts(ca_quidel_transform_j, start = c(2018,01), frequency = 52)
# model with identity H matrix
ca_nonH_j <- create_SSMmodel(ca_quidel_ts_j,
                             covariate = TRUE,
                             trend_degree = 1,
                             cycle_num = 2,
                             cycle_base_period = 52,
                             cycle_type = "common",
                             return_optim_args = TRUE)

ca_nonH_j$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_nonH_j$model$P1) <- 100  ##prior variance of the coefficient

ca_mod_j <- optim(
  fn = ca_nonH_j$calc_loglik,
  par = ca_nonH_j$init_pars,
  model =ca_nonH_j$model,
  control = list(fnscale = -1)
)
ca_fit_j <- ca_nonH_j$update_pars(pars = ca_mod_j$par, model = ca_nonH_j$model)
ca_pred_4_j <- predict(ca_fit_j, n.ahead = 4)
ca_pred_inv_j <- InvBoxCox(as.data.frame(ca_pred_4_j), lambda = 0) - 0.05

#get MSE
ca_sse_j <- (ca_pred_inv_j - ca_test_j)^2
mse1_j <- mean(unlist(ca_sse_j[1])) 
mse2_j <- mean(unlist(ca_sse_j[2]))
#update H matrix
ca_t_j <- create_SSMmodel(ca_quidel_ts_j,
                          covariate = TRUE,
                          init_value_H = c(mse1_j,mse2_j),
                          
                          trend_degree = 1,
                          
                          cycle_num = 2,
                          cycle_base_period = 52,
                          cycle_type = "common",
                          return_optim_args = TRUE)

ca_t_j$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_t_j$model$P1) <- 100

ca_opt_j <- optim(
  fn = ca_t_j$calc_loglik,
  par = ca_t_j$init_pars,
  model = ca_t_j$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



ca_fit_2_j <- ca_t_j$update_pars(pars = ca_opt_j$par, model = ca_t_j$model,init_value_H = c(mse1_j,mse2_j))
ca_pred_2_j <- predict(ca_fit_2_j, n.ahead = 4)
ca_pred_inv_2_j <- InvBoxCox(as.data.frame(ca_pred_2_j), lambda = 0) - 0.05
ca_pred_ts_updated_h_j <- ts(ca_pred_inv_2_j[1], start = c(2018,27), frequency=52)
ca_pred_ts_nonH_j <- ts(ca_pred_inv_j[1], start = c(2018,27), frequency=52)

autoplot(ts(ca_test_j[1], start = c(2018,27), frequency = 52))+
  autolayer(ca_pred_ts_updated_h_j)+
  autolayer(ca_pred_ts_nonH_j)

ca_simulate_incidence_j <- simulate_incidence(ca_fit_2_j,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
ca_sim_inv_j <- InvBoxCox(ca_simulate_incidence_j,lambda = 0)
ca_getCI_j <- quantile(ca_sim_inv_j[1,1,],.025)
ca_getCI_u_j <-quantile(ca_sim_inv_j[1,1,],.975) 

ca_getCI_2_j <- quantile(ca_sim_inv_j[2,1,],.025)
ca_getCI_u2_j <-quantile(ca_sim_inv_j[2,1,],.975) 

ca_getCI_3_j <- quantile(ca_sim_inv_j[3,1,],.025)
ca_getCI_u3_j <-quantile(ca_sim_inv_j[3,1,],.975)

ca_getCI_4_j <- quantile(ca_sim_inv_j[4,1,],.025)
ca_getCI_u4_j <-quantile(ca_sim_inv_j[4,1,],.975) 

ca_pred_df_j <- as.data.frame(ca_pred_ts_updated_h_j)
horizon <- c(1,2,3,4)
ca_pred_df_j <- cbind.data.frame(ca_pred_df_j, horizon,ca_test_j[1])
ggplot()+
  geom_line(mapping = aes( x = ca_pred_df_j$horizon, y = ca_pred_df_j$fit,color = "Observed data"))+
  geom_line(mapping = aes( x = ca_pred_df_j$horizon, y = ca_pred_df_j$wili_ca, color = "Forecasted data"))+
  geom_ribbon(aes(ymin = c(ca_getCI_j,ca_getCI_2_j,ca_getCI_3_j,ca_getCI_4_j), 
                  ymax=c(ca_getCI_u_j, ca_getCI_u2_j,ca_getCI_u3_j, ca_getCI_u4_j),
                  x=ca_pred_df_j$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")

# March 2019 (until week8 of 2019)

ca <- c(states[5], quidel_avail[5])

ca_quidel_m <- full_data[1:438, ca]
ca_quidel_m[1][439:442,] <- NA
ca_test_m <- full_data[431:434,ca]
ca_quidel_transform_m <- BoxCox(ca_quidel_m + 0.05, 0)
ca_quidel_ts_m <- ts(ca_quidel_transform_m, start = c(2018,01), frequency = 52)
# model with identity H matrix
ca_nonH_m <- create_SSMmodel(ca_quidel_ts_m,
                             covariate = TRUE,
                             trend_degree = 1,
                             cycle_num = 2,
                             cycle_base_period = 52,
                             cycle_type = "common",
                             return_optim_args = TRUE)

ca_nonH_m$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_nonH_m$model$P1) <- 100  ##prior variance of the coefficient

ca_mod_m <- optim(
  fn = ca_nonH_m$calc_loglik,
  par = ca_nonH_m$init_pars,
  model =ca_nonH_m$model,
  control = list(fnscale = -1)
)
ca_fit_m <- ca_nonH_m$update_pars(pars = ca_mod_m$par, model = ca_nonH_m$model)
ca_pred_4_m <- predict(ca_fit_m, n.ahead = 4)
ca_pred_inv_m <- InvBoxCox(as.data.frame(ca_pred_4_m), lambda = 0) - 0.05

#get MSE
ca_sse_m <- (ca_pred_inv_m - ca_test_j)^2
mse1_m <- mean(unlist(ca_sse_m[1])) 
mse2_m <- mean(unlist(ca_sse_m[2]))
#update H matrix
ca_t_m <- create_SSMmodel(ca_quidel_ts_m,
                          covariate = TRUE,
                          init_value_H = c(mse1_m,mse2_m),
                          
                          trend_degree = 1,
                          
                          cycle_num = 2,
                          cycle_base_period = 52,
                          cycle_type = "common",
                          return_optim_args = TRUE)

ca_t_m$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(ca_t_m$model$P1) <- 100

ca_opt_m <- optim(
  fn = ca_t_m$calc_loglik,
  par = ca_t_m$init_pars,
  model = ca_t_m$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



ca_fit_2_m <- ca_t_j$update_pars(pars = ca_opt_m$par, model = ca_t_m$model,init_value_H = c(mse1_m,mse2_m))
ca_pred_2_m <- predict(ca_fit_2_m, n.ahead = 4)
ca_pred_inv_2_m <- InvBoxCox(as.data.frame(ca_pred_2_m), lambda = 0) - 0.05
ca_pred_ts_updated_h_m <- ts(ca_pred_inv_2_m[1], start = c(2018,27), frequency=52)
ca_pred_ts_nonH_m <- ts(ca_pred_inv_m[1], start = c(2018,27), frequency=52)

autoplot(ts(ca_test[1], start = c(2018,27), frequency = 52))+
  autolayer(ca_pred_ts_updated_h)+
  autolayer(ca_pred_ts_nonH)

ca_simulate_incidence_m <- simulate_incidence(ca_fit_2_m,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
ca_sim_inv_m <- InvBoxCox(ca_simulate_incidence_m,lambda = 0)
ca_getCI_m <- quantile(ca_sim_inv_m[1,1,],.025)
ca_getCI_u_m <-quantile(ca_sim_inv_m[1,1,],.975) 

ca_getCI_2_m <- quantile(ca_sim_inv_m[2,1,],.025)
ca_getCI_u2_m <-quantile(ca_sim_inv_m[2,1,],.975) 

ca_getCI_3_m <- quantile(ca_sim_inv_m[3,1,],.025)
ca_getCI_u3_m <-quantile(ca_sim_inv_m[3,1,],.975)

ca_getCI_4_m <- quantile(ca_sim_inv_m[4,1,],.025)
ca_getCI_u4_m <-quantile(ca_sim_inv_m[4,1,],.975) 

ca_pred_df_m <- as.data.frame(ca_pred_ts_updated_h_m)
horizon <- c(1,2,3,4)
ca_pred_df_m <- cbind.data.frame(ca_pred_df_m, horizon,ca_test_m[1])
ggplot()+
  geom_line(mapping = aes( x = ca_pred_df_m$horizon, y = ca_pred_df_m$fit,color = "Observed data"))+
  geom_line(mapping = aes( x = ca_pred_df_m$horizon, y = ca_pred_df_m$wili_ca, color = "Forecasted data"))+
  geom_ribbon(aes(ymin = c(ca_getCI_m,ca_getCI_2_m,ca_getCI_3_m,ca_getCI_4_m), 
                  ymax=c(ca_getCI_u_m, ca_getCI_u2_m,ca_getCI_u3_m, ca_getCI_u4_m),
                  x=ca_pred_df_m$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")

