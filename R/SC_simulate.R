# SC

# July 2018 (until week 26 of 2018)
sc <- c(states[37], quidel_avail[37])

sc_quidel <- training_data[1:404, sc]
sc_quidel[1][404,] <- NA
sc_test <- training_data[405:408,sc]
sc_quidel_transform <- BoxCox(sc_quidel + 0.05, 0)
sc_quidel_ts <- ts(sc_quidel_transform, start = c(2018,01), frequency = 52)
# model with identity H matrix
sc_nonH <- create_SSMmodel(sc_quidel_ts,
                           covariate = TRUE,
                           trend_degree = 1,
                           cycle_num = 2,
                           cycle_base_period = 52,
                           cycle_type = "common",
                           return_optim_args = TRUE)

sc_nonH$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_nonH$model$P1) <- 100  ##prior variance of the coefficient

sc_mod <- optim(
  fn = sc_nonH$calc_loglik,
  par = sc_nonH$init_pars,
  model =sc_nonH$model,
  control = list(fnscale = -1)
)
sc_fit <- sc_nonH$update_pars(pars = sc_mod$par, model = sc_nonH$model)
sc_pred_4 <- predict(sc_fit, n.ahead = 4)
sc_pred_inv <- InvBoxCox(as.data.frame(sc_pred_4), lambda = 0) - 0.05

#get MSE
sc_sse <- (sc_pred_inv - sc_test)^2
mse1 <- mean(unlist(sc_sse[1])) 
mse2 <- mean(unlist(sc_sse[2]))
#update H matrix
sc_t <- create_SSMmodel(sc_quidel_ts,
                        covariate = TRUE,
                        init_value_H = c(mse1,mse2),
                        
                        trend_degree = 1,
                        
                        cycle_num = 2,
                        cycle_base_period = 52,
                        cycle_type = "common",
                        return_optim_args = TRUE)

sc_t$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_t$model$P1) <- 100

sc_opt <- optim(
  fn = sc_t$calc_loglik,
  par = sc_t$init_pars,
  model = sc_t$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



sc_fit_2 <- sc_t$update_pars(pars = sc_opt$par, model = sc_t$model,init_value_H = c(mse1,mse2))
sc_pred_2 <- predict(sc_fit_2, n.ahead = 4)
sc_pred_inv_2 <- InvBoxCox(as.data.frame(sc_pred_2), lambda = 0) - 0.05
sc_pred_ts_updated_h <- ts(sc_pred_inv_2[1], start = c(2018,27), frequency=52)
sc_pred_ts_nonH <- ts(sc_pred_inv[1], start = c(2018,27), frequency=52)

autoplot(ts(sc_test[1], start = c(2018,27), frequency = 52))+
  autolayer(sc_pred_ts_updated_h)+
  autolayer(sc_pred_ts_nonH)

sc_simulate_incidence <- simulate_incidence(sc_fit_2,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
sc_sim_inv <- InvBoxCox(sc_simulate_incidence,lambda = 0)
sc_getCI <- quantile(sc_sim_inv[1,1,],.025)
sc_getCI_u <-quantile(sc_sim_inv[1,1,],.975) 

sc_getCI_2 <- quantile(sc_sim_inv[2,1,],.025)
sc_getCI_u2 <-quantile(sc_sim_inv[2,1,],.975) 

sc_getCI_3 <- quantile(sc_sim_inv[3,1,],.025)
sc_getCI_u3 <-quantile(sc_sim_inv[3,1,],.975)

sc_getCI_4 <- quantile(sc_sim_inv[4,1,],.025)
sc_getCI_u4 <-quantile(sc_sim_inv[4,1,],.975) 

sc_pred_df <- as.data.frame(sc_pred_ts_updated_h)
horizon <- c(1,2,3,4)
sc_pred_df <- cbind.data.frame(sc_pred_df, horizon,sc_test[1])
ggplot()+
  geom_line(mapping = aes( x = sc_pred_df$horizon, y = sc_pred_df$fit, color = "Observed"))+
  geom_line(mapping = aes( x = sc_pred_df$horizon, y = sc_pred_df$wili_sc, color  = "Forecasted"))+
  geom_ribbon(aes(ymin = c(sc_getCI,sc_getCI_2,sc_getCI_3,sc_getCI_4), 
                  ymax=c(sc_getCI_u, sc_getCI_u2,sc_getCI_u3, sc_getCI_u4),
                  x=sc_pred_df$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")


# Octorber 2018



sc_quidel_o <- full_data[1:417, sc]
sc_quidel_o[1][418,] <- NA
sc_test_o <- full_data[418:421,sc]
sc_quidel_transform_o <- BoxCox(sc_quidel_o + 0.05, 0)
sc_quidel_ts_o <- ts(sc_quidel_transform_o, start = c(2018,01), frequency = 52)
# model with identity H matrix
sc_nonH_o <- create_SSMmodel(sc_quidel_ts_o,
                             covariate = TRUE,
                             trend_degree = 1,
                             cycle_num = 2,
                             cycle_base_period = 52,
                             cycle_type = "common",
                             return_optim_args = TRUE)

sc_nonH_o$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_nonH_o$model$P1) <- 100  ##prior variance of the coefficient

sc_mod_o <- optim(
  fn = sc_nonH_o$calc_loglik,
  par = sc_nonH_o$init_pars,
  model =sc_nonH_o$model,
  control = list(fnscale = -1)
)
sc_fit_o <- sc_nonH_o$update_pars(pars = sc_mod_o$par, model = sc_nonH_o$model)
sc_pred_4_o <- predict(sc_fit_o, n.ahead = 4)
sc_pred_inv_o <- InvBoxCox(as.data.frame(sc_pred_4_o), lambda = 0) - 0.05

#get MSE
sc_sse_o <- (sc_pred_inv_o - sc_test_o)^2
mse1_o <- mean(unlist(sc_sse_o[1])) 
mse2_o <- mean(unlist(sc_sse_o[2]))
#update H matrix
sc_t_o <- create_SSMmodel(sc_quidel_ts_o,
                          covariate = TRUE,
                          init_value_H = c(mse1_o,mse2_o),
                          
                          trend_degree = 1,
                          
                          cycle_num = 2,
                          cycle_base_period = 52,
                          cycle_type = "common",
                          return_optim_args = TRUE)

sc_t_o$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_t_o$model$P1) <- 100

sc_opt_o <- optim(
  fn = sc_t_o$calc_loglik,
  par = sc_t_o$init_pars,
  model = sc_t_o$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



sc_fit_2_o <- sc_t_o$update_pars(pars = sc_opt_o$par, model = sc_t_o$model,init_value_H = c(mse1_o,mse2_o))
sc_pred_2_o <- predict(sc_fit_2_o, n.ahead = 4)
sc_pred_inv_2_o <- InvBoxCox(as.data.frame(sc_pred_2_o), lambda = 0) - 0.05
sc_pred_ts_updated_h_o <- ts(sc_pred_inv_2_o[1], start = c(2018,27), frequency=52)
sc_pred_ts_nonH_o <- ts(sc_pred_inv_o[1], start = c(2018,27), frequency=52)

autoplot(ts(sc_test[1], start = c(2018,27), frequency = 52))+
  autolayer(sc_pred_ts_updated_h)+
  autolayer(sc_pred_ts_nonH)

sc_simulate_incidence_o <- simulate_incidence(sc_fit_2_o,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
sc_sim_inv_o <- InvBoxCox(sc_simulate_incidence_o,lambda = 0)
sc_getCI_o <- quantile(sc_sim_inv_o[1,1,],.025)
sc_getCI_u_o <-quantile(sc_sim_inv_o[1,1,],.975) 

sc_getCI_2_o <- quantile(sc_sim_inv_o[2,1,],.025)
sc_getCI_u2_o <-quantile(sc_sim_inv_o[2,1,],.975) 

sc_getCI_3_o <- quantile(sc_sim_inv_o[3,1,],.025)
sc_getCI_u3_o <-quantile(sc_sim_inv_o[3,1,],.975)

sc_getCI_4_o <- quantile(sc_sim_inv_o[4,1,],.025)
sc_getCI_u4_o <-quantile(sc_sim_inv_o[4,1,],.975) 

sc_pred_df_o <- as.data.frame(sc_pred_ts_updated_h_o)
horizon <- c(1,2,3,4)
sc_pred_df_o <- cbind.data.frame(sc_pred_df_o, horizon,sc_test_o[1])
ggplot()+
  geom_line(mapping = aes( x = sc_pred_df_o$horizon, y = sc_pred_df_o$fit,color = "Observed"))+
  geom_line(mapping = aes( x = sc_pred_df_o$horizon, y = sc_pred_df_o$wili_sc, color = "Forecasted"))+
  geom_ribbon(aes(ymin = c(sc_getCI_o,sc_getCI_2_o,sc_getCI_3_o,sc_getCI_4_o), 
                  ymax=c(sc_getCI_u_o, sc_getCI_u2_o,sc_getCI_u3_o, sc_getCI_u4_o),
                  x=sc_pred_df_o$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")


# January 2019 (until final week of 2018)


sc_quidel_j <- full_data[1:430, sc]
sc_quidel_j[1][430,] <- NA
sc_test_j <- full_data[431:434,sc]
sc_quidel_transform_j <- BoxCox(sc_quidel_j + 0.05, 0)
sc_quidel_ts_j <- ts(sc_quidel_transform_j, start = c(2018,01), frequency = 52)
# model with identity H matrix
sc_nonH_j <- create_SSMmodel(sc_quidel_ts_j,
                             covariate = TRUE,
                             trend_degree = 1,
                             cycle_num = 2,
                             cycle_base_period = 52,
                             cycle_type = "common",
                             return_optim_args = TRUE)

sc_nonH_j$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_nonH_j$model$P1) <- 100  ##prior variance of the coefficient

sc_mod_j <- optim(
  fn = sc_nonH_j$calc_loglik,
  par = sc_nonH_j$init_pars,
  model =sc_nonH_j$model,
  control = list(fnscale = -1)
)
sc_fit_j <- sc_nonH_j$update_pars(pars = sc_mod_j$par, model = sc_nonH_j$model)
sc_pred_4_j <- predict(sc_fit_j, n.ahead = 4)
sc_pred_inv_j <- InvBoxCox(as.data.frame(sc_pred_4_j), lambda = 0) - 0.05

#get MSE
sc_sse_j <- (sc_pred_inv_j - sc_test_j)^2
mse1_j <- mean(unlist(sc_sse_j[1])) 
mse2_j <- mean(unlist(sc_sse_j[2]))
#update H matrix
sc_t_j <- create_SSMmodel(sc_quidel_ts_j,
                          covariate = TRUE,
                          init_value_H = c(mse1_j,mse2_j),
                          
                          trend_degree = 1,
                          
                          cycle_num = 2,
                          cycle_base_period = 52,
                          cycle_type = "common",
                          return_optim_args = TRUE)

sc_t_j$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_t_j$model$P1) <- 100

sc_opt_j <- optim(
  fn = sc_t_j$calc_loglik,
  par = sc_t_j$init_pars,
  model = sc_t_j$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



sc_fit_2_j <- sc_t_j$update_pars(pars = sc_opt_j$par, model = sc_t_j$model,init_value_H = c(mse1_j,mse2_j))
sc_pred_2_j <- predict(sc_fit_2_j, n.ahead = 4)
sc_pred_inv_2_j <- InvBoxCox(as.data.frame(sc_pred_2_j), lambda = 0) - 0.05
sc_pred_ts_updated_h_j <- ts(sc_pred_inv_2_j[1], start = c(2018,27), frequency=52)
sc_pred_ts_nonH_j <- ts(sc_pred_inv_j[1], start = c(2018,27), frequency=52)

autoplot(ts(sc_test[1], start = c(2018,27), frequency = 52))+
  autolayer(sc_pred_ts_updated_h)+
  autolayer(sc_pred_ts_nonH)

sc_simulate_incidence_j <- simulate_incidence(sc_fit_2_j,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
sc_sim_inv_j <- InvBoxCox(sc_simulate_incidence_j,lambda = 0)
sc_getCI_j <- quantile(sc_sim_inv_j[1,1,],.025)
sc_getCI_u_j <-quantile(sc_sim_inv_j[1,1,],.975) 

sc_getCI_2_j <- quantile(sc_sim_inv_j[2,1,],.025)
sc_getCI_u2_j <-quantile(sc_sim_inv_j[2,1,],.975) 

sc_getCI_3_j <- quantile(sc_sim_inv_j[3,1,],.025)
sc_getCI_u3_j <-quantile(sc_sim_inv_j[3,1,],.975)

sc_getCI_4_j <- quantile(sc_sim_inv_j[4,1,],.025)
sc_getCI_u4_j <-quantile(sc_sim_inv_j[4,1,],.975) 

sc_pred_df_j <- as.data.frame(sc_pred_ts_updated_h_j)
horizon <- c(1,2,3,4)
sc_pred_df_j <- cbind.data.frame(sc_pred_df_j, horizon,sc_test_j[1])
ggplot()+
  geom_line(mapping = aes( x = sc_pred_df_j$horizon, y = sc_pred_df_j$fit,color = "Observed"))+
  geom_line(mapping = aes( x = sc_pred_df_j$horizon, y = sc_pred_df_j$wili_sc, color = "Forecasted"))+
  geom_ribbon(aes(ymin = c(sc_getCI_j,sc_getCI_2_j,sc_getCI_3_j,sc_getCI_4_j), 
                  ymax=c(sc_getCI_u_j, sc_getCI_u2_j,sc_getCI_u3_j, sc_getCI_u4_j),
                  x=sc_pred_df_j$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")

# March 2019 (until week8 of 2019)



sc_quidel_m <- full_data[1:438, sc]
sc_quidel_m[1][439:442,] <- NA
sc_test_m <- full_data[431:434,sc]
sc_quidel_transform_m <- BoxCox(sc_quidel_m + 0.05, 0)
sc_quidel_ts_m <- ts(sc_quidel_transform_m, start = c(2018,01), frequency = 52)
# model with identity H matrix
sc_nonH_m <- create_SSMmodel(sc_quidel_ts_m,
                             covariate = TRUE,
                             trend_degree = 1,
                             cycle_num = 2,
                             cycle_base_period = 52,
                             cycle_type = "common",
                             return_optim_args = TRUE)

sc_nonH_m$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_nonH_m$model$P1) <- 100  ##prior variance of the coefficient

sc_mod_m <- optim(
  fn = sc_nonH_m$calc_loglik,
  par = sc_nonH_m$init_pars,
  model =sc_nonH_m$model,
  control = list(fnscale = -1)
)
sc_fit_m <- sc_nonH_m$update_pars(pars = sc_mod_m$par, model = sc_nonH_m$model)
sc_pred_4_m <- predict(sc_fit_m, n.ahead = 4)
sc_pred_inv_m <- InvBoxCox(as.data.frame(sc_pred_4_m), lambda = 0) - 0.05

#get MSE
sc_sse_m <- (sc_pred_inv_m - sc_test_j)^2
mse1_m <- mean(unlist(sc_sse_m[1])) 
mse2_m <- mean(unlist(sc_sse_m[2]))
#update H matrix
sc_t_m <- create_SSMmodel(sc_quidel_ts_m,
                          covariate = TRUE,
                          init_value_H = c(mse1_m,mse2_m),
                          
                          trend_degree = 1,
                          
                          cycle_num = 2,
                          cycle_base_period = 52,
                          cycle_type = "common",
                          return_optim_args = TRUE)

sc_t_m$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(sc_t_m$model$P1) <- 100

sc_opt_m <- optim(
  fn = sc_t_m$calc_loglik,
  par = sc_t_m$init_pars,
  model = sc_t_m$model,
  init_value_H = "Yes",
  control = list(fnscale = -1)
)



sc_fit_2_m <- sc_t_j$update_pars(pars = sc_opt_m$par, model = sc_t_m$model,init_value_H = c(mse1_m,mse2_m))
sc_pred_2_m <- predict(sc_fit_2_m, n.ahead = 4)
sc_pred_inv_2_m <- InvBoxCox(as.data.frame(sc_pred_2_m), lambda = 0) - 0.05
sc_pred_ts_updated_h_m <- ts(sc_pred_inv_2_m[1], start = c(2018,27), frequency=52)
sc_pred_ts_nonH_m <- ts(sc_pred_inv_m[1], start = c(2018,27), frequency=52)

autoplot(ts(sc_test[1], start = c(2018,27), frequency = 52))+
  autolayer(sc_pred_ts_updated_h)+
  autolayer(sc_pred_ts_nonH)

sc_simulate_incidence_m <- simulate_incidence(sc_fit_2_m,  backcast_steps = 0, forecast_steps = 4, n_sim = 1000) 
sc_sim_inv_m <- InvBoxCox(sc_simulate_incidence_m,lambda = 0)
sc_getCI_m <- quantile(sc_sim_inv_m[1,1,],.025)
sc_getCI_u_m <-quantile(sc_sim_inv_m[1,1,],.975) 

sc_getCI_2_m <- quantile(sc_sim_inv_m[2,1,],.025)
sc_getCI_u2_m <-quantile(sc_sim_inv_m[2,1,],.975) 

sc_getCI_3_m <- quantile(sc_sim_inv_m[3,1,],.025)
sc_getCI_u3_m <-quantile(sc_sim_inv_m[3,1,],.975)

sc_getCI_4_m <- quantile(sc_sim_inv_m[4,1,],.025)
sc_getCI_u4_m <-quantile(sc_sim_inv_m[4,1,],.975) 

sc_pred_df_m <- as.data.frame(sc_pred_ts_updated_h_m)
horizon <- c(1,2,3,4)
sc_pred_df_m <- cbind.data.frame(sc_pred_df_m, horizon,sc_test_m[1])
ggplot()+
  geom_line(mapping = aes( x = sc_pred_df_m$horizon, y = sc_pred_df_m$fit,color = "Observed"))+
  geom_line(mapping = aes( x = sc_pred_df_m$horizon, y = sc_pred_df_m$wili_sc, color = "Forecasted"))+
  geom_ribbon(aes(ymin = c(sc_getCI_m,sc_getCI_2_m,sc_getCI_3_m,sc_getCI_4_m), 
                  ymax=c(sc_getCI_u_m, sc_getCI_u2_m,sc_getCI_u3_m, sc_getCI_u4_m),
                  x=sc_pred_df_m$horizon, fill = "Simulation range"), alpha = 0.3)+
  xlab("Horizon")+
  ylab("ILI")

