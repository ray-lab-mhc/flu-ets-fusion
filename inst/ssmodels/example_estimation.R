library(fluETS)
library(KFAS)
library(mvtnorm)

## A collection of updating functions that are required for the estimation of parameters.
flu_ts <- readRDS(file = "~/Documents/research/epi/flu/flu-ets-fusion/data/ny_flu_wiki.rds")

tm2 <- create_SSMmodel(
  ts_data = flu_ts %>%
    select(wili_jfk, wili_ny_upstate, wiki_common_cold, wiki_cough, wiki_influenza) %>%
    as.matrix(),
  regression_variables = c("christmas_week"),
  regression_data = flu_ts,
  trend_degree = 1,
  cycle_num = 2,
  cycle_base_period = 52,
  cycle_type = "common",
  return_optim_args = TRUE
)

tm2$model$P1inf[] <- 0  ##removes the diffuse initialization
diag(tm2$model$P1) <- 100  ##prior variance of the coefficient

model_fit_pars <- optim(
  fn = tm2$calc_loglik,
  par = tm2$init_pars,
  model = tm2$model,
  control = list(fnscale = -1)
)

model_fit <- tm2$update_pars(pars = model_fit_pars$par, model = tm2$model)

horizon <- 4

# newdata object contains system matrices for future times
# in particular, Z contains covariates which may differ at future times.
newdata_SSMmodel <- create_SSMmodel(
  ts_data = matrix(NA, nrow = horizon, ncol = attr(model_fit, "p")),
  regression_variables = c("christmas_week"),
  regression_data = flu_ts[11:14, ],
  trend_degree = 1,
  cycle_num = 2,
  cycle_base_period = 52,
  cycle_type = "common",
  return_optim_args = FALSE
)$model

newdata_SSMmodel <- tm2$update_pars(pars = model_fit_pars$par, model = newdata_SSMmodel)

sim_inc <- simulate_incidence(model_fit = model_fit, backcast_steps = 2, forecast_steps= 4, n_sim = 5, new_regression_data = flu_ts[11:14, ])


filter_results <- KFS(model_fit)

alpha_t <- t(filter_results$alphahat[nrow(filter_results$alphahat), ])

a1 <- model_fit$T[, , 1] %*% alpha_t
a2 <- model_fit$T[, , 1] %*% a1
a3 <- model_fit$T[, , 1] %*% a2
a4 <- model_fit$T[, , 1] %*% a3

y1 <- newdata_SSMmodel$Z[,,1] %*% a1
y2 <- newdata_SSMmodel$Z[,,2] %*% a2
y3 <- newdata_SSMmodel$Z[,,3] %*% a3
y4 <- newdata_SSMmodel$Z[,,4] %*% a4

predict(model_fit, newdata = newdata_SSMmodel) %>% bind_cols()

rbind(t(y1), t(y2), t(y3), t(y4))



