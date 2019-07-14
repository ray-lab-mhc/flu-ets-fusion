## Cross-validation scores:
## Each model is the combination of four components listed below:
## a) Local level 
## b) Cycle
##   - cycle of frequency 52 
##   - cycle of frequency 52 and frequency 26
##   - cycle of frequency 52, frequency 26, and frequency 13
## c) Dummy variables:
##   - None
##   - Holiday effects: Christmas
## d) Transformation
##   - No
##   - log (lambda = 0)
##   - BoxCox 


score_1 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 37,inits = c(inits_ls_c,-4.139273), update_func = uplik_lc)

score_2 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 37,inits = c(inits_ls_c,-4.139273, log(0.01)), update_func = uplik_model_lcc)

score_3 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 37, inits = c(inits_ls_c,-4.139273, log(0.01), log(0.1)), update_func = uplik_model_lc3)

score_4 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 37,inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_5 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 37,inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_6 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 37,inits = c(inits_d2[1:19]), update_func = uplik_m2)



score_7 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 0,inits = c(inits_d2[1:16]), update_func = uplik_lc)

score_8 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 0,inits = c(inits_d2[1:17]), update_func = uplik_model_lcc)

score_9 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 0,inits = c(inits_d2[1:18]), update_func = uplik_model_lc3)

score_10 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 0,inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_11 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 0,inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_12 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 0,inits = c(inits_d2[1:19]), update_func = uplik_m2)


boxcox <- BoxCox.lambda(flu_ts)
score_13 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = boxcox,inits = c(inits_d2[1:16]), update_func = uplik_lcc_mod)

score_14 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = boxcox,inits = c(inits_d2[1:17]), update_func = uplik_model_lcc)

score_15 <- CV_es(flu_ts, missing = FALSE, is_dummy = "no", dummy = 1, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = boxcox,inits = c(inits_d2[1:18]), update_func = uplik_model_lc3)

score_16 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = boxcox,inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_17 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = boxcox,inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_18 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = flu_ts2[,6:7], window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = boxcox,inits = c(inits_d2[1:19]), update_func = uplik_m2)

## Christmas only
christmas <- ts(data.frame(flu_data2[,10]), start = c(2010, 40),frequency = 52)

score_19 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 37,inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_20 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 37,inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_21 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 37,inits = c(inits_d2[1:19]), update_func = uplik_m2)

score_22 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = 0, inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_23 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = 0, inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_24 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = 0, inits = c(inits_d2[1:19]), update_func = uplik_m2)

score_25 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  1, freq = 52, lamb = boxcox, inits = c(inits_d2[1:17]), update_func = uplik_lcc_d_mod)

score_26 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  2, freq = 52, lamb = boxcox, inits = c(inits_d2[1:18]), update_func = uplik_m1)

score_27 <- CV_es(flu_ts2[,1:5], missing = FALSE, is_dummy = "yes", dummy = christmas, window = 156, horizon = 4, seasonal = "cycle", cycle =  3, freq = 52, lamb = boxcox, inits = c(inits_d2[1:19]), update_func = uplik_m2)



