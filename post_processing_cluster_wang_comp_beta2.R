#uses list from wang_comp_beta2_mycomp.R

#taking list from environment
v <- wang_comp_beta2_list

###POST PROCESSING

post_processing_beta2 <- function(v, n_participants, beta2a, beta2b, alpha, n_visits, delta, visit_frequency, runs){
  #extract beta_hat_trues from v and convert to vector
  beta_hat_trues <- unlist(lapply(v, function(x) x$beta_hat_true), use.names = FALSE)
  
  #extract dataframes of fitted cox models from v and convert to one dataframe (fitted_cox_models)
  fitted_cox_models <- bind_rows(lapply(v, function(x) data.frame(result = x$result)))
  
  #extract dataframes of log rank tests' results from v and convert to one dataframe (logranktests)
  lrts <- bind_rows(lapply(v, function(x) data.frame(lrt = x$lrt)))
  
  #extract proportion of r censored observations in true data from v and convert to one column (prop_r_cens_trues)
  prop_r_cens_true <- data.frame(prop_r_cens_true = unlist((lapply(v, function(x) x$prop_r_cens_true))))
  
  #extract proportion of r censored observations in obs data from v and convert to one column (prop_r_cens_obss)
  prop_r_cens_obs <- data.frame(prop_r_cens_obs = unlist((lapply(v, function(x) x$prop_r_cens_obs))))
  
  #extract proportion of underestimated intervals from v and convert to one column (prop_underest)
  prop_underest <- data.frame(prop_underest = unlist((lapply(v, function(x) x$prop_underest))))
  
  #extract proportion of overestimated intervals from v and convert to one column (prop_overest)
  prop_overest <- data.frame(prop_overest = unlist((lapply(v, function(x) x$prop_overest))))
  
  #add p check column (checks if p is less than alpha) to log rank tests' results dataframe
  lrts <- lrts %>%
    mutate(logrank_p_check = as.numeric(lrt.logrank_p<alpha))
  
  #mean of estimated true betas:
  beta <- mean(beta_hat_trues)
  
  #add estimated true beta, confidence intervals and check if they cover true beta
  fitted_cox_models <- fitted_cox_models %>%
    mutate(beta = beta,
           lci = as.numeric(result.beta_hat.group - qnorm(0.975)*result.se_hat),
           uci = as.numeric(result.beta_hat.group + qnorm(0.975)*result.se_hat),
           check_cover = as.numeric((lci < beta) & (uci > beta)))
  #add bias, ese, mean estimated standard error, coverage probability for each data type
  fitted_cox_models <- fitted_cox_models %>%
    group_by(result.data_type) %>%
    mutate(bias_beta_hat = mean(result.beta_hat.group) - beta,
           emprical_standard_error = sd(result.beta_hat.group),
           mean_estimated_standard_error = mean(result.se_hat),
           coverage_probability = mean(check_cover))
  
  #fuse log rank results and fitted cox models
  full_data <- cbind(lrts, fitted_cox_models, prop_r_cens_true, prop_r_cens_obs, prop_underest, prop_overest)
  
  #calculate the proportion that p is less than alpha per censor/imputation type for log rank test
  full_data <- full_data %>%
    group_by(result.data_type) %>%
    mutate(prop_null_rejected_lrt = mean(logrank_p_check),
           type1error_power_hr_check = as.numeric(result.hr_p < 0.05),
           type1error_power_hr = mean(type1error_power_hr_check))
  
  #add "data_run" column: easier to see which rows correspond to each simulated dataset
  data_run <- data.frame(data_run = rep(1:runs, each = 8))
  #add "beta2a" and "beta2b" column, visit frequency, number of visits, number of participants per arm, delta (wiggle room for visits)
  study_structure <- data.frame(beta2a = beta2a, beta2b = beta2b, visit_frequency = visit_frequency, n_visits = n_visits, participants_per_arm = n_participants, delta = delta)
  full_data <- cbind(data_run, study_structure, full_data)
  
  return(full_data)
}

wang_comp_beta2 <- post_processing_beta2(v = v, n_participants = 25, beta2a = 6.94, beta2b = 11.34, alpha = 0.05, n_visits = 12, delta = 0.5, visit_frequency = "4weeks", runs = 1000)


