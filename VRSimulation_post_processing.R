# VR Manuscript Simulation Post-processing


#post processing function: produces a dataframe where each row corresponds to the results of analysis on one simulated dataset for one "data type" (true TTVR, imputed and interval-censored TTVR for both observed and true data)
#returns a dataframe with information regarding simulation structure, misclassification frequency, frequency of right-censored observation, and RMST, Cox PH model, and LRT performance metrics (power/type one error, bias, 95% coverage probability, estimated and empirical standard errors)
#see README for a description of all variables returned.
post_processing <- function(v, n_participants, beta2a, beta2b, alpha, n_visits, delta, tau, runs){
  #combine all elements of the list in one dataframe (row binding)
  full_data <- do.call(rbind, v)
  
  #add "data_run" column: easier to see which rows correspond to each simulated dataset. For example, if data_run = 50, then that row is a result from the 50th run of the simulation.
  data_run <- data.frame(data_run = rep(1:runs, each = 9))
  #add "beta2a" and "beta2b" column, number of visits, number of participants per arm, delta (the wiggle room for visits)
  study_structure <- data.frame(beta2a = beta2a, beta2b = beta2b, n_visits = n_visits, participants_per_arm = n_participants, delta = delta, tau = tau)
  #add "data_type_description" column
  data_type_description <- data.frame(data_type_description = 
                            rep(c("truth","r_imp_true", "r_imp_obs", "m_imp_true", "m_imp_obs", 
                                  "l_imp_true", "l_imp_obs", "IC_true", "IC_obs"), times = runs))
  
  full_data <- cbind(data_run, data_type_description, study_structure, full_data)
  
  #-----------------------------------------------------------------------------
  # RMST post-processing
  #-----------------------------------------------------------------------------
  #calculate RMST "truth"
  true_RMST_diff <- mean(full_data$diff_RMST[full_data$data_type_description == "truth" &
                                     !is.na(full_data$diff_RMST)])
  true_RMST_arm1 <- mean(full_data$RMST1[full_data$data_type_description == "truth" &
                                      !is.na(full_data$diff_RMST)])
  true_RMST_arm2 <- mean(full_data$RMST2[full_data$data_type_description == "truth" &
                                      !is.na(full_data$diff_RMST)])
  
  #check if confidence intervals cover the true difference in RMSTs, calculate proportion of runs where prop. hazards assumption is violated
  full_data <- full_data %>%
    mutate(check_cover_true_RMST_diff = as.numeric((lci_diff_RMST < true_RMST_diff) & (uci_diff_RMST > true_RMST_diff)),
           ph_violated = as.numeric(prop_haz_p < 0.05))
  
  prop_ph_violated <- mean(full_data$ph_violated)
  mean_propVRunder8_arm0 <- mean(full_data$propVRunder8_arm0)
  mean_propVRunder8_arm1 <- mean(full_data$propVRunder8_arm1)
  
  full_data <- full_data %>%
    mutate(prop_ph_violated = prop_ph_violated,
           mean_propVRunder8_arm0 = mean_propVRunder8_arm0,
           mean_propVRunder8_arm1 = mean_propVRunder8_arm1)
  
  #add diff bias, RMST "truth", probability that confidence intervals cover 0 (type 1 error rate/hypothesis test), mean estimated standard error, empirical standard error, coverage probability
  #add bias and RMST truth for arm 1
  #add bias and RMST truth for arm 2
  full_data <- full_data %>%
    group_by(data_type_description) %>%
    mutate(bias_RMST_diff = mean(diff_RMST[!is.na(diff_RMST)]) - true_RMST_diff,
           true_RMST_diff = true_RMST_diff,
           RMST_CI_covers_0_prob = mean(covers_zero_CI_diff_RMST[!is.na(covers_zero_CI_diff_RMST)]),
           RMST_prop_null_rejected = 1-RMST_CI_covers_0_prob,
           RMST_mean_est_se = mean(se_diff_RMST[!is.na(se_diff_RMST)]),
           RMST_emp_se = sd(diff_RMST[!is.na(diff_RMST)]),
           RMST_CI_covers_true_RMST_diff_prob = mean(check_cover_true_RMST_diff[!is.na(check_cover_true_RMST_diff)]),
           # arm 1
           RMST_bias_arm1 = mean(RMST1[!is.na(RMST1)]) - true_RMST_arm1,
           true_RMST1 = true_RMST_arm1,
           # arm 2
           RMST_bias_arm2 = mean(RMST2[!is.na(RMST2)]) - true_RMST_arm2,
           true_RMST2 = true_RMST_arm2)
  
  
  #add number of iterations (number of unique data_runs) that contained at least one NA
  n_iters_na <- full_data %>%
    group_by(data_run) %>%
    summarise(any_na = any(is.na(diff_RMST)), .groups = "drop") %>%
    summarise(n_iters_na = sum(any_na)) %>%
    pull(n_iters_na)
  
  full_data <- full_data %>%
    mutate(RMST_n_iters_na = n_iters_na)
  
  #-----------------------------------------------------------------------------
  # Cox PH model post-processing
  #-----------------------------------------------------------------------------
  
  #mean of estimated true betas:
  names(full_data)[names(full_data) == "logHR_hat.group"] <- "logHR_hat"
  true_logHR <- mean(full_data$logHR_hat[full_data$data_type_description == "truth"])
  
  full_data <- full_data %>%
    mutate(true_logHR = true_logHR)
  
  #add confidence intervals and check if they cover true log HR
  full_data <- full_data %>%
    mutate(logHR_lci = as.numeric(logHR_hat - qnorm(0.975)*logHR_se_hat),
           logHR_uci = as.numeric(logHR_hat + qnorm(0.975)*logHR_se_hat),
           logHR_CI_check_cover = as.numeric((logHR_lci < true_logHR) & (logHR_uci > true_logHR)))
  
  #add bias, ese, mean estimated standard error, coverage probability for each data type
  full_data <- full_data %>%
    group_by(data_type_description) %>%
    mutate(bias_logHR = mean(logHR_hat) - true_logHR,
           logHR_emp_se = sd(logHR_hat),
           logHR_mean_est_se = mean(logHR_se_hat),
           logHR_covers_true_logHR_prob = mean(logHR_CI_check_cover),
           null_rejected_logHR_check = as.numeric(logHR_p < 0.05),
           logHR_prop_null_rejected = mean(null_rejected_logHR_check))
  
  #-----------------------------------------------------------------------------
  # log-rank test post-processing
  #-----------------------------------------------------------------------------
  
  #add p check column (checks if p is less than alpha) to log rank tests' results dataframe
  full_data <- full_data %>%
    mutate(logrank_p_check = as.numeric(logrank_p<alpha))
  
  #calculate the proportion that p is less than alpha per censor/imputation type for log rank test
  full_data <- full_data %>%
    group_by(data_type_description) %>%
    mutate(LRT_prop_null_rejected = mean(logrank_p_check))
  
  
  #-----------------------------------------------------------------------------
  # cleaning column names
  #-----------------------------------------------------------------------------
  
  full_data <- full_data %>% 
    rename(prop_r_cens_true = `df$prop_r_cens_true[1:9]`, 
           prop_r_cens_obs = `df$prop_r_cens_obs[1:9]`,
           prop_underest = `df$prop_underest[1:9]`,
           prop_overest = `df$prop_overest[1:9]`
           )
  
  
  return(full_data)
}

processed_sim <- post_processing(v = sim_results, n_participants = 100, beta2a = 6.94, beta2b = 11.34, alpha = 0.05, n_visits = 48, delta = 0.5, tau = 8, runs = 1000)





