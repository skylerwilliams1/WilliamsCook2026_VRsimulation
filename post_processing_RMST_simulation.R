#taking list from environment
outRMST_beta2_6.94_9.14_12visits_4weeks_25perarm_tau8 <- RMST_beta2_6.94_9.14_12visits_4weeks_25perarm_tau8
v <- outRMST_beta2_6.94_9.14_12visits_4weeks_25perarm_tau8

###POST PROCESSING

post_processing_beta2_RMST <- function(v, n_participants, beta2a, beta2b, alpha, n_visits, delta, visit_frequency, tau, runs){
  #combine all elements of the list in one dataframe (row binding)
  full_data <- do.call(rbind, v)
  
  #add "data_run" column: easier to see which rows correspond to each simulated dataset
  data_run <- data.frame(data_run = rep(1:runs, each = 9))
  #add "beta2a" and "beta2b" column, visit frequency, number of visits, number of participants per arm, delta (wiggle room for visits)
  study_structure <- data.frame(beta2a = beta2a, beta2b = beta2b, visit_frequency = visit_frequency, n_visits = n_visits, participants_per_arm = n_participants, delta = delta, tau = tau)
  #add "data_type" column
  data_type <- data.frame(data_type = 
                            rep(c("truth","r_imp_true", "r_imp_obs", "m_imp_true", "m_imp_obs", 
                                "l_imp_true", "l_imp_obs", "IC_true", "IC_obs"), times = runs))
  
  full_data <- cbind(data_run, data_type, study_structure, full_data)
  
  #calculate RMST "truth" (only when truth is not NA)
  true_diff <- mean(full_data$diff[full_data$data_type == "truth" &
                                     !is.na(full_data$diff)])
  true_arm1 <- mean(full_data$RMST1[full_data$data_type == "truth" &
                                      !is.na(full_data$diff)])
  true_arm2 <- mean(full_data$RMST2[full_data$data_type == "truth" &
                                      !is.na(full_data$diff)])
  
  #check if lci and uci cover the true RMST diff, calculate proportion of times prop hazards assumption is violated

  
  full_data <- full_data %>%
    mutate(check_cover_true_diff = as.numeric((lci < true_diff) & (uci > true_diff)),
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
    group_by(data_type) %>%
    mutate(bias_diff = mean(diff[!is.na(diff)]) - true_diff,
           true_diff = true_diff,
           covers_0_prob = mean(covers_zero[!is.na(covers_zero)]),
           type1errorrate = 1-covers_0_prob,
           est_se = mean(se[!is.na(se)]),
           emp_se = sd(diff[!is.na(diff)]),
           covers_true_diff_prob = mean(check_cover_true_diff[!is.na(check_cover_true_diff)]),
           # arm 1
           bias_arm1 = mean(RMST1[!is.na(RMST1)]) - true_arm1,
           true_RMST1 = true_arm1,
           # arm 2
           bias_arm2 = mean(RMST2[!is.na(RMST2)]) - true_arm2,
           true_RMST2 = true_arm2)
  
  
  #add number of iterations (number of unique data_runs) that contained at least one NA
  n_iters_na <- full_data %>%
    group_by(data_run) %>%
    summarise(any_na = any(is.na(diff)), .groups = "drop") %>%
    summarise(n_iters_na = sum(any_na)) %>%
    pull(n_iters_na)
  
  full_data <- full_data %>%
    mutate(n_iters_na = n_iters_na)
  
  return(full_data)
}

outRMST_beta2_6.94_9.14_12visits_4weeks_25perarm_tau8 <- post_processing_beta2_RMST(v = v, n_participants = 25, beta2a = 6.94, beta2b = 9.14, alpha = 0.05, n_visits = 12, delta = 0.5, visit_frequency = "4weeks", tau = 8, runs = 1000)








