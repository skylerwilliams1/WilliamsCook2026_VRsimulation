library(survival)
library(tidyverse)
library(icenReg)
library(interval)
library(coin)

one_run_wang_comp_beta2 <- function(n_participants, beta2a, beta2b, alpha, n_visits, delta, visit_frequency, runs){
  
  #function that generates each visit time according to uniform dist (one participant)
  #units of returned vector are weeks
  visit_times <- function(){
    visits <- vector("numeric", n_visits)
    for (i in 2:n_visits){
      visits[i] <- runif(n = 1, min = 4*(i-1) - delta, max = 4*(i-1) + delta)
    }
    return(visits)
  }
  
  #function to simulate data for one group(returns sim_data)
  visits_vl_df <- function(n_visits, delta, n_participants, beta2){
    
    #function creates vector of visit times and corresponding viral loads based on model from Wang et al. 2020
    one_virus <- function(){
      #set fixed effects (these values are from Bing et al. 2020)
      beta1 <- 3.61
      #beta2 <- 6.94
      beta3 <- 2.12
      beta4 <- 1.72
      beta5 <- 0.12
      #random effects (Wang best model included random effects on beta1 and beta3, values taken from Bing paper)
      tao1 <- rnorm(n = 1, mean = 0, sd = 0.75)
      tao3 <- rnorm(n = 1, mean = 0, sd = 0.67)
      #combine random and fixed effects
      theta1 <- beta1 + tao1
      theta2 <- beta2
      theta3 <- beta3 + tao3
      theta4 <- beta4
      theta5 <- beta5
      #calculate viral load for each visit time, save loads in vector
      hloads <- c()
      yloads <- c()
      visits <- visit_times()
      for (i in seq_along(visits)){
        t <- visits[i]
        hload_t <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) #calculates h(), which is the viral load trajectory
        hloads <- c(hloads, hload_t)
        yload_t <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) + #calculates y(), which is h() + random error
          rnorm(1, mean = 0, sd = sqrt(0.26))
        yloads <- c(yloads, yload_t)
      }
      
      #function to calculate h viral load:
      h <- function(t){
        hvl <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t)))
        return(hvl)
      }
      #function to find true time of threshold crossing (threshold is 1000 copies/mL (3 log10(copies/mL)))
      find_true_time_of_crossing <- function(){
        times <- seq(0, visits[n_visits], 0.01)
        h_of_times <- h(times)
        check <- h_of_times > 3
        first_true <- which(check)[1]
        time_of_crossing <- times[first_true]
        return(time_of_crossing)
      }
      
      #function calculate y viral load (h + error):
      y <- function(t){
        yvl <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) +
          rnorm(1, mean = 0, sd = sqrt(0.26))
        return(yvl)
      }
      
      #function to calculate observed time of threshold crossing (threshold is 1000 copies/mL (3 log10(copies/mL)))
      find_obs_time_of_crossing <- function(){
        times <- seq(0, visits[n_visits], 0.01)
        y_of_times <- y(times)
        check <- y_of_times > 3
        first_true <- which(check)[1]
        time_of_crossing <- times[first_true]
        return(time_of_crossing)
      }
      
      true_time <- find_true_time_of_crossing()
      if (is.na(true_time)){
        true_time <- visits[n_visits]
      }
      obs_time <- find_obs_time_of_crossing()
      if (is.na(obs_time)){
        obs_time <- visits[n_visits]
      }
      #function to record L bound of true interval of crossing
      true_interval_L <- function(){
        check_L <-  visits > true_time #checks what visits are past the true time of crossing
        L_index <- which(check_L)[1] #checks which visit is the first visit past the true time of crossing
        true_L <- visits[L_index - 1] #the true L interval is the visit prior
        true_R <- visits[L_index] #the true R interval is the first visit past the true time of crossing
        if (true_time == visits[n_visits]){ #if the true time of crossing is right censored (occurs past the last visit)
          true_L <- visits[n_visits]
        }
        return(true_L)
      }
      
      #function to record R bound of true interval of crossing
      true_interval_R <- function(){
        check_L <-  visits > true_time
        L_index <- which(check_L)[1]
        true_L <- visits[L_index - 1]
        true_R <- visits[L_index]
        if (true_time == visits[n_visits]){
          true_R <- Inf
        }
        return(true_R)
      }
      
      #function to record L bound of observed interval of crossing
      observed_interval_L <- function(){
        check_L <-  yloads > 3
        L_index <- which(check_L)[1]
        
        if (!is.na(L_index)){
          if (L_index > 1){
            obs_L <- visits[L_index - 1]
            obs_R <- visits[L_index]
          }
          if (L_index == 1){
            obs_L <- -1
            obs_R <- -1
          }
        }
        if (is.na(L_index)){
          obs_L <- visits[n_visits]
        }
        return(obs_L)
      }
      
      #function to record R bound of observed interval of crossing
      observed_interval_R <- function(){
        check_L <-  yloads > 3
        L_index <- which(check_L)[1]
        
        if (!is.na(L_index)){
          if (L_index > 1){
            obs_L <- visits[L_index - 1]
            obs_R <- visits[L_index]
          }
          if (L_index == 1){
            obs_L <- -1
            obs_R <- -1
          }
        }
        if (is.na(L_index)){
          obs_R <- Inf
        }
        return(obs_R)
      }
      
      
      true_int_L <- true_interval_L()
      true_int_R <- true_interval_R()
      obs_int_L <- observed_interval_L()
      obs_int_R <- observed_interval_R()
      
      #function to add true status variable based on true interval
      find_true_status <- function(){
        if (true_int_R == Inf){
          true_status <- 0
        }
        else {
          true_status <- 1
        }
        return(true_status)
      }
      
      #function to add observed status variable based on observed interval
      find_obs_status <- function(){
        if (obs_int_R == Inf){
          obs_status <- 0
        }
        else {
          obs_status <- 1
        }
        return(obs_status)
      }
      
      true_status <- find_true_status()
      obs_status <- find_obs_status()
      
      #fuses everything into one vector
      visits_loads <- c(beta2, visits, hloads, yloads, true_time, obs_time, 
                        true_int_L, true_int_R, obs_int_L, obs_int_R, true_status, obs_status)
      return(visits_loads)
    }
    
    #creates df of visits and viral loads for n_participants. makes sure that n_participants are included in each dataset (eliminates and replaces cases where obs viral load is over threshold at time 0)
    virus_df <- c()
    difference <- n_participants
    while (difference > 0){
      for (i in 1:difference){
        indiv <- one_virus()
        if (indiv[7] == -1 | indiv[8] == -1 | indiv[9] == -1 | indiv[10] == -1){
          next
        }
        virus_df <- data.frame(rbind(virus_df, indiv))
      }
      difference <- n_participants - nrow(virus_df)
    }
    
    #function that renames columns of above df 
    label_visits_loads <- function(){
      colnames(virus_df) <- c(
        paste("beta2"),
        paste("visit", 1:n_visits, sep=""),
        paste("hload", 1:n_visits, sep=""),
        paste("yload", 1:n_visits, sep=""),
        paste("true_crossing_time"),
        paste("observed_crossing_time"),
        paste("true_interval_L"),
        paste("true_interval_R"),
        paste("obs_interval_L"),
        paste("obs_interval_R"),
        paste("true_status"),
        paste("obs_status")
      )
      return(virus_df)
    }
    
    #calls above function and saves result
    virus_df <- label_visits_loads()
    sim_data <- virus_df
    
    #returns final df of visit times and viral loads for each participant
    return(sim_data)
    
    
  }
  
  #function for imputing the data
  add_imputation <- function(sim_data){
    sim_data <- sim_data %>%
      mutate(l_imp_obs = obs_interval_L, 
             r_imp_obs = obs_interval_R,
             l_imp_true = true_interval_L,
             r_imp_true = true_interval_R)
    sim_data$r_imp_obs[sim_data$r_imp_obs == Inf] <- sim_data$l_imp_obs[sim_data$r_imp_obs == Inf]
    sim_data$r_imp_true[sim_data$r_imp_true == Inf] <- sim_data$l_imp_true[sim_data$r_imp_true == Inf]
    sim_data <- sim_data %>%
      mutate(m_imp_obs = rowMeans(select(.,r_imp_obs, l_imp_obs)),
             m_imp_true = rowMeans(select(.,r_imp_true, l_imp_true)))
    return(sim_data)
  }
  
  #function for calculating proportion of right-censored data and adding to sim_data
  proportion_right_censored <- function(sim_data){
    prop_r_cens_true <- mean(sim_data$true_status == 0)
    prop_r_cens_obs <- mean(sim_data$obs_status == 0)
    sim_data <- sim_data %>%
      mutate(prop_r_cens_true = prop_r_cens_true,
             prop_r_cens_obs = prop_r_cens_obs)
    return(sim_data)
  }
  
  #function for calculating the proportion of misclassified intervals in L and R directions
  proportion_misclassified <- function(sim_data){
    prop_underest <- mean(sim_data$obs_interval_L < sim_data$true_interval_L)
    prop_overest <- mean(sim_data$obs_interval_L > sim_data$true_interval_L)
    sim_data <- sim_data %>%
      mutate(prop_underest = prop_underest,
             prop_overest = prop_overest)
    return(sim_data)
  }
  
  #function to simulate data for two groups
  one_run_data <- function(){
    sim_data1 <- visits_vl_df(beta2 = beta2a, n_participants = n_participants, n_visits = n_visits, delta = 0.5)
    sim_data1$group <- 1
    sim_data2 <- visits_vl_df(beta2 = beta2b, n_participants = n_participants, n_visits = n_visits, delta = 0.5)
    sim_data2$group <- 2
    sim_data <- rbind(sim_data1, sim_data2)
    sim_data <- add_imputation(sim_data)
    sim_data <- proportion_right_censored(sim_data)
    sim_data <- proportion_misclassified(sim_data)
    return(sim_data)
  }
  
  #function to fit cox model on true h crossing time in order to estimate true beta
  fit_true <- function(sim_data){
    cox_fit <- coxph(Surv(true_crossing_time, true_status) ~ group, data = sim_data)
    beta_hat_true <- coef(cox_fit)
    return(beta_hat_true)
  }
  
  #function for fitting and evaluating cox PH models on different types of data 
  fit_cox_ph <- function(sim_data){
    
    #fit cox model: right imputation on true (h) data
    cox_fit <- coxph(Surv(r_imp_true, true_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 1 #'r_imp_true'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result1 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: right imputation on observed (y) data
    cox_fit <- coxph(Surv(r_imp_obs, obs_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 2 #'r_imp_obs'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result2 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p) 
    
    #fit cox model: mid imputation on true (h) data
    cox_fit <- coxph(Surv(m_imp_true, true_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 3 #'m_imp_true'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result3 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: mid imputation on observed (y) data
    cox_fit <- coxph(Surv(m_imp_obs, obs_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 4 #'m_imp_obs'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result4 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: left imputation on true (h) data
    cox_fit <- coxph(Surv(l_imp_true, true_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 5 #'l_imp_true'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result5 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: left imputation on observed (y) data
    cox_fit <- coxph(Surv(l_imp_obs, obs_status) ~ group, data = sim_data)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    data_type <- 6 #'l_imp_obs'
    hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
    result6 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: interval censored on true (h) data
    cox_fit <- ic_sp(Surv(true_interval_L, true_interval_R, type = "interval2") ~ group, data = sim_data, bs_samples = 100)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    se_hat <- se_hat[1,1]
    data_type <- 7 #'ic_true'
    ses   <- sqrt(diag(vcov(cox_fit)))
    zvals <- beta_hat / ses
    hr_p <- 2 * (1 - pnorm(abs(zvals)))
    result7 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #fit cox model: interval censored on observed (y) data
    cox_fit <- ic_sp(Surv(obs_interval_L, obs_interval_R, type = "interval2") ~ group, data = sim_data, bs_samples = 100)
    beta_hat <- coef(cox_fit)
    se_hat <- sqrt(cox_fit$var)
    se_hat <- se_hat[1,1]
    data_type <- 8 #'ic_obs'
    ses   <- sqrt(diag(vcov(cox_fit)))
    zvals <- beta_hat / ses
    hr_p <- 2 * (1 - pnorm(abs(zvals)))
    result8 <- c(data_type = data_type, beta_hat = beta_hat, se_hat = se_hat, hr_p = hr_p)
    
    #bind the results
    result <- rbind(result1, result2, result3, result4, result5, result6, result7, result8)
    
    return(result)
  }
  
  #function for executing log rank test on different types of data
  do_log_rank <- function(sim_data){
    #test on right imputation of true (h) data
    test1 <- logrank_test(Surv(r_imp_true, true_status) ~ factor(group), data = sim_data)
    statistic1 <- statistic(test1)
    p1 <- pvalue(test1)
    lrt1 <- c(logrank_stat = statistic1, logrank_p = p1)
    #test on mid imputation of true (h) data
    test2 <- logrank_test(Surv(m_imp_true, true_status) ~ factor(group), data = sim_data)
    statistic2 <- statistic(test2)
    p2 <- pvalue(test2)
    lrt2 <- c(logrank_stat = statistic2, logrank_p = p2)
    #test on left imputation of true (h) data
    test3 <- logrank_test(Surv(l_imp_true, true_status) ~ factor(group), data = sim_data)
    statistic3 <- statistic(test3)
    p3 <- pvalue(test3)
    lrt3 <- c(logrank_stat = statistic3, logrank_p = p3)
    #test on interval censored true (h) data
    test4 <- ictest(Surv(true_interval_L, true_interval_R, type = "interval2") ~ group, data = sim_data, icontrol = icfitControl(maxit = 1000000))
    statistic4 <- test4$statistic
    p4 <- test4$p.value
    lrt4 <- c(logrank_stat = statistic4, logrank_p = p4)
    #test on right imputation of obs (y) data
    test5 <- logrank_test(Surv(r_imp_obs, obs_status) ~ factor(group), data = sim_data)
    statistic5 <- statistic(test5)
    p5 <- pvalue(test5)
    lrt5 <- c(logrank_stat = statistic5, logrank_p = p5)
    #test on mid imputation of obs (y) data
    test6 <- logrank_test(Surv(m_imp_obs, obs_status) ~ factor(group), data = sim_data)
    statistic6 <- statistic(test6)
    p6 <- pvalue(test6)
    lrt6 <- c(logrank_stat = statistic6, logrank_p = p6)
    #test on left imputation of obs (y) data
    test7 <- logrank_test(Surv(l_imp_obs, obs_status) ~ factor(group), data = sim_data)
    statistic7 <- statistic(test7)
    p7 <- pvalue(test7)
    lrt7 <- c(logrank_stat = statistic7, logrank_p = p7)
    #test on interval censored obs (y) data
    test8 <- ictest(Surv(obs_interval_L, obs_interval_R, type = "interval2") ~ group, data = sim_data, icontrol = icfitControl(maxit = 1000000))
    statistic8 <- test8$statistic
    p8 <- test8$p.value
    lrt8 <- c(logrank_stat = statistic8, logrank_p = p8)
    
    #bind results
    lrt <- rbind(lrt1, lrt5, lrt2, lrt6, lrt3, lrt7, lrt4, lrt8)
    
    return(lrt)
  }
  
  ####one run: simulate data, analysis, combine everything into a list
  sim_data <- one_run_data()
  
  #combine results of estimated logHR, cox models, and log rank tests for all runs
  one_run <- list(beta_hat_true = fit_true(sim_data), 
                  result = fit_cox_ph(sim_data),
                  lrt = do_log_rank(sim_data), 
                  prop_r_cens_true = sim_data$prop_r_cens_true[1:8],
                  prop_r_cens_obs = sim_data$prop_r_cens_obs[1:8],
                  prop_underest = sim_data$prop_underest[1:8],
                  prop_overest = sim_data$prop_overest[1:8],
                  study_structure = data.frame(beta2a = beta2a, beta2b = beta2b, visit_frequency = visit_frequency, n_visits = n_visits, participants_per_arm = n_participants, delta = delta)
  )
  
  return(one_run)
  
}
####ONE RUN FUNCTION ENDS HERE

####FOR LOOP  
#repeats one run of fitting all model types on simulated data and saves all runs

set.seed(9302002)
out_beta2_6.94_11.34_12visits_4weeks_25perarm <- list()
for (i in 1:1000){
  out_beta2_6.94_11.34_12visits_4weeks_25perarm[[i]] <- one_run_wang_comp_beta2(n_participants = 25, beta2a = 6.94, beta2b = 11.34, alpha = 0.05, n_visits = 12, delta = 0.5, visit_frequency = "4weeks", runs = 1)
}


