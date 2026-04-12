# VR Manuscript Simulation
## simulates time to viral rebound in an ATI trial according to Wang et al. (2020) model 
## fits the Cox PH model, log-rank test, and RMST to the data
## runs 1000 times for each simulation setting (one sample size/visit frequency/effect size combination)


# load required packages
library(survival)
library(tidyverse)
library(icenReg)
library(interval)
library(coin)
library(survRM2)

#-------------------------------------------------------------------------------
# Required arguments for VR_simulation:
# tau (numeric): restriction time for RMST
# n_participants (numeric): one-arm sample size for each simulated dataset.
# beta2a and beta2b (numeric): value for "delay in viral rebound" parameter from Wang et al. (2020) model
# alpha (numeric): significance level for hypothesis tests
# delta (numeric): one half of the allowable clinical visit window (in weeks) for each viral load measurement
# n_visits (numeric: 48, 24, or 12): number of viral load measurements over 48 weeks (proxy for the frequency of viral load measurements)
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Functions for data generation
#-------------------------------------------------------------------------------

#function that generates each visit time according to uniform dist (one participant)
#returns vector of visit times (in weeks)
visit_times <- function(n_visits, delta){
  visits <- vector("numeric", n_visits)
  if (n_visits == 48){
    for (i in 2:n_visits){
      visits[i] <- runif(n = 1, min = (i-1) - delta, max = (i-1) + delta)
    }
  }
  else if (n_visits == 24){
    for (i in 2:n_visits){
      visits[i] <- runif(n = 1, min = 2*(i-1) - delta, max = 2*(i-1) + delta)
    }
  }
  else if (n_visits == 12){
    for (i in 2:n_visits){
      visits[i] <- runif(n = 1, min = 4*(i-1) - delta, max = 4*(i-1) + delta)
    }
  }
  
  return(visits)
}

#function to calculate h (true) viral load at visit time, t
#returns a single viral load
h <- function(t, theta1, theta2, theta3, theta4, theta5){
  hvl <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t)))
  return(hvl)
}

#function to find true time of threshold crossing (threshold is 1000 copies/mL (3 log10(copies/mL)))
#returns the (approximate) time at which h = 3
find_true_time_of_crossing <- function(visits, n_visits, theta1, theta2, theta3, theta4, theta5){
  times <- seq(0, visits[n_visits], 0.01)
  h_of_times <- h(t = times, theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4, theta5 = theta5)
  check <- h_of_times > 3
  first_true <- which(check)[1]
  time_of_crossing <- times[first_true]
  return(time_of_crossing)
}

#function calculate y (observed) viral load (h + error) at visit time, t
#returns a single viral load
y <- function(t, theta1, theta2, theta3, theta4, theta5){
  yvl <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) +
    rnorm(1, mean = 0, sd = sqrt(0.26))
  return(yvl)
}

#function to calculate observed time of threshold crossing (threshold is 1000 copies/mL (3 log10(copies/mL)))
#returns the (approximate) time at which y = 3
find_obs_time_of_crossing <- function(visits, n_visits, theta1, theta2, theta3, theta4, theta5){
  times <- seq(0, visits[n_visits], 0.01)
  y_of_times <- y(t = times, theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4, theta5 = theta5)
  check <- y_of_times > 3
  first_true <- which(check)[1]
  time_of_crossing <- times[first_true]
  return(time_of_crossing)
}

#function to record L bound of true interval of crossing, assuming true time of crossing is captured by [L, R]. (L and R are in the set of a participant's visit times)
#returns L
true_interval_L <- function(visits, true_time, n_visits){
  check_L <-  visits > true_time #checks what visits are past the true time of crossing
  L_index <- which(check_L)[1] #checks which visit is the first visit past the true time of crossing
  true_L <- visits[L_index - 1] #the true L interval is the visit prior
  true_R <- visits[L_index] #the true R interval is the first visit past the true time of crossing
  if (true_time == visits[n_visits]){ #if the true time of crossing is right censored (occurs past the last visit)
    true_L <- visits[n_visits]
  }
  return(true_L)
}

#function to record R bound of true interval of crossing, assuming true time of crossing is captured by [L, R]. (L and R are in the set of a participant's visit times)
#returns R
true_interval_R <- function(visits, true_time, n_visits){
  check_L <-  visits > true_time
  L_index <- which(check_L)[1]
  true_L <- visits[L_index - 1]
  true_R <- visits[L_index]
  if (true_time == visits[n_visits]){
    true_R <- Inf
  }
  return(true_R)
}

#function to record L bound of observed interval of crossing, assuming observed time of crossing is captured by [L, R]. (L and R are in the set of a participant's visit times)
#returns L
observed_interval_L <- function(visits, yloads, n_visits){
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

#function to record R bound of observed interval of crossing, assuming observed time of crossing is captured by [L, R]. (L and R are in the set of a participant's visit times)
#returns R
observed_interval_R <- function(visits, yloads, n_visits){
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

#function to find true time to viral rebound survival status
# if the time of viral rebound according to h is right-censored, then func. returns status = 0. If it is not, then status = 1. 
find_true_status <- function(true_int_R){
  if (true_int_R == Inf){
    true_status <- 0
  }
  else {
    true_status <- 1
  }
  return(true_status)
}

#function to find observed time to viral rebound survival status
# if the time of viral rebound according to y is right-censored, then func. returns status = 0. If it is not, then status = 1. 
find_obs_status <- function(obs_int_R){
  if (obs_int_R == Inf){
    obs_status <- 0
  }
  else {
    obs_status <- 1
  }
  return(obs_status)
}

#function to create vector of visit times and corresponding viral loads based on model from Wang et al. 2020
#returns vector of visit times, viral loads produced by h at all visit times, viral loads produced by y at all visit times, 
# time to viral rebound using h, time to viral rebound using y, [L, R] for both measures of time to viral rebound, survival status for both measures of viral rebound
one_virus <- function(n_visits, beta2, delta){
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
  visits <- visit_times(n_visits = n_visits, delta = delta)
  for (i in seq_along(visits)){
    t <- visits[i]
    hload_t <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) #calculates h(), which is the viral load trajectory
    hloads <- c(hloads, hload_t)
    yload_t <- (theta1)*(t/(t+exp(theta2-theta3*t))) + theta4*(1/(1+exp(theta5*t))) + #calculates y(), which is h() + random error
      rnorm(1, mean = 0, sd = sqrt(0.26))
    yloads <- c(yloads, yload_t)
  }
  
  true_time <- find_true_time_of_crossing(visits = visits, n_visits = n_visits, theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4, theta5 = theta5)
  if (is.na(true_time)){
    true_time <- visits[n_visits]
  }
  obs_time <- find_obs_time_of_crossing(visits = visits, n_visits = n_visits, theta1 = theta1, theta2 = theta2, theta3 = theta3, theta4 = theta4, theta5 = theta5)
  if (is.na(obs_time)){
    obs_time <- visits[n_visits]
  }
  
  true_int_L <- true_interval_L(visits = visits, true_time = true_time, n_visits = n_visits)
  true_int_R <- true_interval_R(visits = visits, true_time = true_time, n_visits = n_visits)
  obs_int_L <- observed_interval_L(visits = visits, yloads = yloads, n_visits = n_visits)
  obs_int_R <- observed_interval_R(visits = visits, yloads = yloads, n_visits = n_visits)
  
  true_status <- find_true_status(true_int_R = true_int_R)
  obs_status <- find_obs_status(obs_int_R = obs_int_R)
  
  #fuses everything into one vector
  visits_loads <- c(beta2, visits, hloads, yloads, true_time, obs_time, 
                    true_int_L, true_int_R, obs_int_L, obs_int_R, true_status, obs_status)
  return(visits_loads)
}

#function to simulate data for one trial arm
#returns sim_data, a dataframe containing the information returned by one_virus for one trial arm. 
visits_vl_df <- function(n_visits, delta, n_participants, beta2){
  
  #creates df of visits and viral loads for n_participants. makes sure that n_participants are included in each dataset (eliminates and replaces cases where obs viral load is over threshold at time 0)
  virus_df <- c()
  difference <- n_participants
  while (difference > 0){
    for (i in 1:difference){
      indiv <- one_virus(n_visits = n_visits, beta2 = beta2, delta = delta)
      if (indiv[7] == -1 | indiv[8] == -1 | indiv[9] == -1 | indiv[10] == -1){
        next
      }
      virus_df <- data.frame(rbind(virus_df, indiv))
    }
    difference <- n_participants - nrow(virus_df)
  }
  
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
  
  sim_data <- virus_df
  
  #returns final df of visit times and viral loads for each participant
  return(sim_data)
  
}

#function for imputing (L, M and R-imputation) the data
#returns sim_data, with left, midpoint, and right imputed times to viral rebound (for both h and y) added for all participants
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
#returns sim_data, with the proportion of right-censored observations for data generated by both h and y added.
proportion_right_censored <- function(sim_data){
  prop_r_cens_true <- mean(sim_data$true_status == 0)
  prop_r_cens_obs <- mean(sim_data$obs_status == 0)
  sim_data <- sim_data %>%
    mutate(prop_r_cens_true = prop_r_cens_true,
           prop_r_cens_obs = prop_r_cens_obs)
  return(sim_data)
}

#function for calculating the proportion of misclassified intervals in L and R directions
#returns sim_data, with the addition of:
# 1) the proportion of times the interval containing the time to viral rebound according to y is HIGHER than the interval containing the time to viral rebound by h (prop_overest), and
# 2) the proportion of times the interval containing the time to viral rebound according to y is LOWER than the interval containing the time to viral rebound by h (prop_underest)
proportion_misclassified <- function(sim_data){
  prop_underest <- mean(sim_data$obs_interval_L < sim_data$true_interval_L)
  prop_overest <- mean(sim_data$obs_interval_L > sim_data$true_interval_L)
  sim_data <- sim_data %>%
    mutate(prop_underest = prop_underest,
           prop_overest = prop_overest)
  return(sim_data)
}

#function to simulate data for two groups
#returns sim_data, which now has the visit times, viral loads by h and y at each of the visit times, the TTVR according to h and y, the survival status according to h and y,
# the interval [L, R] containing the TTVR according to h and y, the right, left, and midpoint-imputed TTVR according to h and y, the proportion of the dataset with misclassified intervals of TTVR, and the proportion of the dataset with right-censored TTVR
# for each subject in the two-arm ATI trial. 
# sim_data also includes if the TTVR according to h is less than 8 for each subject. 
one_run_data <- function(beta2a, beta2b, n_participants, n_visits, delta, visit_frequency){
  sim_data1 <- visits_vl_df(beta2 = beta2a, n_participants = n_participants, n_visits = n_visits, delta = delta)
  sim_data1$group <- 1
  sim_data2 <- visits_vl_df(beta2 = beta2b, n_participants = n_participants, n_visits = n_visits, delta = delta)
  sim_data2$group <- 2
  sim_data <- rbind(sim_data1, sim_data2)
  sim_data <- add_imputation(sim_data)
  sim_data <- proportion_right_censored(sim_data)
  sim_data <- proportion_misclassified(sim_data)
  sim_data <- sim_data %>%
    mutate(VRunder8 = as.numeric(true_crossing_time < 8))
  return(sim_data)
}


#-------------------------------------------------------------------------------
# Function to test the proportional hazards assumption for each generated dataset (sim_data)
#-------------------------------------------------------------------------------

#hypothesis test for proportional hazards (fits cox model on true_crossing_time)
test_prop_hazards <- function(df){
  
  #truth
  cox_fit <- coxph(Surv(true_crossing_time, true_status) ~ group, data = df)
  schoenfeld <- cox.zph(cox_fit)
  prop_hazards_p <- schoenfeld$table["group", "p"]
  
  return(prop_hazards_p)
  
}

#-------------------------------------------------------------------------------
# Functions for the RMST
#-------------------------------------------------------------------------------

#function for creating a matrix for curve_area function (defined below), restricted by tau.
#returns int_prob_matrix_restrict: matrix with four rows, each column is an interval on the turnbull curve for the simulated data set.
#row 1: left side of interval
#row 2: right side of interval
#row 3: drop in survival probability (turnbull intervals only)
#row 4: survival probability (non-turnbull intervals only)
int_matrix <- function(df, tau){
  fit <- ic_np(cbind(df$left, df$right), data = df)
  T_bull_Intervals <- fit$T_bull_Intervals
  p_hat <- fit$p_hat
  #check if last Turnbull interval is equal to infinity
  if (T_bull_Intervals[2,ncol(T_bull_Intervals)] == Inf){
    last_int_inf <- 1
  } else {
    last_int_inf <- 0
  }
  
  #find largest finite time of observation in data set, t_max
  times <- c(df$left, df$right)
  t_max <- max(times[is.finite(times)])
  
  #check if the last turnbull interval is half infinite. if it is, check if tau > t_max. if so, return error message.
  if (last_int_inf == 1){
    if (tau > t_max){
      stop("tau, the time of restriction, must be less than or equal to the largest time of observation.")
    }
  }
  
  #check if tau < 0. if not, return error message.
  if (tau <= 0){
    stop("tau, the time of restriction, must be greater than 0.")
  }
  
  #bind turnbull intervals and p_hat
  #delete p_hat values equal to 0 and corresponding turnbull intervals from T_bull_Intervals
  #first row is left side of turnbull intervals, second row is right side, third row is corresponding p_hat values
  int_prob_matrix <- rbind(T_bull_Intervals, p_hat)
  int_prob_matrix <- int_prob_matrix[, int_prob_matrix[3,] != 0, drop = FALSE]
  #add fourth row corresponding to survival function
  int_prob_matrix <- rbind(int_prob_matrix, rep(NA, ncol(int_prob_matrix)))
  #add non-turnbull interval columns
  new_cols <- c(0,int_prob_matrix[1,1], NA, 1)
  if (ncol(int_prob_matrix) != 1){
    for (i in 1:(ncol(int_prob_matrix)-1)){
      if (!isTRUE(all.equal(int_prob_matrix[2,i], int_prob_matrix[1, i+1]))){
        new_col <- c(int_prob_matrix[2,i], int_prob_matrix[1, i+1], NA, NA)
        new_cols <- cbind(new_cols, new_col)
      }
      else if (all.equal(int_prob_matrix[2,i], int_prob_matrix[1, i+1])){
        next
      }
    }
  } else if (ncol(int_prob_matrix) == 1){
    new_cols <- as.matrix(new_cols)
  }
  #fill in survival function (4th row of int_prob_matrix_full_sort) for non-turnbull intervals (new_cols)
  ps <- c()
  for (i in 1:ncol(new_cols)){
    ps[i] <- sum(int_prob_matrix[3, int_prob_matrix[1, ] < new_cols[1,i]])
  }
  
  new_cols[4,] <- 1 - ps
  
  #combine turnbull and regular intervals
  
  int_prob_matrix_full <- cbind(int_prob_matrix, new_cols)
  
  #sort them in time order, add final interval column if last turnbull interval doesn't extend to t_max
  
  int_prob_matrix_full_sort <- int_prob_matrix_full[, order(int_prob_matrix_full[1, ])]
  if (int_prob_matrix_full_sort[2, ncol(int_prob_matrix_full_sort)] < tau){
    int_prob_matrix_full_sort <- cbind(int_prob_matrix_full_sort, c(int_prob_matrix_full_sort[2,ncol(int_prob_matrix_full_sort)], tau, NA, 0))
  }
  
  #restrict int_prob_matrix_full_sort based on tau
  int_prob_matrix_restrict <- NULL
  for (i in (1:ncol(int_prob_matrix_full_sort))){
    if (tau > int_prob_matrix_full_sort[1,i] && tau <= int_prob_matrix_full_sort[2,i]){
      int_prob_matrix_restrict <- int_prob_matrix_full_sort[,1:i]
      break
    }
  }
  
  if (is.null(int_prob_matrix_restrict)) {
    print(int_prob_matrix_full_sort)
    stop("tau does not fall within any interval in int_prob_matrix_full_sort.")
  }
  
  int_prob_matrix_restrict <- as.matrix(int_prob_matrix_restrict)
  
  #check if tau is in a turnbull or non-turnbull interval
  #if tau is in a non-turnbull interval, change the right side of the interval to tau.
  #if tau is in a turnbull interval, find the proportion of the total interval that the interval to tau corresponds to.
  #then, find that proportion of the p_hat row (third row of int_prob_matrix_full_sort, the difference in probability across the turnbull interval).
  #then, replace the p_hat entry for that column with the new p_hat, and replace the right side of the interval with tau. 
  if (is.na(int_prob_matrix_restrict[3, ncol(int_prob_matrix_restrict)])){
    int_prob_matrix_restrict[2,ncol(int_prob_matrix_restrict)] <- tau
  } else if (!is.na(int_prob_matrix_restrict[3, ncol(int_prob_matrix_restrict)])){
    int_prop <- (tau - int_prob_matrix_restrict[1,ncol(int_prob_matrix_restrict)])/(int_prob_matrix_restrict[2,ncol(int_prob_matrix_restrict)] - int_prob_matrix_restrict[1,ncol(int_prob_matrix_restrict)])
    p_hat_new <- int_prop * int_prob_matrix_restrict[3,ncol(int_prob_matrix_restrict)]
    int_prob_matrix_restrict[3,ncol(int_prob_matrix_restrict)] <- p_hat_new
    int_prob_matrix_restrict[2,ncol(int_prob_matrix_restrict)] <- tau
  }
  
  return(int_prob_matrix_restrict)
  
}

#function for calculating full area under the turnbull curve (using linear smoothing technique) given int_prob_matrix_restrict
#returns area under the curve
curve_area <- function(int_prob_matrix_restrict){
  
  #RECTANGLES
  
  #vector of areas of rectangles under survival curve (width of non-turnbull interval * survival probability on this interval)
  rect_areas <- c()
  for (i in 1:ncol(int_prob_matrix_restrict)){
    if (!is.na(int_prob_matrix_restrict[4,i])){
      height <- int_prob_matrix_restrict[4,i]
      width <- int_prob_matrix_restrict[2,i] - int_prob_matrix_restrict[1,i]
      rect_areas <- c(rect_areas, height * width)
      
    } 
  }
  
  #TRAPEZOIDS (area = 0.5*(b1 + b2)*width)
  #calculate b1, b2, and width for all intervals (turnbull and non-turnbull)
  #apply the trapezoid area formula for all intervals
  #delete the indices of the trap_areas vector that correspond to non-turnbull intervals (because these are rectangles)
  int_prob_matrix_restrict[is.na(int_prob_matrix_restrict)] <- 0
  width <- c()
  if (ncol(int_prob_matrix_restrict) == 1){
    b1 <- c(1)
    if (int_prob_matrix_restrict[3,1] == 0){
      b2 <- c(1)
    } else {
      b2 <- c(1-int_prob_matrix_restrict[3,1])
    }
  } else {
    b1 <- c(1)
    b2 <- 1 - int_prob_matrix_restrict[3,1]
    for (i in 2:ncol(int_prob_matrix_restrict)){
      b1 <- c(b1, b1[length(b1)] - int_prob_matrix_restrict[3,i-1])
      b2 <- c(b2, b1[length(b1)] - int_prob_matrix_restrict[3, i])
    }
  }
  
  for (i in 1:ncol(int_prob_matrix_restrict)){
    width <- c(width, int_prob_matrix_restrict[2, i] - int_prob_matrix_restrict[1,i])
  }
  trap_areas <- 0.5*(b1 + b2)*width
  is <- c()
  for (i in 1:ncol(int_prob_matrix_restrict)){
    if (int_prob_matrix_restrict[3,i] == 0){
      is <- c(is, i)
    }
  }
  trap_areas <- trap_areas[-is]
  
  
  #TOTAL AREA UNDER CURVE
  
  area_under_curve <- sum(rect_areas, trap_areas)
  return(area_under_curve)
}

#function for calculating RMST for one trial arm with interval-censored data.
#returns RMST
IC_RMST <- function(df, tau){
  
  
  int_prob_matrix_restrict <- int_matrix(df = df, tau = tau)
  if (is.null(int_prob_matrix_restrict)) {
    stop("int_matrix() failed to return a matrix.")
  }
  rmst <- curve_area(int_prob_matrix_restrict = int_prob_matrix_restrict)
  return(rmst)
  
}    

#function for hypothesis test for difference in two RMSTs with interval-censored data
## 1. calculate the difference in RMSTs between two arms
## 2. construct percentile bootstrap 95% confidence interval (sigma is the boostrap standard error)
## 3. reject the null if confidence interval does not contain 0
#returns dataframe with RMST of both groups, difference in RMSTs, standard error of difference in RMSTs, 95% confidence interval bounds of difference in RMSTs, and if confidence interval contains 0 (saves 1 if it does, 0 if not)
diff_IC_RMST <- function(df, tau){
  df1 <- df[df$group == 0,]
  df2 <- df[df$group == 1,]
  RMST1 <- IC_RMST(df = df1, tau = tau)
  RMST2 <- IC_RMST(df = df2, tau = tau)
  diff <- RMST2 - RMST1
  # bootstrap
  boot_diffs <- c()
  for (i in 1:1000){
    repeat{
      result <- tryCatch({
        df1_boot_sample <- df1[sample(1:nrow(df1), size = nrow(df1), replace = TRUE),]
        RMST1_boot <- IC_RMST(df = df1_boot_sample, tau = tau)
        break
      }, error = function(e){
        message("Error in IC_RMST(df = df1_boot_sample) at i = ", i, ": ", e$message)
        
      })
      
    }
    repeat{
      result <- tryCatch({
        df2_boot_sample <- df2[sample(1:nrow(df2), size = nrow(df2), replace = TRUE),]
        RMST2_boot <- IC_RMST(df = df2_boot_sample, tau = tau)
        break
      }, error = function(e){
        message("Error in IC_RMST(df = df2_boot_sample) at i = ",i, ": ", e$message)
      })
      
    }
    boot_diffs[i] <- RMST2_boot - RMST1_boot
  }
  boot_se <- sd(boot_diffs)
  lci <- unname(quantile(boot_diffs, probs=0.025))
  uci <- unname(quantile(boot_diffs, probs=0.975))
  check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
  RMST_df <- data.frame(RMST1 = RMST1,
                        RMST2 = RMST2,
                        diff_RMST = diff,
                        se_diff_RMST = boot_se,
                        lci_diff_RMST = lci,
                        uci_diff_RMST = uci,
                        covers_zero_CI_diff_RMST = check_cover_0)
  return(RMST_df)
}

#function to run RMST hypothesis tests for all 8 "data types" (R, M, L-imputed and interval-censored data from h and y) and save in one dataframe
#returns dataframe with results of RMST hypothesis test on one simulated dataset (RMST of both groups, difference in RMSTs, standard error of difference in RMSTs, 95% confidence interval bounds of difference in RMSTs, and if confidence interval contains 0 (saves 1 if it does, 0 if not))
# for all 8 data types and for the true TTVR.
#returned dataframe also includes the proportion of misclassifed intervals, right-censored observations, and the proportion of times the TTVR according to h < 8 for each group (all from sim_data). 
#note: the function rmst2() will throw an error if tau is greater than the largest observed event time. In this rare case, NA is saved for that set of RMST results. 
RMST_tests <- function(df, tau){
  
  #change group labeling for rmst2()
  df$group[df$group == 1] <- "0"
  df$group[df$group == 2] <- "1"
  
  #create dfs for IC_RMST function
  sim_data_IC_obs <- df[,c("group", "obs_interval_L", "obs_interval_R")]
  colnames(sim_data_IC_obs)[colnames(sim_data_IC_obs) == "obs_interval_L"] <- "left"
  colnames(sim_data_IC_obs)[colnames(sim_data_IC_obs) == "obs_interval_R"] <- "right"
  #cat("df for IC func obs:\n")
  #print(sim_data_IC_obs)
  sim_data_IC_true <- df[,c("group", "true_interval_L", "true_interval_R")]
  colnames(sim_data_IC_true)[colnames(sim_data_IC_true) == "true_interval_L"] <- "left"
  colnames(sim_data_IC_true)[colnames(sim_data_IC_true) == "true_interval_R"] <- "right"
  #cat("df for IC func true:\n")
  #print(sim_data_IC_true)
  
  # true (exact) crossing times (for bias calculation)
  truth <- tryCatch(
    rmst2(time = df$true_crossing_time, status = df$true_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(truth)){
    truth_RMST_df <- data.frame(RMST1 = NA,
                                RMST2 = NA,
                                diff_RMST = NA,
                                se_diff_RMST = NA,
                                lci_diff_RMST = NA,
                                uci_diff_RMST = NA,
                                covers_zero_CI_diff_RMST = NA)
  } else{
    RMST1 <- truth$RMST.arm0$rmst[1]
    RMST2 <- truth$RMST.arm1$rmst[1]
    true_diff <- truth$unadjusted.result[1]
    se <- sqrt((truth$RMST.arm1$rmst[2])^2 + (truth$RMST.arm0$rmst[2])^2)
    lci <- truth$unadjusted.result[4]
    uci <- truth$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    truth_RMST_df <- data.frame(RMST1 = RMST1,
                                RMST2 = RMST2,
                                diff_RMST = true_diff,
                                se_diff_RMST = se,
                                lci_diff_RMST = lci,
                                uci_diff_RMST = uci,
                                covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  
  # R imputed, true
  #print(df$r_imp_true)
  r_true <- tryCatch(
    rmst2(time = df$r_imp_true, status = df$true_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(r_true)){
    r_true_RMST_df <- data.frame(RMST1 = NA,
                                 RMST2 = NA,
                                 diff_RMST = NA,
                                 se_diff_RMST = NA,
                                 lci_diff_RMST = NA,
                                 uci_diff_RMST = NA,
                                 covers_zero_CI_diff_RMST = NA)
  } else{
    RMST1 <- r_true$RMST.arm0$rmst[1]
    RMST2 <- r_true$RMST.arm1$rmst[1]
    diff <- r_true$unadjusted.result[1]
    se <- sqrt((r_true$RMST.arm1$rmst[2])^2 + (r_true$RMST.arm0$rmst[2])^2)
    lci <- r_true$unadjusted.result[4]
    uci <- r_true$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    r_true_RMST_df <- data.frame(RMST1 = RMST1,
                                 RMST2 = RMST2,
                                 diff_RMST = diff,
                                 se_diff_RMST = se,
                                 lci_diff_RMST = lci,
                                 uci_diff_RMST = uci,
                                 covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  # R imputed, obs
  
  r_obs <- tryCatch(
    rmst2(time = df$r_imp_obs, status = df$obs_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(r_obs)){
    r_obs_RMST_df <- data.frame(RMST1 = NA,
                                RMST2 = NA,
                                diff_RMST = NA,
                                se_diff_RMST = NA,
                                lci_diff_RMST = NA,
                                uci_diff_RMST = NA,
                                covers_zero_CI_diff_RMST = NA)
  } else{
    RMST1 <- r_obs$RMST.arm0$rmst[1]
    RMST2 <- r_obs$RMST.arm1$rmst[1]
    diff <- r_obs$unadjusted.result[1]
    se <- sqrt((r_obs$RMST.arm1$rmst[2])^2 + (r_obs$RMST.arm0$rmst[2])^2)
    lci <- r_obs$unadjusted.result[4]
    uci <- r_obs$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    r_obs_RMST_df <- data.frame(RMST1 = RMST1,
                                RMST2 = RMST2,
                                diff_RMST = diff,
                                se_diff_RMST = se,
                                lci_diff_RMST = lci,
                                uci_diff_RMST = uci,
                                covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  # M imputed, true
  #print(df$m_imp_true)
  m_true <- tryCatch(
    rmst2(time = df$m_imp_true, status = df$true_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(m_true)){
    m_true_RMST_df <- data.frame(RMST1 = NA,
                                 RMST2 = NA,
                                 diff_RMST = NA,
                                 se_diff_RMST = NA,
                                 lci_diff_RMST = NA,
                                 uci_diff_RMST = NA,
                                 covers_zero_CI_diff_RMST = NA)
  } else {
    RMST1 <- m_true$RMST.arm0$rmst[1]
    RMST2 <- m_true$RMST.arm1$rmst[1]
    diff <- m_true$unadjusted.result[1]
    se <- sqrt((m_true$RMST.arm1$rmst[2])^2 + (m_true$RMST.arm0$rmst[2])^2)
    lci <- m_true$unadjusted.result[4]
    uci <- m_true$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    m_true_RMST_df <- data.frame(RMST1 = RMST1,
                                 RMST2 = RMST2,
                                 diff_RMST = diff,
                                 se_diff_RMST = se,
                                 lci_diff_RMST = lci,
                                 uci_diff_RMST = uci,
                                 covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  # M imputed, obs
  #print(df$m_imp_obs)
  m_obs <- tryCatch(
    rmst2(time = df$m_imp_obs, status = df$obs_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(m_obs)){
    m_obs_RMST_df <- data.frame(RMST1 = NA,
                                RMST2 = NA,
                                diff_RMST = NA,
                                se_diff_RMST = NA,
                                lci_diff_RMST = NA,
                                uci_diff_RMST = NA,
                                covers_zero_CI_diff_RMST = NA)
  } else {
    RMST1 <- m_obs$RMST.arm0$rmst[1]
    RMST2 <- m_obs$RMST.arm1$rmst[1]
    diff <- m_obs$unadjusted.result[1]
    se <- sqrt((m_obs$RMST.arm1$rmst[2])^2 + (m_obs$RMST.arm0$rmst[2])^2)
    lci <- m_obs$unadjusted.result[4]
    uci <- m_obs$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    m_obs_RMST_df <- data.frame(RMST1 = RMST1,
                                RMST2 = RMST2,
                                diff_RMST = diff,
                                se_diff_RMST = se,
                                lci_diff_RMST = lci,
                                uci_diff_RMST = uci,
                                covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  # L imputed, true
  
  l_true <- tryCatch(
    rmst2(time = df$l_imp_true, status = df$true_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(l_true)){
    l_true_RMST_df <- data.frame(RMST1 = NA,
                                 RMST2 = NA,
                                 diff_RMST = NA,
                                 se_diff_RMST = NA,
                                 lci_diff_RMST = NA,
                                 uci_diff_RMST = NA,
                                 covers_zero_CI_diff_RMST = NA)
  } else {
    RMST1 <- l_true$RMST.arm0$rmst[1]
    RMST2 <- l_true$RMST.arm1$rmst[1]
    diff <- l_true$unadjusted.result[1]
    se <- sqrt((l_true$RMST.arm1$rmst[2])^2 + (l_true$RMST.arm0$rmst[2])^2)
    lci <- l_true$unadjusted.result[4]
    uci <- l_true$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    l_true_RMST_df <- data.frame(RMST1 = RMST1,
                                 RMST2 = RMST2,
                                 diff_RMST = diff,
                                 se_diff_RMST = se,
                                 lci_diff_RMST = lci,
                                 uci_diff_RMST = uci,
                                 covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  #L imputed, obs
  
  l_obs <- tryCatch(
    rmst2(time = df$l_imp_obs, status = df$obs_status, arm = df$group, tau = tau),
    error = function(e) {
      cat("Error:", e$message, "\n")
      NULL
    }
  )
  if (is.null(l_obs)){
    l_obs_RMST_df <- data.frame(RMST1 = NA,
                                RMST2 = NA,
                                diff_RMST = NA,
                                se_diff_RMST = NA,
                                lci_diff_RMST = NA,
                                uci_diff_RMST = NA,
                                covers_zero_CI_diff_RMST = NA)
  } else {
    RMST1 <- l_obs$RMST.arm0$rmst[1]
    RMST2 <- l_obs$RMST.arm1$rmst[1]
    diff <- l_obs$unadjusted.result[1]
    se <- sqrt((l_obs$RMST.arm1$rmst[2])^2 + (l_obs$RMST.arm0$rmst[2])^2)
    lci <- l_obs$unadjusted.result[4]
    uci <- l_obs$unadjusted.result[7]
    check_cover_0 <- as.numeric((lci < 0) & (uci > 0))
    l_obs_RMST_df <- data.frame(RMST1 = RMST1,
                                RMST2 = RMST2,
                                diff_RMST = diff,
                                se_diff_RMST = se,
                                lci_diff_RMST = lci,
                                uci_diff_RMST = uci,
                                covers_zero_CI_diff_RMST = check_cover_0)
  }
  
  # IC data, true
  IC_true_RMST_df <- diff_IC_RMST(df = sim_data_IC_true, tau = tau)
  # IC data, observed
  IC_obs_RMST_df <- diff_IC_RMST(df = sim_data_IC_obs, tau = tau)
  
  full_RMST_results <- rbind(truth_RMST_df, r_true_RMST_df, r_obs_RMST_df, m_true_RMST_df, 
                             m_obs_RMST_df, l_true_RMST_df, l_obs_RMST_df, 
                             IC_true_RMST_df, IC_obs_RMST_df)
  
  full_RMST_results <- cbind(full_RMST_results, df$prop_r_cens_true[1:9], df$prop_r_cens_obs[1:9],
                             df$prop_underest[1:9], df$prop_overest[1:9])
  
  #calculate the proportion of rebound times < 8 in both groups
  propVRunder8_arm0 <- mean(df$VRunder8[df$group == "0"])
  propVRunder8_arm0 <- rep(propVRunder8_arm0, times = nrow(full_RMST_results))
  propVRunder8_arm1 <- mean(df$VRunder8[df$group == "1"])
  propVRunder8_arm1 <- rep(propVRunder8_arm1, times = nrow(full_RMST_results))
  
  full_RMST_results <- full_RMST_results %>%
    mutate(propVRunder8_arm0 = propVRunder8_arm0,
           propVRunder8_arm1 = propVRunder8_arm1)
  
  
  return(full_RMST_results)
}

#-------------------------------------------------------------------------------
# Functions for the Cox PH model
#-------------------------------------------------------------------------------

#function for fitting Cox PH models on true TTVR (generated by h) and imputed and interval-censored data generated by h and y (8 "data types"-- not including the true TTVR)
#returns dataframe with log hazard ratio estimate, esimated standard error, and the p-value associated with the log HR for all 8 data types
fit_cox_ph <- function(sim_data){
  
  #fit cox model: true TTVR (generated by h)
  cox_fit <- coxph(Surv(true_crossing_time, true_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 0 #'truth'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result0 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: right imputation on true (h) data
  cox_fit <- coxph(Surv(r_imp_true, true_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 1 #'r_imp_true'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result1 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: right imputation on observed (y) data
  cox_fit <- coxph(Surv(r_imp_obs, obs_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 2 #'r_imp_obs'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result2 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: mid imputation on true (h) data
  cox_fit <- coxph(Surv(m_imp_true, true_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 3 #'m_imp_true'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result3 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: mid imputation on observed (y) data
  cox_fit <- coxph(Surv(m_imp_obs, obs_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 4 #'m_imp_obs'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result4 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: left imputation on true (h) data
  cox_fit <- coxph(Surv(l_imp_true, true_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 5 #'l_imp_true'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result5 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: left imputation on observed (y) data
  cox_fit <- coxph(Surv(l_imp_obs, obs_status) ~ group, data = sim_data)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  data_type <- 6 #'l_imp_obs'
  hr_p <- summary(cox_fit)$coefficients[, "Pr(>|z|)"]
  result6 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: interval censored on true (h) data
  cox_fit <- ic_sp(Surv(true_interval_L, true_interval_R, type = "interval2") ~ group, data = sim_data, bs_samples = 100)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  se_hat <- se_hat[1,1]
  data_type <- 7 #'ic_true'
  ses   <- sqrt(diag(vcov(cox_fit)))
  zvals <- beta_hat / ses
  hr_p <- 2 * (1 - pnorm(abs(zvals)))
  result7 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #fit cox model: interval censored on observed (y) data
  cox_fit <- ic_sp(Surv(obs_interval_L, obs_interval_R, type = "interval2") ~ group, data = sim_data, bs_samples = 100)
  beta_hat <- coef(cox_fit)
  se_hat <- sqrt(cox_fit$var)
  se_hat <- se_hat[1,1]
  data_type <- 8 #'ic_obs'
  ses   <- sqrt(diag(vcov(cox_fit)))
  zvals <- beta_hat / ses
  hr_p <- 2 * (1 - pnorm(abs(zvals)))
  result8 <- c(data_type = data_type, logHR_hat = beta_hat, logHR_se_hat = se_hat, logHR_p = hr_p)
  
  #bind the results
  result <- rbind(result0, result1, result2, result3, result4, result5, result6, result7, result8)
  
  return(result)
}

#-------------------------------------------------------------------------------
# Functions for the log-rank test
#-------------------------------------------------------------------------------

#function for executing log rank test on imputed and interval-censored data generated by h and y (8 "data types")
#returns dataframe with log-rank test statistic and p-value for all 8 data types
do_log_rank <- function(sim_data){
  
  #test on true (h) TTVR
  test0 <- logrank_test(Surv(true_crossing_time, true_status) ~ factor(group), data = sim_data)
  statistic0 <- statistic(test0)
  p0 <- pvalue(test0)
  lrt0 <- c(logrank_stat = statistic0, logrank_p = p0)
  
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
  lrt <- rbind(lrt0, lrt1, lrt5, lrt2, lrt6, lrt3, lrt7, lrt4, lrt8)
  
  return(lrt)
}

#-------------------------------------------------------------------------------
# One run of the full simulation
#-------------------------------------------------------------------------------

#function that returns proportional hazards test result and RMST, log-rank test, and Cox PH model results for one simulated dataset
VRSimulation_one_run <- function(tau, n_participants, beta2a, beta2b, alpha, n_visits, delta){
  
  #simulate data
  sim_data <- one_run_data(beta2a = beta2a, beta2b = beta2b, n_participants = n_participants, n_visits = n_visits, 
                           delta = delta)
  
  #prop hazard test on sim_data
  prop_haz <- test_prop_hazards(df = sim_data)
  
  #RMST results
  one_run_RMST <- RMST_tests(df = sim_data, tau = tau)
  
  #Cox PH results
  one_run_CPH <- fit_cox_ph(sim_data = sim_data)
  
  #log-rank results
  one_run_LRT <- do_log_rank(sim_data = sim_data)
  
  #combine RMST results with log-rank and Cox PH results
  one_run <- cbind(one_run_RMST, one_run_CPH, one_run_LRT)
  
  #combine one_run df with prop. hazards result
  one_run <- one_run %>%
    mutate(prop_haz_p = prop_haz)
  
  
  return(one_run)
  
}

#-------------------------------------------------------------------------------
# 1000 runs of the simulation: do this for each simulation setting of interest
#-------------------------------------------------------------------------------

sim_results <- list()
i <- 1
while (i <= 1000) {
  cat("Running iteration", i, "\n")
  
  sim_results[[i]] <- VRSimulation_one_run(tau = 8, n_participants = 100, beta2a = 6.94, beta2b = 11.34, alpha = 0.05, n_visits = 48, delta = 0.5)
  
  i <- i + 1
}





