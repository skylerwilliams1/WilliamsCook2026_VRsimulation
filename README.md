# WilliamsCook2026_VRsimulation

VRSimulation.R and VRSimulation_post_processing.R include the simulation code used to produce the results in the Williams and Cook viral rebound manuscript. Find the full explanation of our simulation in the publication. Below, we briefly describe how to run the simulation and post processing, although both scripts are commented with further detail and explanations.  Notably, our code to calculate the RMST and execute the test of the difference in RMST with interval-censored data is included in VRSimulation.R. We also offer a description for each variable in the final dataset that the simulation and post-processing returns. 

## VRSimulation.R

### Required arguments:

- tau: the restriction time (in weeks) for the RMST. We used tau = 8.
- n_participants: the one-arm sample size. We used 100 and 25.
- beta2a: the delay in viral rebound parameter for trial arm one. We used 6.94.
- beta2b: the delay in viral rebound parameter for trial arm two. We used 6.94, 9.14, and 11.34.
- alpha: the desired significance level for hypothesis tests. We used 0.05.
- delta: one half of the clinical visit window for viral load monitoring (in weeks). We used 0.5.
- n_visits: the number of uniformly distributed viral load monitoring visits over 48 weeks. Our function only takes 48, 24, and 12, but can be easily altered to take other values.

### Instructions:

- To run one simulation setting with 1000 runs, simply load all functions in order and run the while loop at the bottom. Specify the arguments as you wish.
- The code will return a list with 1000 items, where each item contains the result of 1 run (1 dataset generated and analyzed). Save this list and proceed to post-processing.

## VRSimulation_post_processing.R

### Required arguments:

The same as VRSimulation.R, in addition to:

- v: the list returned by VRSimulation.R
- runs: the number of runs used in VRSimulation.R. We used 1000.

### Instructions:

- Simply load and run the post_processing function with the arguments specified as appropriate.
- A data frame will be returned with the full results.

## Explanation of variables in returned dataset

- data_run: run number of 1000 run simulaton.
- data_type_description: shorthand to describe the data that was used for analysis (imputed or interval-censored, "true" or "observed" data).
- data_type: a numerical code for data_type_description.
- tau: the restriction time for the RMST. We used tau = 8.
- participants_per_arm: the one-arm sample size.
- beta2a: the delay in viral rebound parameter for trial arm one.
- beta2b: the delay in viral rebound parameter for trial arm two.
- alpha: the desired significance level for hypothesis tests. 
- delta: one half of the clinical visit window for viral load monitoring (in weeks).
- n_visits: the number of uniformly distributed viral load monitoring visits over 48 weeks. 
- RMST1/RMST2: the RMST of arm one/two.
- diff_RMST: RMST2 - RMST1
- se_diff_RMST: the estimated standard error of diff_RMST.
- RMST_mean_est_se. The mean of se_diff_RMST within each data type.
- RMST_emp_se: the standard deviation of diff_RMST within each data type. 
- lci_diff_RMST: the lower bound of the 95% confidence interval of true_RMST_diff.
- uci_diff_RMST: the upper bound of the 95% confidence interval of true_RMST_diff.
- true_RMST_diff: the mean RMST_diff within the "truth" data type: the exact (not imputed or interval-censored) time of viral rebound from "true" data. In other words, a best estimate of the true difference in RMSTs across the two arms of the trial.
- true_RMST1: the mean RMST1 within the "truth" data type.
- true_RMST2: the mean RMST2 within the "truth" data type.
- covers_zero_CI_diff_RMST: 0 or 1, where 1 indicates that the 95% confidence interval of true_RMST_diff covers 0.
- check_cover_true_RMST_diff: 0 or 1, where 1 indicates that the 95% confidence interval of true_RMST_diff covers the true_RMST_diff.
- RMST_CI_covers_0_prob: the proportion of runs where covers_zero_CI_diff_RMST = 1 within each data type.
- RMST_CI_covers_true_RMST_diff_prob: the proportion of runs where check_cover_true_RMST_diff = 1 within each data type. The coverage probability of the diff_RMST 95% confidence interval.
- RMST_prop_null_rejected: 1 - RMST_CI_covers_0_prob. The proportion of runs where the null hypothesis that the two RMSTs are equal between groups is rejected. 
- bias_RMST_diff: the mean bias of the RMST_diff within each data_type.
- RMST_bias_arm1: the bias of RMST1 within each data type.
- RMST_bias_arm2: the bias of RMST2 within each data type.
- RMST_n_iters_na: the number of runs where tau > the largest recorded event time. RMST analysis was not completed in these runs. 
- prop_r_cens_true: the proportion of right-censored observations in "true" data.
- prop_r_cens_obs: the proportion of right-censored observations in "observed" data.
- prop_underest: the proportion of underestimated intervals in "observed" data.
- prop_overest: the proportion of overestimated intervals in "observed" data.
- propVRunder8_arm0: the proportion of events that occurred under 8 weeks in the first arm of the trial.
- propVRunder8_arm1: the proportion of events that occurred under 8 weeks in the second arm of the trial.
- mean_propVRunder8_arm0: the mean proportion of events that occurred under 8 weeks in the first arm of the trial.
- mean_propVRunder8_arm1: the mean proportion of events that occurred under 8 weeks in the second arm of the trial.
- logHR_hat: estimated log hazard ratio from Cox proportional hazards model.
- logHR_se_hat: estimated standard error of logHR_hat.
- logHR_p: p-value of logHR_hat.
- null_rejected_logHR_check: 0 or 1, where 1 indicates logHR_p < 0.05.
- logHR_prop_null_rejected: the proportion of runs where null_rejected_logHR_check = 1 within each data type. 
- true_logHR: the mean of the logHR_hat within the "truth" data type: the exact (not imputed or interval-censored) time of viral rebound from "true" data. In other words, a best estimate of the true log hazard ratio across the two arms of the trial.
- logHR_lci: the lower bound of the 95% confidence interval of true_logHR.
- logHR_uci: the upper bound of the 95% confidence interval of true_logHR.
- logHR_CI_check_cover: 0 or 1, where 1 indicates that the 95% confidence interval of true_logHR covers true_logHR.
- logHR_covers_true_logHR_prob: the proportion of runs where logHR_CI_check_cover = 1. The coverage probability of the logHR_true 95% confidence intervals.
- bias_logHR: the mean bias of logHR_hat within each data type.
- logHR_emp_se: the standard deviation of logHR_hat within each data type.
- logHR_mean_est_se: the mean estimated standard of logHR_hat within each data type. 
- logrank_stat: log-rank test statistic.
- logrank_p: corresponding p-value of logrank_stat.
- logrank_p_check: 0 or 1, where 1 indicates logrank_p < 0.05.
- LRT_prop_null_rejected: the proportion of runs where logrank_p_check = 1 within each data type.
- prop_haz_p: p-value of test for proportional hazards (Schoenfeld residuals test)
- ph_violated: 0 or 1, where 1 indicates that prop_haz_p < 0.05.
- prop_ph_violated: the proportion of runs in the entire simulation in which the proportional hazards assumption was violated.





