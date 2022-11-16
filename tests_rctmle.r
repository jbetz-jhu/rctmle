rm(list = ls())

set.seed(12345)

source(
  file.path(
    # "https://raw.githubusercontent.com/jbetz-jhu/rctmle/main/rctmle.r"
    "C:", "Users", "josh-jhbc", "Documents", "rctmle", "rctmle.r"
  )
)


n_per_arm <- 150
n_covariates <- 5
n_non_monotone <- 15
visit_times <- c(0, 0.5, 1)
visits <- length(visit_times)


### Create Data - Null Treatment ###############################################
sim_data <-
  sim_lme_trial(
    bl_covariates =
      iid_centered_mvn(
        n = n_per_arm*2,
        n_covariates = n_covariates,
        sigma_x = diag(n_covariates)
      ),
    visit_times = visit_times,
    mean_outcomes = rep(10, visits),
    outcome_cov = 
      list(
        c("x_1" = 1, "x_2" = 0.5, "x_3" = -0.25),
        c("x_1" = 1, "x_2" = 0.5, "x_3" = -0.25),
        c("x_1" = 1, "x_2" = 0.5, "x_3" = -0.25)
      ),
    re_variance = c(0.5),
    re_correlation = NULL,
    residual_sd = 1,
    pr_dropout = rep(0.1, visits),
    dropout_cov = NULL
  )



# Introduce non-monotone missingness
non_monotone_rows <- which(sim_data$last_observed > 1)
sim_data$y_obs_1[sample(x = non_monotone_rows, size = n_non_monotone)] <- NA
sim_data$y_obs_2[sample(x = non_monotone_rows, size = n_non_monotone)] <- NA


# Map the outcome to (0, 1) interval
y_max <- ceiling(max(unlist(sim_data[, paste0("y_",  1:visits)])))
y_min <- floor(min(unlist(sim_data[, paste0("y_",  1:visits)])))

# Create Probability-Scale Outcome
sim_data[, paste0("y_pr_", 1:visits)] <-
  sapply(
    X = sim_data[, paste0("y_", 1:visits)],
    FUN = function(x) compress_range(x = x, min = y_min, max = y_max)
  )

# Create Binary Outcome
sim_data[, paste0("y_bin_", 1:visits)] <-
  sapply(
    X = sim_data[, paste0("y_pr_", 1:visits)],
    FUN = function(x) rbinom(n = length(x), size = 1, prob = x)
  )

# Introduce Missingness
sim_data[, paste0("y_pr_obs_", 1:visits)] <-
  sim_data[, paste0("y_pr_", 1:visits)]*
  sapply(
    X = sim_data[, paste0("y_obs_",  1:visits)],
    FUN = function(x) ifelse(test = is.na(x), yes = NA, no = 1)
  )

sim_data[, paste0("y_bin_obs_", 1:visits)] <-
  sim_data[, paste0("y_bin_", 1:visits)]*
  sapply(
    X = sim_data[, paste0("y_obs_",  1:visits)],
    FUN = function(x) ifelse(test = is.na(x), yes = NA, no = 1)
  )




### Set up Propensity & IPW Formulas ###########################################

propensity_score_formula <-
  tx ~ x_1 + x_2 + x_3 + x_4

inverse_weight_formulas <-
  list(
    ~ x_1 + x_2 + x_3 + x_4 + tx,
    ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1,
    ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1 + y_obs_2
  )

outcome_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1,
    y_obs_3 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1 + y_obs_2
  )




### 1.1 Continuous - No Imputation #############################################
test_1_1_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_monotone = FALSE,
    verbose = TRUE
  )


test_1_1_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_monotone = FALSE,
    verbose = TRUE
  )



### 1.2 Continuous - Gaussian Imputation #######################################
impute_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
  )

test_1_2_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_monotone = TRUE,
    impute_model = "gaussian",
    impute_formulas = 
      list(
        y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
        y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
      ),
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )


test_1_2_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_monotone = TRUE,
    impute_model = "gaussian",
    impute_formulas = 
      list(
        y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
        y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
      ),
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )




### 1.3 Continuous - PMM Imputation ############################################
test_1_3_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_monotone = TRUE,
    impute_model = "pmm",
    impute_formulas = 
      list(
        y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
        y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
      ),
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )


test_1_3_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_monotone = TRUE,
    impute_model = "pmm",
    impute_formulas = 
      list(
        y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
        y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
      ),
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )


### 1.1 - 1.3 Results ##########################################################
test_1_1_diff$ate
test_1_2_diff$ate
test_1_3_diff$ate

test_1_1_ratio$ate
test_1_2_ratio$ate
test_1_3_ratio$ate
