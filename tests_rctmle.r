rm(list = ls())

set.seed(12345)

source(
  file.path(
    "https://raw.githubusercontent.com/jbetz-jhu/rctmle/main/rctmle.r"
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




### 1. Continous - Impute Parametric - Gaussian/Identity #######################

outcome_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1,
    y_obs_3 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1 + y_obs_2
  )

impute_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_obs_1
  )

test_1_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )

test_1_diff$ate


test_1_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = gaussian(link = "identity")
      ),
    verbose = TRUE
  )

test_1_ratio$ate


### 2. Continous - Impute Parametric - Gaussian/Log ###########################
test_2_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "gaussian",
    imputation_args = 
      list(
        family = gaussian(link = "log")
      ),
    verbose = TRUE
  )

test_2_diff$ate


test_2_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "gaussian",
    imputation_args = 
      list(
        family = gaussian(link = "log")
      ),
    verbose = TRUE
  )

test_2_ratio$ate


### 3. Continous - Impute PMM - Gaussian/Log ###################################
test_3_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasipoisson(link = "log")
      ),
    verbose = TRUE
  )

test_3_diff$ate


test_3_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "gaussian",
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasipoisson(link = "log")
      ),
    verbose = TRUE
  )

test_3_ratio$ate


### 4. Logistic - Impute PMM - Quasibinomial/Logit #############################
test_4_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasibinomial(link = "logit")
      ),
    verbose = TRUE
  )

test_4_diff$ate


test_4_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasibinomial(link = "logit")
      ),
    verbose = TRUE
  )

test_4_ratio$ate


### 5. Logistic - Impute PMM - Quasi: Logit, Var = mu(1-mu) ####################
test_5_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasi(link = "logit", variance = "mu(1-mu)")
      ),
    verbose = TRUE
  )

test_5_diff$ate



test_5_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = quasi(link = "logit", variance = "mu(1-mu)")
      ),
    verbose = TRUE
  )

test_5_ratio$ate



### 6. Logistic - Impute PMM - Beta/Logit ######################################
test_6_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = betar,
        link = "logit"
      ),
    verbose = TRUE
  )

test_6_diff$ate


test_6_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "pmm",
    imputation_args = 
      list(
        family = betar,
        link = "logit"
      ),
    verbose = TRUE
  )

test_6_ratio$ate


### 7. Logistic - Impute Beta/Logit ############################################
test_7_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "beta",
    imputation_args = 
      list(
        family = betar,
        link = "logit"
      ),
    verbose = TRUE
  )

test_7_diff$ate


test_7_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "logistic",
    outcome_range = c(y_min, y_max),
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "beta",
    imputation_args = 
      list(
        family = betar,
        link = "logit"
      ),
    verbose = TRUE
  )

test_7_diff$ate


### 8. Binomial - Impute Binomial/Logit ########################################
outcome_formulas <-
  list(
    y_bin_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_bin_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_bin_obs_1,
    y_bin_obs_3 ~ x_1 + x_2 + x_3 + x_4 + tx + y_bin_obs_1 + y_bin_obs_2
  )

impute_formulas <-
  list(
    y_bin_obs_1 ~ x_1 + x_2 + x_3 + x_4 + tx,
    y_bin_obs_2 ~ x_1 + x_2 + x_3 + x_4 + tx + y_bin_obs_1
  )

inverse_weight_formulas <-
  list(
    ~ x_1 + x_2 + x_3 + x_4 + tx,
    ~ x_1 + x_2 + x_3 + x_4 + tx + y_bin_obs_1,
    ~ x_1 + x_2 + x_3 + x_4 + tx + y_bin_obs_1 + y_bin_obs_2
  )


test_8_diff <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "binomial",
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "binomial",
    imputation_args = 
      list(
        family = binomial,
        link = "logit"
      ),
    verbose = TRUE
  )

test_8_diff$ate


test_8_ratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "binomial",
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "binomial",
    imputation_args = 
      list(
        family = binomial,
        link = "logit"
      ),
    verbose = TRUE
  )

test_8_ratio$ate


test_8_oddsratio <-
  rctmle(
    data = sim_data,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "binomial",
    estimand = "oddsratio",
    impute_formulas = impute_formulas,
    impute_model = "binomial",
    imputation_args = 
      list(
        family = binomial,
        link = "logit"
      ),
    verbose = TRUE
  )

test_8_oddsratio$ate
