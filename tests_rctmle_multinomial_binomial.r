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

# Cut into Ordinal Outcome
outcome_breaks <- c(-Inf, 9, 10, 12, Inf)

sim_data_no_absorb <- sim_data

rm(sim_data)

sim_data_no_absorb$y_1 <- 
  cut(x = sim_data_no_absorb$y_1, breaks = outcome_breaks, labels = 1:4)
sim_data_no_absorb$y_2 <- 
  cut(x = sim_data_no_absorb$y_2, breaks = outcome_breaks, labels = 1:4)
sim_data_no_absorb$y_3 <- 
  cut(x = sim_data_no_absorb$y_3, breaks = outcome_breaks, labels = 1:4)

sim_data_no_absorb$y_obs_1 <- 
  cut(x = sim_data_no_absorb$y_obs_1, breaks = outcome_breaks, labels = 1:4)
sim_data_no_absorb$y_obs_2 <- 
  cut(x = sim_data_no_absorb$y_obs_2, breaks = outcome_breaks, labels = 1:4)
sim_data_no_absorb$y_obs_3 <- 
  cut(x = sim_data_no_absorb$y_obs_3, breaks = outcome_breaks, labels = 1:4)


# Introduce non-monotone missingness
non_monotone_rows <- which(sim_data_no_absorb$last_observed > 1)
sim_data_no_absorb$y_obs_1[sample(x = non_monotone_rows, size = n_non_monotone)] <- NA
sim_data_no_absorb$y_obs_2[sample(x = non_monotone_rows, size = n_non_monotone)] <- NA

# Make outcome level 4 absorbing
sim_data_absorb <- sim_data_no_absorb

absorbed <- with(sim_data_absorb, which(y_1 == 4))
sim_data_absorb[absorbed, c("y_2", "y_3")] <- 4

absorbed <- with(sim_data_absorb, which(y_2 == 4))
sim_data_absorb[absorbed, c("y_3")] <- 4

absorbed <- with(sim_data_absorb, which(y_obs_1 == 4))
sim_data_absorb[absorbed, c("y_obs_2", "y_obs_3")] <- 4

absorbed <- with(sim_data_absorb, which(y_obs_2 == 4))
sim_data_absorb[absorbed, c("y_obs_3")] <- 4

# Make final outcome binary
sim_data_no_absorb$y_3 <-
  with(sim_data_no_absorb, 1*(y_3 == 1 | y_3 == 2))

sim_data_no_absorb$y_obs_3 <-
  with(sim_data_no_absorb, 1*(y_obs_3 == 1 | y_obs_3 == 2))

sim_data_absorb$y_3 <-
  with(sim_data_absorb, 1*(y_3 == 1 | y_3 == 2))

sim_data_absorb$y_obs_3 <-
  with(sim_data_absorb, 1*(y_obs_3 == 1 | y_obs_3 == 2))



### Set up Propensity & IPW Formulas ###########################################

propensity_score_formula <-
  tx ~ s(x_1) + x_2 + x_3 + x_4

inverse_weight_formulas <-
  list(
    ~ s(x_1) + x_2 + x_3 + x_4 + tx,
    ~ s(x_1) + x_2 + x_3 + x_4 + tx + y_obs_1,
    ~ s(x_1) + x_2 + x_3 + x_4 + tx + y_obs_1 + y_obs_2
  )




### 1. Multinomial-Binomial - No absorbing state ###############################

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
    data = sim_data_no_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial_binomial",
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_1_diff$ate


test_1_ratio <-
  rctmle(
    data = sim_data_no_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial_binomial",
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_1_ratio$ate


test_1_ratio <-
  rctmle(
    data = sim_data_no_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial_binomial",
    estimand = "oddsratio",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_1_ratio$ate




### 2. Multinomial-Binomial - With absorbing state ###############################

test_2_diff <-
  rctmle(
    data = sim_data_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial-binomial",
    absorbing_state = 4,
    absorbing_outcome = 0,
    estimand = "difference",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_2_diff$ate


test_2_ratio <-
  rctmle(
    data = sim_data_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial-binomial",
    absorbing_state = 4,
    absorbing_outcome = 0,
    estimand = "ratio",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_2_ratio$ate


test_2_ratio <-
  rctmle(
    data = sim_data_absorb,
    propensity_score_formula = propensity_score_formula,
    inverse_weight_formulas = inverse_weight_formulas,
    outcome_formulas = outcome_formulas,
    outcome_type = "multinomial_binomial",
    estimand = "oddsratio",
    impute_formulas = impute_formulas,
    impute_model = "multinomial",
    verbose = TRUE
  )

test_2_ratio$ate




### 3. Error Handling Checks ####################################################
