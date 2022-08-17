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


### Continuous Outcome - Null Treatment ########################################
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




### TESTS - GAUSSIAN ###########################################################
# Use outcome on original scale
impute_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx + y_obs_1
  )

# Use y_obs_1, ..., y_obs_T
impute_columns <-
  non_monotone(
    data = sim_data,
    y_columns = paste0("y_obs_", 1:visits)
  )

# Gaussian - Identity Link
impute_gam(
    data = sim_data,
    impute_columns = impute_columns,
    impute_formulas = impute_formulas,
    verbose = TRUE,
    model = "gaussian",
    family = gaussian(link = "identity")
  )

# Gaussian - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "gaussian",
  family = gaussian(link = "log")
)

# Gaussian - Inverse Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "gaussian",
  family = gaussian(link = "inverse")
)


# Test Smooth Terms & Interactions
impute_formulas <-
  list(
    y_obs_1 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx,
    y_obs_2 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + s(y_obs_1, by = tx)
  )

# Gaussian - Identity Link: Smooth Terms & Interactions
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "gaussian",
  family = gaussian(link = "identity")
)




## TESTS - PMM ################################################################
# Use outcome on original scale
impute_formulas <-
  list(
    y_obs_1 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx,
    y_obs_2 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx + y_obs_1
  )

# Gaussian - Identity Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "identity")
)

# Gaussian - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "log")
)

# Gaussian - Inverse Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "inverse")
)

# Quasipoisson - Identity Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = quasipoisson(link = "identity")
)

# Quasipoisson - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = quasipoisson(link = "log")
)

# Gamma - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = Gamma(link = "log")
)


# Test Smooth Terms & Interactions
impute_formulas <-
  list(
    y_obs_1 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx,
    y_obs_2 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + s(y_obs_1, by = tx)
  )

# PMM - Quasipoisson Log Link: Smooth Terms & Interactions
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = quasipoisson(link = "log")
)


# Use outcome on probability scale
impute_formulas <-
  list(
    y_pr_obs_1 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx,
    y_pr_obs_2 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx + y_obs_1
  )


# Use y_pr_obs_1, ..., y_pr_obs_T
impute_columns <-
  non_monotone(
    data = sim_data,
    y_columns = paste0("y_pr_obs_", 1:visits)
  )

# Gaussian - Identity Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "identity")
)

# Gaussian - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "log")
)

# Gaussian - Inverse Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = gaussian(link = "inverse")
)

# Quasi - Logit Link: Variance = mu(1-mu)
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = quasi(link = "logit", variance =  "mu(1-mu)")
)

# Beta - Logit Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = betar(link = "logit")
)


# Test Smooth Terms & Interactions
impute_formulas <-
  list(
    y_pr_obs_1 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx,
    y_pr_obs_2 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + s(y_pr_obs_1, by = tx)
  )

# PMM - Quasipoisson Log Link: Smooth Terms & Interactions
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "pmm",
  family = quasibinomial(link = "logit")
)




## TESTS - Beta ################################################################
# Use outcome on probability scale
impute_formulas <-
  list(
    y_pr_obs_1 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx,
    y_pr_obs_2 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx + y_pr_obs_1
  )

# Use y_pr_obs_1, ..., y_pr_obs_T
impute_columns <-
  non_monotone(
    data = sim_data,
    y_columns = paste0("y_pr_obs_", 1:visits)
  )

# Beta - Logit Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "beta",
  family = betar(link = "logit")
)

# Beta - Probit Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "beta",
  family = betar(link = "probit")
)

# Beta - Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "beta",
  family = betar(link = "log")
)

# Beta - Complementary Log-Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "beta",
  family = betar(link = "cloglog")
)


# Test Smooth Terms & Interactions
impute_formulas <-
  list(
    y_pr_obs_1 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx,
    y_pr_obs_2 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx + s(y_pr_obs_1)
  )

# Beta - Logit Link: Smooth Terms & Interactions
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "beta",
  family = betar(link = "cloglog")
)



## TESTS - Binomial ############################################################
# Use outcome on probability scale
impute_formulas <-
  list(
    y_bin_obs_1 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx,
    y_bin_obs_2 ~ x_1 + x_2 + x_3 + x_4 + x_5 + tx + y_bin_obs_1
  )

# Use y_pr_obs_1, ..., y_pr_obs_T
impute_columns <-
  non_monotone(
    data = sim_data,
    y_columns = paste0("y_bin_obs_", 1:visits)
  )

# Beta - Logit Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "binomial",
  family = binomial(link = "logit")
)

# Beta - Probit Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "binomial",
  family = binomial(link = "probit")
)

### Log Link does not converge: known problem with log link
# # Beta - Log Link
# impute_gam(
#   data = sim_data,
#   impute_columns = impute_columns,
#   impute_formulas = impute_formulas,
#   verbose = TRUE,
#   model = "binomial",
#   family = binomial(link = "log")
# )

# Beta - Complementary Log-Log Link
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "binomial",
  family = binomial,
  link = "cloglog"
)



# Test Smooth Terms & Interactions
impute_formulas <-
  list(
    y_bin_obs_1 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx,
    y_bin_obs_2 ~ s(x_1) + x_2 + x_3 + x_4 + x_5 + tx + y_bin_obs_1
  )

# Beta - Logit Link: Smooth Terms & Interactions
impute_gam(
  data = sim_data,
  impute_columns = impute_columns,
  impute_formulas = impute_formulas,
  verbose = TRUE,
  model = "binomial",
  family = binomial,
  link = "logit"
)