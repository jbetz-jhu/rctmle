library(mvtnorm) # for rmvnorm in iid_centered_mvn()
library(mgcv) # for impute_gam
library(nnet) # for impute_monotone_multinomial
library(boot) # for tmle - bootstrap inference

### compress_range #############################################################
compress_range <- 
  function(
    x,
    min,
    max,
    enforce_range = TRUE
  ) {
    
    if(sum(x > max, na.rm = TRUE) + sum(x < min, na.rm = TRUE) > 0) {
      warning("x contains values outside of specified range")
      if(enforce_range){ 
        warning("values outside of specified range replaced with boundaries")
        x[which(x > max)] <- max
        x[which(x < min)] <- min
      }
    }
    
    (x - min)/(max - min)
  }




### decompress_range ###########################################################
decompress_range <-
  function(
    x,
    min,
    max
  ) {
    x*(max - min) + min
  }




### iid_centered_mvn ###########################################################
iid_centered_mvn <-
  function(
    n,
    n_covariates,
    sigma_x
  ) {
    setNames(
      object =
        data.frame(
          mvtnorm::rmvnorm(
            n = n,
            sigma = sigma_x
          )
        ),
      nm = paste0("x_", 1:n_covariates)
    )
  }




### sim_lme_trial ##############################################################
# Create a data.frame of n simulated participants with p baseline covariates and
#   T continuous outcomes, based on a linear mixed effects (LME) model.
#   Simulated data can incorporate missing data due to study dropout.
#
# n: number of participants
# T: times of outcomes: (T-1) intermediate and one primary outcome
# p: number of baseline covariates
#
# NOTE: First the complete data vector is constructed, and then missingness is
# applied in the last step - missing and complete datasets are returned. If MAR
# mechanism is used, missingness should only depend on previous outcomes and
# baseline covariates.
#
# How to specify outcome and missingness regression models: Suppose the complete
# data vector for each participant consists of 5 baseline covariates and 3
# outcomes. For simplicity, the covariates are named X1, X2, ..., X5.
#
# If outcome_cov is not specified, the covariates will be unrelated to the
# outcome. Elements should be a named vector, with names corresponding to the
# names of baseline covariates or treatment ("tx"). These are used to create
# linear predictors, and the outcome is sampled from a normal distribution given
# the linear predictors and random effects.
#
# Models for dropout are specified similarly. Probabilities are supplied for
# the probability of dropping out at a given visit when all covariates are zero.
# Logistic regression coefficients can be supplied for baseline covariates,
# treatment assignment ("tx"), or outcomes (e.g. "y_1", "y_2", "y_3", etc.),
# however, for the MAR mechanism, missingness should only depend on previously
# measured variables.
# 
# input:
# bl_covariates: (n x p data.frame) of continuous or binary indicators
#   representing covariates measured prior to baseline with no missing data.
# mean_outcomes: (numeric of length T)) - the mean outcome at each assessment
#   for an individual with all covariates and random effects equal to 0.
# outcome_cov: (list of length (T-1)) - regression coefficients
#   for associations between prior data and outcomes, with names corresponding
#   to the columns for baseline covariates, 'treatment', or prior outcomes
#   ('y_1', 'y_2', ...). If NULL (default), outcomes will have no association
#   with either baseline covariates, treatment assignment, or intermediate
#   outcomes.  If a single element in list is NULL, no covariates will have
#   an association with outcomes at that visit
# re_variance: (numeric of length 1-2) - the variance of the random slopes and
#   intercepts.
# re_correlation: correlation between random slopes and random intercepts: only
#   necessary when both random slopes and intercepts are supplied (re_variance
#   has length of 2).
# residual_sd: residual standard deviation
# pr_dropout: (1 x T probability) rate of dropout PRIOR TO each visit for those
#   with the default level of all covariates. If NULL (default), no dropout is
#   introduced.
# dropout_cov: (list of length (T-1)) - logistic regression coefficients
#   for associations between dropout and data (covariates, treatment
#   assignment, or outcomes) at each visit, with names corresponding to the
#   appropriate variables. If NULL (default), dropout will have no association
#   with either baseline covariates, treatment assignment, or intermediate
#   outcomes. If a single element in list is NULL, no covariates will have
#   an association with dropout at that visit.
#
# NOTE: For missingness models, outcomes are CENTERED, but not scaled: If an
#   outcome is included in the model, the intercept is the probability of
#   missingness for individuals with covariates equal to 0 and the population
#   mean value for the outcome.

sim_lme_trial <-
  function(
    bl_covariates,
    visit_times,
    mean_outcomes,
    outcome_cov = NULL,
    re_variance,
    re_correlation = NULL,
    residual_sd,
    pr_dropout = NULL,
    dropout_cov = NULL
  ) {
    
    n_random_effects <-
      nrow(as.matrix(re_variance))
    
    stopifnot(length(visit_times) == length(mean_outcomes))
    
    if(n_random_effects == 1) {
      # Only random intercept is supplied
      re_covariance <- re_variance
    } else if(n_random_effects == 2) {
      # Slope/intercept supplied
      if(is.null(re_correlation)) { # Uncorrelated random effects
        re_correlation <- diag(n_random_effects)
      } else if (nrow(as.matrix(re_correlation)) == 1) { # Scalar correlation
        re_correlation <-
          diag(2) +
          re_correlation*(matrix(data = 1, nrow = 2, ncol = 2) - diag(2))
      } 
      
      re_covariance <- 
        (matrix(sqrt(re_variance), ncol = 1) %*% 
           matrix(sqrt(re_variance), nrow = 1))*re_correlation
      
    } else {
      stop(paste0("Only random intercept, random slope, or random intercept ",
                  "and slope are implemented."))
    }
    
    
    n_outcomes <- length(x = visit_times)
    n_obs <- nrow(bl_covariates)
    
    # 1. Treatment is assigned
    study_data <-
      data.frame(
        id = 1:nrow(bl_covariates),
        bl_covariates,
        tx =
          sample(x = c(rep(0, ceiling(n_obs/2)),
                       rep(1, n_obs - ceiling(n_obs/2))),
                 size = n_obs),
        last_observed = NA,
        dropout_i = NA
      )
    
    # 2. Determine Distribution of outcomes
    
    # 2.1 Determine Fixed Effects
    out_cov_lp <- matrix(0, nrow = n_obs, ncol = n_outcomes)
    
    if(!is.null(outcome_cov)){
      for(i in 1:n_outcomes) {
        if(!is.null(outcome_cov[[i]])) {
          out_cov_lp[, i] <-
            as.matrix(study_data[, names(outcome_cov[[i]])]) %*% 
            outcome_cov[[i]]
        }
      }
    }
    
    # 2.2 Simulate Random effects are simulated
    # Sample random effects
    random_effects <-
      rmvnorm(
        n = n_obs,
        sigma = as.matrix(re_covariance)
      )
    
    # Create (n x T) matrix of disturbances due to random effects
    if(ncol(random_effects) == 1){ # Random Intercept
      random_trajectories <-
        kronecker(random_effects,
                  matrix(1, ncol = length(visit_times)))
      colnames(random_effects) <- c("random_intercept")
    } else { # Centered Random Effects at Mean Visit Time
      random_trajectories <-
        random_effects %*%
        rbind(1, visit_times - mean(visit_times))
      colnames(random_effects) <- c("random_intercept", "random_slope")
    }
    
    # 2.3. Add fixed effects and random effects to residual
    residuals <-
      rmvnorm(n = n_obs,
              mean = matrix(0, ncol = n_outcomes),
              sigma = diag(residual_sd^2, nrow = n_outcomes)) # Residual
    colnames(residuals) <- paste0("residual_", 1:n_outcomes)
    
    study_data[, paste0("y_", 1:n_outcomes)] <-
      study_data[, paste0("y_obs_", 1:n_outcomes)] <-
      # Note: In order to make the dropout model easier to use, the mean
      # outcome is added in later.
      matrix(data = 0, nrow = n_obs, ncol = n_outcomes) +
      out_cov_lp + # Linear Predictor
      random_trajectories + # Random Effects
      residuals
    
    study_data <-
      data.frame(study_data, random_effects, residuals)
    
    # 3. dropout is assigned
    for(i in 1:n_outcomes) {
      
      if(is.null(pr_dropout[[i]])){
      } else if(pr_dropout[[i]] > 0 & pr_dropout[[i]] < 1) {
        
        if(is.null(dropout_cov[[i]])) {
          drop_cov_lp <- qlogis(1 - pr_dropout[[i]])
        } else {
          drop_cov_lp <- qlogis(1 - pr_dropout[[i]]) -
            as.matrix(study_data[names(dropout_cov[[i]])]) %*% 
            dropout_cov[[i]]
        }
        
        study_data$dropout_i <- 
          rlogis(n = n_obs, location = drop_cov_lp, scale = 1) < 0
        
        dropout.rows <-
          with(study_data, which(is.na(last_observed) & dropout_i))
        
        study_data[dropout.rows,
                   paste0("y_obs_", i:n_outcomes)] <- NA
        
        study_data$last_observed[dropout.rows] <- (i - 1)
      }
    }
    
    study_data$last_observed[which(is.na(study_data$last_observed))] <- 
      n_outcomes
    
    study_data$dropout_i <- NULL
    
    # 4. Add in mean outcome
    study_data[, paste0("y_", 1:n_outcomes)] <-
      study_data[, paste0("y_", 1:n_outcomes)] + 
      matrix(data = mean_outcomes, nrow = n_obs, ncol = n_outcomes,
             byrow = TRUE)
    
    study_data[, paste0("y_obs_", 1:n_outcomes)] <-
      study_data[, paste0("y_obs_", 1:n_outcomes)] + 
      matrix(data = mean_outcomes, nrow = n_obs, ncol = n_outcomes,
             byrow = TRUE)
    
    study_data
  }




### impute_covariates_mean_mode ################################################
# Uses mean for numeric covariates, mode for non-numeric covariates.
impute_covariates_mean_mode <-
  function(data, x_columns) {
    
    Mode <- function(x) {
      ux <- unique(x)
      ux[which.max(tabulate(match(x, ux)))]
    }
    
    any_missing <- function(x){
      any(is.na(x))
    }
    
    impute_mean <- function(x) {
      ifelse(
        test = is.na(x),
        yes = mean(x, na.rm = TRUE),
        no = as.numeric(x)
      )
    }
    
    impute_mode <- function(x) {
      ifelse(
        test = is.na(x),
        yes = Mode(x),
        no = x
      )
    }
    
    numeric_vars <-
      names(which(sapply(X = data[x_columns], FUN = is.numeric)))
    
    non_numeric_vars <-
      setdiff(x_columns, numeric_vars)
    
    for(i in numeric_vars){
      data[, i] <- impute_mean(data[, i])
    }
    
    for(i in non_numeric_vars){
      data[, i] <- impute_mode(data[, i])
    }
    
    return(data)
  }




### non_monotone ###############################################################
# determine which missing values are non-monotone
non_monotone <- function(data, y_columns){
  n_outcomes <- length(y_columns)
  
  non_monotone <-
    matrix(
      data = NA,
      nrow = nrow(data),
      ncol = length(y_columns) - 1
    )
  
  for(i in 1:ncol(non_monotone)){
    
    non_monotone[, i] <-
      is.na(data[, y_columns[i]]) &
      if(i == ncol(non_monotone)) {
        !is.na(data[, y_columns[(i + 1):n_outcomes]])
      } else{
        rowSums(
          !is.na(data[, y_columns[(i + 1):n_outcomes]])
        ) > 0
      }
  }
  
  return(
    setNames(
      object = data.frame(non_monotone),
      nm = head(y_columns, - 1)
    )
  )
}




### absorbing_state_check ######################################################
absorbing_state_check <-
  function(
    data,
    y_columns,
    outcome_type,
    absorbing_state,
    absorbing_outcome
  ) {
    
    if(outcome_type != "multinomial-binomial"){
      stop("Absorbing state only for multinomial-binomial outcomes.")
    } else if(is.null(absorbing_state) | is.null(absorbing_outcome)){
      stop("'absorbing_state' and 'absorbing_outcome' must be supplied.")
    }
    
    final_outcome <- tail(x = y_columns, 1)
    intermediate_outcomes <- setdiff(x = y_columns, y = final_outcome)
    
    if(length(intermediate_outcomes > 0)){
      
      first_absorbing <-
        apply(
          X = data[, head(x = y_columns, -1)],
          FUN = function(x) which(x %in% absorbing_state)[1],
          MARGIN = 1
        )
      
      last_non_absorbing <-
        apply(
          X = data[, head(x = y_columns, -1)],
          FUN = function(x) which(!(x %in% c(NA, absorbing_state)))[1],
          MARGIN = 1
        )
      
      absorbing_conflicts <-
        which(last_non_absorbing > first_absorbing)
      
      if(length(absorbing_conflicts) > 0){
        stop("Transition out of absorbing state in rows: ",
             paste0(absorbing_conflicts, collapse = ", "))
      }
      
      # Carry forward absorbing state
      for(i in 1:(length(y_columns) - 1)){
        rows <- which(first_absorbing == i)
        data[rows, intermediate_outcomes] <- absorbing_state
        data[rows, final_outcome] <- absorbing_outcome
      }
    }
    
    return(data)
  }




### impute_outcomes_to_monotone ################################################
impute_outcomes_to_monotone <-
  function(
    ...
  ) {
    
    arg_list <- as.list(substitute(list(...)))[-1L]
    
    if(arg_list$model %in% c("gaussian", "beta", "binomial", "pmm")) {
      do.call(
        what = impute_gam,
        args = arg_list[formalArgs(impute_gam)] 
      )
    } else if(arg_list$model == "multinomial") {
      if(is.null(arg_list$absorbing_state)){
        do.call(
          what = impute_multinomial,
          args = arg_list[formalArgs(impute_multinomial)]
        )
      } else {
        do.call(
          what = impute_multinomial_absorbing,
          args = arg_list[formalArgs(impute_multinomial_absorbing)] 
        )
      }
    } else{
      stop("`model` not recognized.")
    }
  }




### impute_gam #################################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# GAMs using the `mgcv` package. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes - This is checked
# by tmle_precheck().
#
# data: data.frame to impute
# impute_columns: a data.frame indicating which rows to impute for each outcome
# impute_formulas: a list of formulas, specifying imputation models.
# verbose: logical - should imputation models and parameters be returned?
# model: character - indicating how imputations should be generated, either from
#   a probability model (gaussian/beta/binomial) or using predictive mean
#   matching (pmm) to resample from observed values.
# impute_family/impute_link: family and link - see ?mgcv::gam()
# donors: numeric - number of candidates to resample from at random in PMM:
#   ignored for parametric imputation.
impute_gam <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    stochastic = TRUE,
    propagate = TRUE,
    model,
    family,
    donors = 10,
    verbose = FALSE,
    outcome_range = NULL,
    ...
  ) {
    if(!propagate){
      original_data <- data
    }
    
    # Get names of variables to be imputed from LHS of formulas
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- imputation_parameters <- list()
    
    for(i in 1:length(impute_formulas)){
      
      # Fit Imputation Model among observed cases
      impute_model <-
        mgcv::gam(
          formula = impute_formulas[[i]],
          data =
            if(propagate){
              data
            } else {
              original_data
            },
          family = family,
          ...
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      pred_interval <-
        predict(
          object = impute_model,
          newdata = data[which(impute_columns[, i]),],
          type = "response",
          se.fit = (model == "gaussian" & stochastic)
        )
      
      if(!stochastic){
        
        if(!is.null(outcome_range)){
          pred_interval <-
            pmin(pmax(pred_interval, outcome_range[1]), outcome_range[2])
        }
        
        data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
          pred_interval
        
      } else if(model == "gaussian"){
        
        if(stochastic){
          pred_interval <-
            with(
              pred_interval,
              data.frame(fit, se.fit)
            )
          
          # Sample from Prediction Interval
          pred_interval$df <- impute_model$df.residual
          pred_interval$sigma_sq <- impute_model$sig2
          pred_interval$sd <- with(pred_interval, sqrt(se.fit^2 + sigma_sq))
          pred_interval$lcl <-
            with(pred_interval, fit + qt(p = 0.025, df = df)*sd)
          pred_interval$ucl <-
            with(pred_interval, fit + qt(p = 0.975, df = df)*sd)
          
          # Fill in imputed values
          imputations <-
            with(
              pred_interval,
              rnorm(
                n = n_to_impute[i],
                mean = fit,
                sd = sd
              )
            )
          
          imputations <-
            pmin(pmax(imputations, outcome_range[1]), outcome_range[2])
          
          data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
            imputations
        }
      } else {
        pred_interval <-
          data.frame(fit = pred_interval)
        
        
        if(model == "beta") {
          pred_interval$m <-
            impute_model$family$getTheta(trans = TRUE)
          pred_interval$alpha <- 
            with(pred_interval, fit*m)
          pred_interval$beta <- 
            with(pred_interval, (1 - fit)*m)
          
          data[which(impute_columns[,i]), impute_formula_lhs[i]] <-
            rbeta(
              n = n_to_impute[i],
              shape1 = pred_interval$alpha,
              shape2 = pred_interval$beta,
            )
          
        } else if(model == "binomial"){
          data[which(impute_columns[,i]), impute_formula_lhs[i]] <-
            rbinom(
              n = n_to_impute[i],
              size = 1,
              prob = pred_interval$fit
            )
        } else if(model == "pmm"){
          # Get all non-missing values for candidates
          candidates <-
            na.omit(data[, impute_formula_lhs[i]])
          
          # Find distance between fit and candidates 
          candidate_distance <-
            kronecker(
              X = candidates,
              Y = matrix(data = 1, ncol = length(pred_interval$fit))
            ) - 
            kronecker(
              X = matrix(data = 1, nrow = length(candidates)),
              Y = matrix(data = pred_interval$fit, nrow = 1)
            )
          
          selected_candidates <-
            apply(
              X = candidate_distance^2,
              MARGIN = 2,
              FUN = function(x, donors)
                sample(x = which(rank(x) < donors), size = 1),
              donors = donors
            )
          
          data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
            candidates[selected_candidates]
        }
      }
      
      if(verbose) imputation_parameters[[i]] <- pred_interval
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models,
          imputation_parameters = imputation_parameters
        )
      )
    } else{
      return(data)
    }
  }




### impute_multinomial #########################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# multinomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
# Use `impute_multinomial_absorbing` if the outcome scale has an
# absorbing state (e.g. Death) that should be preserved
impute_multinomial <-
  function(
    data = data,
    impute_columns = impute_columns,
    impute_formulas,
    verbose = FALSE
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_to_impute <- colSums(impute_columns)
    
    imputation_models <- list()
    
    for(i in 1:length(impute_formulas)){
      
      outcome_levels <- levels(data[, impute_formula_lhs[i]])
      
      # Fit Imputation Model among observed cases
      impute_model <-
        nnet::multinom(
          formula = impute_formulas[[i]],
          data = data,
          trace = FALSE
        )
      
      if(verbose) imputation_models[[i]] <- impute_model
      
      # Fill in imputed values
      imputations <-
        apply(
          X = 
            predict(
              object = impute_model,
              newdata = data[which(impute_columns[, i]), ],
              type = "probs"
            ),
          MARGIN = 1,
          FUN = function(x, outcome_levels) 
            outcome_levels[
              which(rmultinom(n = 1, size = 1, prob = x) == 1)
            ],
          outcome_levels = outcome_levels
        )
      
      data[which(impute_columns[, i]), impute_formula_lhs[i]] <-
        imputations
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }




### impute_multinomial_absorbing ################################################
# INTERNAL FUNCTION: Impute non-monotone missing data to monotone using
# multinomial GLM with user-specified link. `impute_columns` and the LHS of
# `impute_formulas` must match temporal ordering of outcomes -
# This is checked by tmle_precheck().
# Use `impute_multinomial_absorbing` if the outcome scale has an
# absorbing state (e.g. Death) that should be preserved
impute_multinomial_absorbing <-
  function(
    data = data,
    impute_columns,
    y_columns,
    impute_formulas,
    absorbing_state,
    verbose = FALSE
  ) {
    
    impute_formula_lhs <-
      sapply(
        X = impute_formulas,
        FUN = function(x) all.vars(update(x, . ~ 0))
      )
    
    n_outcomes <- ncol(impute_columns)
    n_to_impute <- colSums(impute_columns)
    all_outcome_levels <-
      unique(
        unlist(
          lapply(
            X = data[, names(impute_columns)],
            FUN = levels
          )
        )
      )
    
    imputation_models <- list()
 
    for(i in 1:length(impute_formulas)){
      outcome_levels <- levels(data[, impute_formula_lhs[i]])
      non_absorbing <- setdiff(x = outcome_levels, y = absorbing_state)
      
      outcome_i <- which(y_columns == impute_formula_lhs[i])
      
      # Range of imputations depends on the next observed value:
      # If next observed outcome is absorbing state, use full scale.
      # If next observed outcome is NOT absorbing state, impute from
      # non-absorbing states to avoid inconsistency
      
      next_observed_outcome <-
        apply(
          X = as.matrix(data[, tail(x = y_columns, -outcome_i)]),
          MARGIN = 1,
          FUN = function(x, table = all_outcome_levels)
            min(match(x = x, table = table))
        )

      impute_full <-
        (impute_columns[, impute_formula_lhs[i]]) &
        (next_observed_outcome %in% absorbing_state)
      
      impute_non_absorbing <-
        (impute_columns[, impute_formula_lhs[i]]) &
        !(next_observed_outcome %in% absorbing_state)
      
      
      if(sum(impute_full) > 0){
        
        impute_model <-
          nnet::multinom(
            formula = impute_formulas[[i]],
            data = data,
            trace = FALSE
          )
        
        if(verbose)
          imputation_models[[i]] <-
            list(full_range = impute_model)
        
        imputations <-
          apply(
            X = predict(
              object = impute_model,
              newdata = data[which(impute_full), ],
              type = "probs"
            ),
            MARGIN = 1,
            FUN = function(x, outcome_levels) 
              outcome_levels[
                which(rmultinom(n = 1, size = 1, prob = x) == 1)
              ],
            outcome_levels = outcome_levels
          )
        
        
        # Carry absorbing state forward
        if(i < length(impute_formulas)) {
          new_absorbed_rows <-
            which(impute_full)[which(imputations %in% absorbing_state)]
          
          data[new_absorbed_rows, tail(x = y_columns, -outcome_i)] <-
            absorbing_state
          
          impute_columns[new_absorbed_rows, -c(1:outcome_i)] <- TRUE
        }
        
        data[which(impute_full), impute_formula_lhs[i]] <-
          imputations
      }
      
      if(sum(impute_non_absorbing) > 0){
        
        # Subset to those not absorbed in current or next outcome
        non_absorbed <- 
          data[
            which(!data[, y_columns[outcome_i]] %in% absorbing_state), 
          ]
        
        for(j in y_columns[1:outcome_i]){
          non_absorbed[, j] <- droplevels(non_absorbed[, j])
        }
        
        impute_model <-
          nnet::multinom(
            formula = impute_formulas[[i]],
            data = non_absorbed,
            trace = FALSE
          )
        
        if(verbose)
          imputation_models[[i]] <-
            list(non_absorbing = impute_model)
        
        imputations <-
          predict(
            object = impute_model,
            newdata = data[which(impute_non_absorbing), ],
            type = "probs"
          )
        
        
        
        
        if(!is.matrix(imputations)){
          # If only 2 classes, multinom produces vector, not matrix
          # Need probability for each class
          if(length(outcome_levels) == 2) {
            imputations <- cbind(1 - imputations, imputations)
          } else {
            imputations <- matrix(data = imputations, nrow = 1)
          }
        }
        
        imputations <-
          apply(
            X = imputations,
            MARGIN = 1,
            FUN = function(x, outcome_levels) 
              outcome_levels[
                which(rmultinom(n = 1, size = 1, prob = x) == 1)
              ],
            outcome_levels =
              levels(non_absorbed[, impute_formula_lhs[i]])
          )
        
        data[which(impute_non_absorbing), impute_formula_lhs[i]] <-
          imputations
      }
    }
    
    if(verbose){
      return(
        list(
          imputed_data = data,
          imputation_models = imputation_models
        )
      )
    } else{
      return(data)
    }
  }




### compute_inverse_weights ####################################################
compute_inverse_weights <-
  function(
    data,
    inverse_weight_formulas,
    inverse_weight_stratified = FALSE,
    inverse_weight_prediction_tx = "counterfactual",
    y_columns,
    tx_column,
    verbose = FALSE,
    absorbing_state = NULL
  ) {
    
    n_outcomes <- length(y_columns)
    
    ipw_models <- list()
    
    missing_outcomes <- 
      colSums(is.na(data[, y_columns]))
    
    # No missing outcomes
    if (all(missing_outcomes == 0)) {
      # No missing data
      ip_weights_tx_1 <- ip_weights_tx_0 <-
        matrix(
          data = 1,
          nrow = nrow(data),
          ncol = (n_outcomes + 1)
        )
      
    } else {
      
      ip_weights_tx_1 <- ip_weights_tx_0 <-
        matrix(
          data = NA,
          nrow = nrow(data),
          ncol = (n_outcomes + 1)
        )
      
      ip_weights_tx_1[, 1] <- ip_weights_tx_0[, 1] <- 1
      
      for(i in 1:n_outcomes){
        
        # If no dropouts or no new dropouts, carry weight forward
        # Find those previously absorbed or censored
        if(i > 1){
          prev_absorbed <-
            data[, y_columns[(i - 1)]] %in% absorbing_state
          
          prev_censored <-
            is.na(data[, y_columns[(i - 1)]])
          
          carry_weights_forward <- 
            missing_outcomes[i] == 0 |
            diff(missing_outcomes[(i - 1):i]) == 0
          
        } else {
          prev_absorbed <- prev_censored <-
            rep(FALSE, nrow(data))
          
          carry_weights_forward <- 
            missing_outcomes[i] == 0
        }
        
        uncensored <- which(!prev_censored)
        

        if(carry_weights_forward) {
          ip_weights_tx_1[uncensored, (i + 1)] <-
            ip_weights_tx_1[uncensored, i]
          
          ip_weights_tx_0[uncensored, (i + 1)] <-
            ip_weights_tx_0[uncensored, i]
        } else {
          
          if(is.null(absorbing_state)){
            uncensored_unabsorbed <- uncensored
          } else {
            uncensored_unabsorbed <-
              which(!(prev_absorbed|prev_censored))
            
            # Carry weights forward for those absorbed
            if(sum(prev_absorbed) > 0){
              ip_weights[which(prev_absorbed), (i + 1)] <-
                ip_weights[which(prev_absorbed), i]
            }
          }
          
          if(inverse_weight_stratified){
            
            uncensored_unabsorbed_tx_1 <-
              intersect(uncensored_unabsorbed, which(data[, tx_column] == 1))
            
            uncensored_unabsorbed_tx_0 <-
              intersect(uncensored_unabsorbed, which(data[, tx_column] == 0))
            
            ipw_model_tx_1 <-
              stats::glm(
                formula = inverse_weight_formulas[[i]],
                data = data[uncensored_unabsorbed_tx_1, ],
                family = binomial
              )
            
            ipw_model_tx_0 <-
              stats::glm(
                formula = inverse_weight_formulas[[i]],
                data = data[uncensored_unabsorbed_tx_0, ],
                family = binomial
              )
            
            if(inverse_weight_prediction_tx == "counterfactual") {
              pr_obs_fit_tx_1 <-
                predict(
                  object = ipw_model_tx_1,
                  newdata = 
                    within(
                      data = data[uncensored_unabsorbed,],
                      expr = {
                        eval(parse(text = paste0(tx_column, " = 1")))
                      }
                    ),
                  type = "response"
                )
              
              pr_obs_fit_tx_0 <-
                predict(
                  object = ipw_model_tx_0,
                  newdata = 
                    within(
                      data = data[uncensored_unabsorbed,],
                      expr = {
                        eval(parse(text = paste0(tx_column, " = 0")))
                      }
                    ),
                  type = "response"
                )
              
              ip_weights_tx_1[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_1[uncensored_unabsorbed, i]*
                (1/pr_obs_fit_tx_1)
              
              ip_weights_tx_0[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_0[uncensored_unabsorbed, i]*
                (1/pr_obs_fit_tx_0)
            
            } else if(inverse_weight_prediction_tx == "observed") {
              ip_weights_tx_1[uncensored_unabsorbed_tx_1, (i + 1)] <-
                ip_weights_tx_1[uncensored_unabsorbed_tx_1, i]*
                (1/ipw_model_tx_1$fitted.values)
              
              ip_weights_tx_1[uncensored_unabsorbed_tx_0, (i + 1)] <-
                ip_weights_tx_1[uncensored_unabsorbed_tx_0, i]*
                (1/ipw_model_tx_0$fitted.values)
              
              ip_weights_tx_0[, (i + 1)] <-
                ip_weights_tx_1[, (i + 1)]
            }

            if(verbose) {
              ipw_models[[i]] <- 
                list(
                  ipw_model_tx_1 = ipw_model_tx_1,
                  ipw_model_tx_0 = ipw_model_tx_0
                )
            }
            
          } else {
            ipw_model <-
              # mgcv::gam(
              stats::glm(
                formula = inverse_weight_formulas[[i]],
                data = data[uncensored_unabsorbed, ],
                family = binomial
              )
            
            if(verbose) ipw_models[[i]] <- ipw_model
            
            if(inverse_weight_prediction_tx == "counterfactual") {
              pr_obs_fit_tx_1 <-
                predict(
                  object = ipw_model,
                  newdata = 
                    within(
                      data = data[uncensored_unabsorbed,],
                      expr = {
                        eval(parse(text = paste0(tx_column, " = 1")))
                      }
                    ),
                  type = "response"
                )
              
              pr_obs_fit_tx_0 <-
                predict(
                  object = ipw_model,
                  newdata = 
                    within(
                      data = data[uncensored_unabsorbed,],
                      expr = {
                        eval(parse(text = paste0(tx_column, " = 0")))
                      }
                    ),
                  type = "response"
                )
              
              ip_weights_tx_1[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_1[uncensored_unabsorbed, i]*(1/pr_obs_fit_tx_1) 
              
              ip_weights_tx_0[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_0[uncensored_unabsorbed, i]*(1/pr_obs_fit_tx_0) 
              
            } else if(inverse_weight_prediction_tx == "observed") {
              pr_obs_fit <-
                predict(
                  object = ipw_model,
                  newdata = data[uncensored_unabsorbed,],
                  type = "response"
                )
              
              ip_weights_tx_1[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_1[uncensored_unabsorbed, i]*(1/pr_obs_fit) 
              
              ip_weights_tx_0[uncensored_unabsorbed, (i + 1)] <-
                ip_weights_tx_0[uncensored_unabsorbed, i]*(1/pr_obs_fit) 
            }
          }
        }
      }
    }
    
    if(verbose){
      return(
        list(
          # Drop leading column
          ip_weights_tx_1 = ip_weights_tx_1[, -1],
          ip_weights_tx_0 = ip_weights_tx_0[, -1],
          ipw_models = ipw_models
        )
      )
    } else{
      return(
        # Drop leading column
        list(
          ip_weights_tx_1 = ip_weights_tx_1[, -1],
          ip_weights_tx_0 = ip_weights_tx_0[, -1]
        )
      )
    }
  }




### tmle_get_args ##############################################################
# Take arguments from call to TMLE: determine which variables are baseline
# covariates (x_columns), which is the treatment indicator (tx_column),
# and which are outcomes (y_columns) based on the formulas supplied.
# If imputation formulas are supplied, auxilliary covariates (in imputation
# model but not in propensity score, inverse weight, or outcome models) are 
# aslo returned (imputation_x).
tmle_get_args <-
  function(
    data,
    propensity_score_formula,
    inverse_weight_formulas,
    outcome_formulas,
    impute_formulas = NULL
  ) {
    
    # Get names of outcomes from left-hand-side of formulas
    y_columns <- 
      sapply(
        X = outcome_formulas,
        FUN =
          function(x) all.vars(update(x, . ~ 0))
      )
    
    n_outcomes <- length(y_columns)
    
    # Get name of treatment indicator from left-hand-side of formula
    tx_column <-
      all.vars(update(propensity_score_formula, . ~ 0))
    
    
    # Get names of predictors from formulas
    x_columns <-
      unique(
        do.call(
          what = c,
          args = 
            lapply(
              X = c(propensity_score_formula,
                    outcome_formulas,
                    inverse_weight_formulas),
              FUN = function(x, data)
                names(get_all_vars(formula = x, data = data)),
              data = data
            )
        )
      )
    
    # Outcomes may be in inverse weight or outcome formulas: remove them.
    x_columns <-
      setdiff(
        x = x_columns,
        y = c(y_columns, tx_column)
      )
    
    
    if(!is.null(impute_formulas)) {
      imputation_x <-
        do.call(
          what = c,
          args = 
            lapply(
              X = impute_formulas,
              FUN = function(x, data)
                names(get_all_vars(formula = x, data = data)),
              data = data
            )
        )
      
      imputation_x <-
        setdiff(
          x = imputation_x,
          y = c(x_columns, tx_column, y_columns)
        )
    } else {
      imputation_x <- NULL
    }
    
    return(
      list(
        x_columns = x_columns,
        tx_column = tx_column,
        y_columns = y_columns,
        imputation_x = imputation_x
      )
    )
  }




### tmle_precheck ##############################################################
# Check arguments from TMLE for potential problems: missing values in baseline
# covariates or treatment, treatment variable is not binary, missing data that
# is not monotone and imputation arguments are not specified, or the imputation
# arguments are not correctly specified.
# If any issues are identified, the function halts with specific error messages,
# and returns the indices of intermediate outcomes to be imputed if no errors
# are found.
tmle_precheck <-
  function(
    data,
    x_columns,
    y_columns,
    tx_column,
    impute_covariates = FALSE,
    impute_monotone = NULL,
    imputation_x = NULL,
    impute_formulas = NULL,
    impute_model = NULL,
    inverse_weight_prediction_tx,
    outcome_type,
    outcome_prediction_tx,
    outcome_range = NULL,
    absorbing_state = NULL
  ) {
    
    # Check for missing covariates
    imputation_x <-
      setdiff(
        x = imputation_x,
        y = c(x_columns, tx_column, y_columns)
      )
    
    n_missing_covariates <-
      colSums(is.na(data[, c(x_columns, tx_column, imputation_x)]))
    
    
    if(sum(n_missing_covariates) > 0 & !impute_covariates) {
      n_missing_covariates <-
        n_missing_covariates[n_missing_covariates > 0]
      stop(
        "Missing covariate/treatment data: ",
        paste(
          paste0(names(n_missing_covariates),
                 " (", n_missing_covariates, ")"),
          collapse = ", "
        ),
        " - impute_covariates must be set to TRUE."
      )
    }
    
    if(!(inverse_weight_prediction_tx %in% c("counterfactual", "observed"))){
      stop("Only supported values of `inverse_weight_prediction_tx` are: ",
           "\"counterfactual\" and \"observed\".")
    }
    
    if(!(outcome_prediction_tx %in% c("counterfactual", "observed"))){
      stop("Only supported values of `outcome_prediction_tx` are: ",
           "\"counterfactual\" and \"observed\".")
    }
    
    # Check for treatment indicator that is not in {0, 1}
    if(sum(!(data[, tx_column] %in% 0:1)) > 0){
      stop(paste0("Treatment column `", tx_column, "` must be binary."))
    }
    
    # Check outcome specification vs. data
    if(outcome_type == "logistic"){
      if(is.null(outcome_range)){
        stop("Outcome range not specified for logistic outcome model:")
      }
      
      out_of_range <-
        colSums(x = data[, y_columns] > outcome_range[2], na.rm = TRUE) +
        colSums(x = data[, y_columns] < outcome_range[1], na.rm = TRUE)
      
      if(sum(out_of_range) > 0) {
        stop("Outcomes out of range [0, 1]: ",
             paste(names(which(out_of_range > 0)), collapse = ", "))
      }
    } else if (outcome_type == "binomial") {
      out_of_range <-
        apply(
          X = data[, y_columns],
          MARGIN = 2,
          FUN = function(x) sum(!x %in% c(0, 1, NA))
        )
      
      if(sum(out_of_range) > 0) {
        stop("Outcomes out of range {0, 1}: ",
             paste(names(which(out_of_range > 0)), collapse = ", "))
      }
    } else if(outcome_type == "multinomial-binomial"){
      if(
        !all(
          sapply(
            X = data[ , head(y_columns, -1)],
            FUN = class
          ) == "factor"
        )
      ){
        stop("Intermediate multinomial outcomes must be type 'factor'")
      }
      
      out_of_range <-
        !(data[, tail(y_columns, 1)] %in% c(0, 1, NA))
      
      if(sum(out_of_range) > 0 ){
        stop("`outcome_type` == 'multinomial-binomial' with final outcome ",
             "values outside {0, 1} in rows ",
             paste(which(out_of_range), collapse = ", "))
      }
    } else{
      if(!is.null(absorbing_state)){
        stop("`absorbing_state` specified without ",
             "`outcome_type` == 'multinomial-binomial'")
      }
    }
    
    non_monotone_outcomes <-
      non_monotone(
        data = data,
        y_columns = y_columns
      )
    
    # Check imputation parameters
    if(sum(non_monotone_outcomes) > 0){
      
      # Non-monotone missingness - No handling of missingness specified
      if(length(impute_monotone) != 1){
        stop("Missingness is not monotone: ",
             "`impute_monotone` must be either FALSE or TRUE.")
      } else if(!(impute_monotone %in% c(FALSE, TRUE))){
        stop("Missingness is not monotone: ",
             "`impute_monotone` must be FALSE or TRUE.")
      }
      
      if(impute_monotone){
        n_non_monotone <- colSums(non_monotone_outcomes)
        
        # Find which variables must be imputed: Ensure temporal order
        non_monotone_vars <- names(which(n_non_monotone > 0))
        
        imputation_rhs <- 
          sapply(
            X = impute_formulas,
            FUN =
              function(x) all.vars(update(x, . ~ 0))
          )
        
        imputation_rhs_missing <- 
          setdiff(x = non_monotone_vars, y = imputation_rhs)
        
        if(length(imputation_rhs_missing) > 0) {
          stop("No model specified for outcomes with non-monotone missingness: ",
               paste(imputation_rhs_missing, collapse = ", "))
        }
        
        if(is.null(impute_formulas)){
          stop("Missingness is not monotone and no imputation formula is supplied.")
        }
        
        if(!impute_model %in%
           c("gaussian", "binomial", "beta", "pmm", "multinomial")) {
          stop("Invalid specification for `impute_model`:",
               impute_model)
        }
        
        # Find which variables must be imputed: Ensure temporal order
        non_monotone_vars <- 
          names(which(colSums(non_monotone_outcomes) > 0))
        
        imputation_rhs <- 
          sapply(
            X = impute_formulas,
            FUN =
              function(x) all.vars(update(x, . ~ 0))
          )
        
        imputation_rhs_missing <- 
          setdiff(x = non_monotone_vars, y = imputation_rhs)
        
        if(length(imputation_rhs_missing) > 0) {
          stop("No model specified for outcomes with non-monotone missingness: ",
               paste(imputation_rhs_missing, collapse = ", "))
        }
        
        # Subset formulas to variables with non-monotone missingness
        # Order in temporal sequence if needed
        impute_formulas <-
          impute_formulas[
            match(x = intersect(x = imputation_rhs, y = non_monotone_vars),
                  table = y_columns)
          ]
        
        
      } else {
        impute_columns = NULL
        impute_formulas = NULL
      }
    }
    
    return(
      list(
        impute_columns = non_monotone_outcomes,
        impute_formulas = impute_formulas
      )
    )
  }




### tmle_get_formulas ##########################################################
# Check specified formulas: missingness formulas can be supplied as a
# right-hand-side formula - use y_columns to add appropriate outcome.
# Since constructing the TMLE algorithm requires constructing new variables,
# new variable names are added to the data.frame to avoid over-writing the
# original values (`..y1`, `..y2`, `..y3` instead of `y1`, `y2`, `y3`). The
# outcome formulas need to be adjusted accordingly.
tmle_get_formulas <-
  function(
    y_columns,
    inverse_weight_formulas,
    outcome_formulas
  ){
    outcome_formulas_tx_1 <- outcome_formulas_tx_0 <- list()
    
    for(i in 1:length(inverse_weight_formulas)){
      # Cut off LHS (if supplied) - use y_columns to find response
      inverse_weight_formulas[[i]] <-
        reformulate(
          termlabels =
            Reduce(paste, deparse(tail(inverse_weight_formulas[[i]], 1)[[1]])),
          response =
            paste("!is.na(", y_columns[i], ")", collapse = "", sep = "")
        )
    }

    for(i in 1:length(outcome_formulas)){
      outcome_formulas_tx_1[[i]] <-
        reformulate(
          termlabels =
            Reduce(paste, deparse(tail(outcome_formulas[[i]], 1)[[1]])),
          response =
            paste0("..", y_columns[i], "_tx_1")
        )
      
      outcome_formulas_tx_0[[i]] <-
        reformulate(
          termlabels =
            Reduce(paste, deparse(tail(outcome_formulas[[i]], 1)[[1]])),
          response =
            paste0("..", y_columns[i], "_tx_0")
        )
    }
    
    return(
      list(
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas_tx_1 = outcome_formulas_tx_1,
        outcome_formulas_tx_0 = outcome_formulas_tx_0
      )
    )
  }




### tmle_compute ###############################################################
# NOTE: this is an internal function, meant to be called by tmle(). It has been
# stripped down to speed up bootstrapping at the expense of error handling that
# is done by tmle().
tmle_compute <-
  function(
    data,
    y_columns, tx_column, x_columns,
    propensity_score_formula,
    inverse_weight_formulas,
    inverse_weight_stratified = FALSE,
    inverse_weight_prediction_tx = "counterfactual",
    outcome_formulas_tx_1,
    outcome_formulas_tx_0,
    outcome_type,
    estimand,
    outcome_prediction_tx = "counterfactual",
    outcome_range = NULL,
    absorbing_state = NULL,
    absorbing_outcome = NULL,
    impute_covariates = FALSE,
    impute_monotone = FALSE,
    impute_formulas = NULL,
    impute_model = NULL,
    imputation_args = NULL,
    verbose = FALSE,
    max_abs_weight = 20
  ) {
    
    if(verbose) data_original <- data
    
    n_outcomes <- length(y_columns)
  
    impute_columns <-
      non_monotone(data = data, y_columns = y_columns)
    
    if(impute_covariates){
      data <-
        impute_covariates_mean_mode(
          data = data, x_columns = x_columns
        )
    }
    
    if(impute_monotone){
      imputed_outcomes <-
        do.call(
          what = impute_outcomes_to_monotone,
          args =
            c(
              list(
                model = impute_model,
                data = data,
                impute_columns = impute_columns,
                impute_formulas = impute_formulas,
                absorbing_state = absorbing_state,
                y_columns = y_columns,
                verbose = verbose
              ),
              imputation_args
            )
        )
      
      if(verbose){
        data <- imputed_outcomes$imputed_data
      } else {
        data <- imputed_outcomes
      }
    } else {
      
      for(i in 1:ncol(impute_columns)){
        
        censor_outcome <- names(impute_columns)[i]
        
        data[which(impute_columns[, i]),
             y_columns[which(y_columns == censor_outcome):length(y_columns)]
        ] <- NA
      }
      
      imputed_outcomes <- NULL
    }
    
    # Fit propensity score model & compute propensity score
    propensity_model <-
      # mgcv::gam(
      stats::glm(
        formula = propensity_score_formula,
        data = data,
        family = binomial
      )
    
    pr_tx_1 <- propensity_model$fitted.values
    pr_tx_0 <- (1 - pr_tx_1)
    
    # Compute inverse probability of censoring weights
    inverse_weights <-
      compute_inverse_weights(
        data = data,
        inverse_weight_formulas = inverse_weight_formulas,
        y_columns = y_columns,
        tx_column = tx_column,
        inverse_weight_stratified = inverse_weight_stratified,
        inverse_weight_prediction_tx = inverse_weight_prediction_tx,
        absorbing_state = absorbing_state,
        verbose = verbose
      )
    
    ipw_tx_1 <- inverse_weights$ip_weights_tx_1
    ipw_tx_0 <- inverse_weights$ip_weights_tx_0
    
    ipw_tx_1 <-
      ipw_tx_1*kronecker(
        X = matrix(data = 1, ncol = n_outcomes),
        Y = 1/pr_tx_1
      )
    
    ipw_tx_0 <-
      ipw_tx_0*kronecker(
        X = matrix(data = 1, ncol = n_outcomes),
        Y = 1/pr_tx_0
      )
    
    # Truncate Weights
    ipw_tx_1 <-
      apply(
        X = ipw_tx_1,
        MARGIN = 2,
        FUN = function(x) pmin(x, max_abs_weight)
      )
    
    ipw_tx_0 <-
      apply(
        X = ipw_tx_0,
        MARGIN = 2,
        FUN = function(x) pmin(x, max_abs_weight)
      )
    
    if(any(colSums(!is.na(ipw_tx_1)) == 0)){
      stop("Error in IP weights columns: ",
           paste0(which(colSums(!is.na(ipw_tx_1)) == 0), collapse = ", "))
    }
    
    if(any(colSums(!is.na(ipw_tx_0)) == 0)){
      stop("Error in IP weights columns: ",
           paste0(which(colSums(!is.na(ipw_tx_0)) == 0), collapse = ", "))
    }
    
    # Compute sequence of regression fits:
    regression_sequence_tx_1 <- regression_sequence_tx_0 <- list()
    
    # Create indicators for those not previously censored or absorbed
    uncensored <-
      !is.na(data[, head(y_columns, -1)])
    
    absorbed <-
      t(
        apply(
          X = data[, head(y_columns, -1)],
          MARGIN = 1,
          FUN = function(x, state = absorbing_state)
            cummax(x %in% state)
        )
      )

    
    # Add IPW, uncensored indicators, outcome model fits 
    fit_y_columns_tx_1 <- paste0("..", y_columns, "_tx_1")
    fit_y_columns_tx_0 <- paste0("..", y_columns, "_tx_0")
    data[, c(fit_y_columns_tx_1, fit_y_columns_tx_0)] <- NA
    data[, c(tail(x = fit_y_columns_tx_1, 1),
             tail(x = fit_y_columns_tx_0, 1))] <-
      data[, tail(x = y_columns, 1)]
    
    data <-
      cbind(
        data, 
        setNames(object = data.frame(ipw_tx_1),
                 nm = paste0("..ipw_tx_1_", 1:ncol(ipw_tx_1))),
        setNames(object = data.frame(ipw_tx_0),
                 nm = paste0("..ipw_tx_0_", 1:ncol(ipw_tx_0))),
        setNames(object = data.frame(cbind(TRUE, uncensored & !absorbed)),
                 nm = paste0("..u_", 1:ncol(ipw_tx_1)))
    )
    
    if(outcome_type %in% "gaussian"){
      glm_family <- gaussian
    } else{
      glm_family <- quasibinomial
    }
    
    for(i in n_outcomes:1){
      uncensored <- which(data[, paste0("..u_", i)])
      
      outcome_regression_tx_1 <-
        eval(
          parse(
            text = 
              paste0(
                "stats::glm(formula = ", 
                Reduce(paste, deparse(outcome_formulas_tx_1[[ i ]])), ", ",
                "family = glm_family, ",
                "data = data, subset = ..u_", i, ", ",
                "weights = ..ipw_tx_1_", i, ")"
              )
          )
        )
      
      outcome_regression_tx_0 <-
        eval(
          parse(
            text = 
              paste0(
                "stats::glm(formula = ", 
                Reduce(paste, deparse(outcome_formulas_tx_0[[ i ]])), ", ",
                "family = glm_family, ",
                "data = data, subset = ..u_", i, ", ",
                "weights = ..ipw_tx_0_", i, ")"
              )
          )
        )
      
      if(verbose) regression_sequence_tx_1[[i]] <- outcome_regression_tx_1
      if(verbose) regression_sequence_tx_0[[i]] <- outcome_regression_tx_0
      
      if(i > 1) {
        data[uncensored, fit_y_columns_tx_1[(i - 1)]] <-
          predict(
            object = outcome_regression_tx_1,
            newdata = 
              within(
                data = subset(data, eval(parse(text = paste0("..u_", i)))),
                expr = {
                  if(outcome_prediction_tx == "counterfactual") {
                    eval(parse(text = paste0(tx_column, " = 1")))
                  }
                }
              ),
            type = "response"
          )
        
        data[uncensored, fit_y_columns_tx_0[(i - 1)]] <-
          predict(
            object = outcome_regression_tx_0,
            newdata = 
              within(
                data = subset(data, eval(parse(text = paste0("..u_", i)))),
                expr = {
                  if(outcome_prediction_tx == "counterfactual") {
                    eval(parse(text = paste0(tx_column, " = 0")))
                  }
                }
              ),
            type = "response"
          )
        
        if(!is.null(absorbing_state)){
          absorbed_i <- which(data[, y_columns[(i-1)]] %in% absorbing_state)
          if(length(absorbed_i) > 0){
            data[absorbed_i, fit_y_columns_tx_1[(i - 1)]] <- absorbing_outcome
            data[absorbed_i, fit_y_columns_tx_0[(i - 1)]] <- absorbing_outcome
          }
        }
      }
    }
    
    # Compute ATE
    y_tx_1 <-
      predict(
        object = outcome_regression_tx_1,
        newdata = 
          within(
            data = data,
            expr = {eval(parse(text = paste0(tx_column, " = 1")))}
          ),
        type = "response"
      )
    
    y_tx_0 <-
      predict(
        object = outcome_regression_tx_0,
        newdata = within(
          data = data,
          expr = {eval(parse(text = paste0(tx_column, " = 0")))}
        ),
        type = "response"
      )
    
    if(outcome_type == "logistic"){
      y_tx_1 <- 
        decompress_range(
          x = y_tx_1,
          min = outcome_range[1],
          max = outcome_range[2]
        )
      
      y_tx_0 <- 
        decompress_range(
          x = y_tx_0,
          min = outcome_range[1],
          max = outcome_range[2]
        )
    }
    
    if(estimand == "difference"){
      ate <- mean(y_tx_1) - mean(y_tx_0)
    } else if (estimand == "ratio"){
      ate <- mean(y_tx_1)/mean(y_tx_0)
    } else if (estimand == "oddsratio"){
      ate <- mean(y_tx_1)*(1 - mean(y_tx_0))/(mean(y_tx_0)*(1 - mean(y_tx_1)))
    }
    
    
    if(verbose){
      return(
        list(
          "ATE: Estimate" = ate,
          "E[Y|A=1]: Estimate" = mean(y_tx_1),
          "E[Y|A=0]: Estimate" = mean(y_tx_0),
          y_tx_1 = y_tx_1,
          y_tx_0 = y_tx_0,
          data = data,
          outcome_range = outcome_range,
          outcome_formulas_tx_1 = outcome_formulas_tx_1,
          outcome_formulas_tx_0 = outcome_formulas_tx_0,
          outcome_models_tx_1 = regression_sequence_tx_1,
          outcome_models_tx_0 = regression_sequence_tx_0,
          ipw_tx_1 = ipw_tx_1,
          ipw_tx_0 = ipw_tx_0,
          inverse_weights = inverse_weights,
          inverse_weight_formulas = inverse_weight_formulas,
          propensity_model = propensity_model,
          pr_tx_1 = pr_tx_1,
          pr_tx_0 = pr_tx_0,
          imputed_outcomes = imputed_outcomes,
          impute_formulas = impute_formulas,
          imputation_args = imputation_args,
          data_original = data_original
        )
      )
    } else {
      return(
        c("ATE: Estimate" = ate,
          "E[Y|A=1]: Estimate" = mean(y_tx_1),
          "E[Y|A=0]: Estimate" = mean(y_tx_0))
        )
    }
  }




### tmle_boot_wrap #############################################################
# wrapper to be used with boot()
tmle_boot_wrap <-
  function(
    data, indices = NULL, tmle_args
  ) {
  
    if(!is.null(indices)) {
      data <- data[indices, ]
    }

    do.call(
      what =
        function(x = data, ...) {
          tmle_compute(data = x, ...)
        },
      args = tmle_args
    )
  }





### tmle #######################################################################
rctmle <-
  function(
    data,
    propensity_score_formula,
    inverse_weight_formulas,
    inverse_weight_stratified = FALSE,
    inverse_weight_prediction_tx = "counterfactual",
    outcome_formulas,
    outcome_type,
    estimand = "difference",
    outcome_prediction_tx = "counterfactual",
    outcome_range = NULL,
    absorbing_state = NULL,
    absorbing_outcome = NULL,
    impute_covariates = FALSE,
    impute_monotone = FALSE,
    impute_formulas = NULL,
    impute_model = NULL,
    imputation_args = NULL,
    ci = FALSE,
    verbose = FALSE,
    bootstrap_n = 10000,
    bootstrap_type = c("bca", "norm", "basic", "perc")[1],
    alpha = 0.05,
    ...
  ) {
    
    if(!is.null(outcome_type)) outcome_type <- tolower(outcome_type)
    if(!is.null(impute_model)) impute_model <- tolower(impute_model)
    
    # arg_list <- as.list(substitute(list(...)))[-1L]
    arg_list <- list()
    
    if(ci) {
      if(!bootstrap_type %in% c("bca", "norm", "basic", "perc")) {
        stop("Unrecognized bootstrap method: ", bootstrap_type)
      }
    }
    
    tmle_args <-
      tmle_get_args(
        data = data,
        propensity_score_formula = propensity_score_formula,
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas = outcome_formulas,
        impute_formulas = impute_formulas
      )
    
    
    # When absorbing state is specified, check for consistency and carry forward
    # absorbing state values
    if(!is.null(absorbing_state)){
      checked_data <-
        absorbing_state_check(
          data = data,
          y_columns = tmle_args$y_columns,
          outcome_type = outcome_type,
          absorbing_state = absorbing_state,
          absorbing_outcome = absorbing_outcome
        )
      
      environment(checked_data) <-
        environment(data)
      
      data <- checked_data
    }
    
    precheck_results <-
      tmle_precheck(
        data = data,
        x_columns = tmle_args$x_columns,
        tx_column = tmle_args$tx_column,
        y_columns = tmle_args$y_columns,
        impute_covariates = impute_covariates,
        imputation_x = tmle_args$imputation_x,
        impute_formulas = impute_formulas,
        impute_model = impute_model,
        impute_monotone = impute_monotone,
        inverse_weight_prediction_tx = inverse_weight_prediction_tx,
        outcome_type = outcome_type,
        outcome_prediction_tx = outcome_prediction_tx,
        outcome_range = outcome_range,
        absorbing_state = absorbing_state
      )
    
    if(outcome_type == "logistic"){
      data[, tmle_args$y_columns] <-
        apply(
          X = data[, tmle_args$y_columns],
          MARGIN = 2,
          FUN = compress_range,
          min = outcome_range[1],
          max = outcome_range[2],
          enforce_range = TRUE
        )
    }
    
    # Construct left-hand-side of formulas
    tmle_formulas <-
      tmle_get_formulas(
        y_columns = tmle_args$y_columns,
        inverse_weight_formulas = inverse_weight_formulas,
        outcome_formulas = outcome_formulas
      )
    
    
    # Re-set environments for formulas updated within functions
    for (i in 1:length(tmle_formulas$inverse_weight_formulas))
      environment(tmle_formulas$inverse_weight_formulas[[i]]) <-
      environment(inverse_weight_formulas[[i]])
    
    for (i in 1:length(outcome_formulas)){
      environment(tmle_formulas$outcome_formulas_tx_1[[i]]) <-
        environment(outcome_formulas[[i]])
      
      environment(tmle_formulas$outcome_formulas_tx_0[[i]]) <-
        environment(outcome_formulas[[i]])
    }

    if(impute_monotone){    
      for (i in 1:length(precheck_results$impute_formulas))
        environment(precheck_results$impute_formulas[[i]]) <-
          environment(impute_formulas[[i]])
    }

    
    arg_list <-
      c(
        list(
          y_columns = tmle_args$y_columns,
          tx_column = tmle_args$tx_column,
          x_columns = tmle_args$x_columns,
          propensity_score_formula = propensity_score_formula,
          inverse_weight_formulas =
            tmle_formulas$inverse_weight_formulas,
          inverse_weight_stratified =
            inverse_weight_stratified,
          inverse_weight_prediction_tx = inverse_weight_prediction_tx,
          outcome_formulas_tx_1 =
            tmle_formulas$outcome_formulas_tx_1,
          outcome_formulas_tx_0 =
            tmle_formulas$outcome_formulas_tx_0,
          outcome_prediction_tx =
            outcome_prediction_tx,
          outcome_type = outcome_type,
          estimand = estimand,
          outcome_range = outcome_range,
          absorbing_state = absorbing_state,
          absorbing_outcome = absorbing_outcome,
          impute_covariates = impute_covariates,
          impute_monotone = impute_monotone,
          impute_formulas = impute_formulas,
          imputation_args = imputation_args,
          impute_model = impute_model,
          verbose = verbose
        ),
        # impute_model, imputation_args, ...
        arg_list
      )
    
    tmle_result <-
      tmle_boot_wrap(
        data = data,
        tmle_args = arg_list
      )
    
    if(ci){
      arg_list$verbose <- FALSE
      
      tmle_boot <- 
        boot(
          data = data,
          statistic = tmle_boot_wrap,
          R = bootstrap_n,
          tmle_args =
            arg_list
        )
      
      tmle_boot_se <- 
        apply(
          X = tmle_boot$t,
          MAR = 2,
          FUN = sd
        )
      
      tmle_boot_ci_ate <-
        boot.ci(
          boot.out = tmle_boot,
          conf  = (1 - alpha),
          type = bootstrap_type,
          index = 1
        )
      
      tmle_boot_ci_y0_a_1 <-
        boot.ci(
          boot.out = tmle_boot,
          conf  = (1 - alpha),
          type = bootstrap_type,
          index = 2
        )
      
      tmle_boot_ci_y0_a_0 <-
        boot.ci(
          boot.out = tmle_boot,
          conf  = (1 - alpha),
          type = bootstrap_type,
          index = 3
        )
      
      lcl_ucl_ate <-
        tail(x = tmle_boot_ci_ate[bootstrap_type][[1]][1,], 2)
      lcl_ucl_y0_a_1 <-
        tail(x = tmle_boot_ci_y0_a_1[bootstrap_type][[1]][1,], 2)
      lcl_ucl_y0_a_0 <-
        tail(x = tmle_boot_ci_y0_a_0[bootstrap_type][[1]][1,], 2)
      
      return_result <-
        c("ATE: Estimate" = as.numeric(tmle_result[1]),
          "ATE: SE" = tmle_boot_se[1],
          "ATE: LCL" = lcl_ucl_ate[1],
          "ATE: UCL" = lcl_ucl_ate[2],
          "E[Y|A=1]: Estimate" = as.numeric(tmle_result[2]),
          "E[Y|A=1]: SE" = tmle_boot_se[2],
          "E[Y|A=1]: LCL" = lcl_ucl_y0_a_1[1],
          "E[Y|A=1]: UCL" = lcl_ucl_y0_a_1[2],
          "E[Y|A=0]: Estimate" = as.numeric(tmle_result[3]),
          "E[Y|A=0]: SE" = tmle_boot_se[2],
          "E[Y|A=0]: LCL" = lcl_ucl_y0_a_0[1],
          "E[Y|A=0]: UCL" = lcl_ucl_y0_a_0[2]
        )
      
      if(verbose) {
        return_result <-
          c(
            return_result,
            list(
              tmle_boot = tmle_boot
            ),
            tmle_result[setdiff(x = names(tmle_result), y = names(return_result))]
          )
      }
      
      return(return_result)
    } else {
      return(tmle_result)
    }
  }


boot_p_value <-
  function(
    boot_object,
    ci_method = "bca",
    index = 1,
    alpha_max = 1,
    alpha_min = 10^-5,
    tolerance = 0.001,
    max_evaluations = 40,
    verbose = FALSE
  ){
    
    converged <- FALSE
    continue <- TRUE
    i <- n_rejected <- n_fail_to_reject <- 0
    
    all_ci_results <-
      data.frame(
        alpha = rep(NA, max_evaluations),
        lcl = NA,
        ucl = NA,
        rejected = NA
      )
    
    current_min <- alpha_min
    current_max <- alpha_max
    current_alpha <- mean(c(alpha_max, alpha_min))
    
    while(continue) {
      i <- i +1
      
      all_ci_results$alpha[i] <- current_alpha
      
      ci_result <-
        boot.ci(
          boot.out = boot_object,
          conf = 1 - current_alpha,
          type = ci_method,
          index = index
        )
      
      all_ci_results[i, c("lcl", "ucl")] <- 
        tail(ci_result[[ci_method]][1,], 2)
      
      rejected <- all_ci_results$rejected[i] <-
        with(all_ci_results, sign(lcl[i]) == sign(ucl[i]))
      
      if(rejected) {
        # Decrease Alpha
        n_rejected <- n_rejected + 1
        current_max <- min(c(current_alpha, current_max))
        current_alpha <-
          current_alpha - c(current_alpha - current_min)/2
        
      } else {
        # Increase Alpha
        n_fail_to_reject <- n_fail_to_reject + 1
        current_min <- max(c(current_alpha, current_min))
        current_alpha <-
          current_alpha + c(current_max - current_alpha)/2
      }
      
      
      # Evaluate convergence
      if(i > 1){
        diff_alpha <- abs(diff(all_ci_results$alpha[c(i-1, i)]))
        
        if(n_fail_to_reject > 4 & n_rejected > 4 & diff_alpha <= tolerance) {
          converged <- TRUE
          continue <- FALSE
        }
      }
      
      if(i >= max_evaluations){
        continue <- FALSE
        warning("Iteration limit reached without convergence.")
      }
    }
    
    boot_p_value <- min(subset(all_ci_results, rejected)$alpha)
    
    if(verbose){
      return(
        list(
          boot_p_value = boot_p_value,
          converged = converged,
          boot_object = boot_object,
          ci_method = ci_method,
          index = index,
          alpha_max = alpha_max,
          alpha_min = alpha_min,
          tolerance = tolerance,
          max_evaluations = 40,
          n_evaluations = i,
          all_ci_results = all_ci_results
        )
      )
    } else {
      return(boot_p_value)
    }
  }
