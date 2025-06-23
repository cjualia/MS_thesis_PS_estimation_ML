# Loading packages

library(caret)
library(randomForest)
library(gbm)
library(rpart)
library(MatchIt)
library(ipw)
library(cobalt)
library(ggplot2)
library(patchwork)
library(survey)

#-------------------------------------------------------------------------------
# SIMULATION

# Generating simulated data for a specific scenario.
# This function defines scenarios A-E, but user can define their own just by 
# changing coeficients, covariables, and formulas for treatment and outcome 
# modeling. 
simulation <- function(size, scenario){
  
  # covariates
  w1 <- rnorm(size,0,1) 
  w2 <- rlnorm(size,0,0.5)
  w3 <- rnorm(size,0,10)
  w4 <- rbinom(size,1,0.5)
  w5 <- rbinom(size,1,0.2)
  w6 <- sample(1:5, size, replace=T, prob = c(0.5,0.3,0.1,0.05,0.05)) 
 
  # complex terms and interactions
  w7 <- sin(w1)
  w8 <- w2^2
  w9 <- w3*w4
  w10 <- w4*w5
  
  # scenarios 
  
  if (scenario == "A"){
    # true propensity score: probability of assignment of treatment
    tps <- (1 + exp(-(-3.5 + 1*w1 + 1*w2 + 0.1*w3 + 2*w4 + 2*w5 + 0.4*w6)))^ -1
    # outcome
    Y0 <- -5 + 0.5*w1 + 0.5*w2 + 0.05*w3 + 1*w4 + 1*w5 + 0.2*w6 
  }
  if (scenario == "B"){
    tps <- (1 + exp(-(-3.5 + 1*w1 + 1*w2 + 0.1*w3 + 2*w4 + 2*w5 + 0.8*w6 + 1*w10)))^ -1
    Y0 <- -3.5 + 0.4*w1 + 0.03*w2 + 0.03*w3 + 0.75*w4 + 0.75*w5 + 0.2*w6 + 0.4*w7
    + 0.02*w8 + 0.04*w9 + 0.5*w10
  }
  if (scenario == "C"){
    tps <- (1 + exp(-(-1 + 0.8*w1 + 0.06*w2 + 0.06*w3 + 1.5*w4 + 1.5*w5 + 0.4*w6 
                      + 0.8*w7 + 0.04*w8 + 0.08*w9 + 1*w10)))^ -1
    Y0 <- -3.5 + 0.4*w1 + 0.03*w2 + 0.03*w3 + 0.75*w4 * 0.75*w5 + 0.2*w6 + 0.4*w7
    + 0.02*w8 + 0.04*w9 + 0.5*w10
  }
  if (scenario == "D"){
    tps <- (1 + exp(-(-1.3 + 0.8*w1 + 0.06*w2 + 0.06*w3 + 1.5*w4 + 1.5*w5 + 0.4*w6 
                      + 0.8*w7 + 0.04*w8 + 0.08*w9 + 1*w10)))^ -1
    Y0 <- -3.7 + 0.4*w1 + 0.03*w3 + 0.75*w4 * 0.75*w5 + 0.2*w6 + 0.4*w7
    + 0.04*w9 + 0.5*w10
  }
  if (scenario == "E"){
    tps <- (1 + exp(-(-1.3 + 0.8*w1 + 0.06*w2 + 0.06*w3 + 1.5*w4 + 1.5*w5 + 0.4*w6 
                      + 0.8*w7 + 0.04*w8 + 0.08*w9 + 1*w10)))^ -1
    Y0 <- -3.7 + 0.4*w1 + 0.03*w3 + 0.75*w4 * 0.75*w5 + 0.2*w6 + 1*w7
    + 2*w9 + 1.5*w10
  }
  
  # counterfactual Y1 corresponds to the outcome when treated,
  # therefore, the effect of the treatment is added (here it's set as constant, -1):
  effect <- -1
  Y1 <- Y0 + effect
  
  # treatment assignment:
  unif1 <- runif(size,0,1) # Uniform continuous distribution (0,1)
  Tr <- ifelse(tps > unif1, 1, 0) # When PS > number, Tr=1 (treated); else, Tr=0 (non-treated)
  
  # continuous outcome Y
  Y    <- Tr*Y1 + (1-Tr)*Y0
  
  # simulated dataset
  sim <- as.data.frame(cbind(w1, w2, w3 ,w4, w5, w6, w7, w8, w9, w10, Tr, Y, tps))
  
  sim <- sim[sim$tps >= 0.05 & sim$tps <= 0.95, ]
  
  return(sim)
}

# Monte Carlo simulations of each scenario
# executes `simulation` nrep times. nrep is the number of Monte Carlo simulations
rep_simulation <- function(nrep, size, scenario){
  all_sim <- list()
  for (r in seq(nrep)){
    df <- simulation(size, scenario)
    all_sim[[r]] <- df
  }
  return(all_sim)
}

#-------------------------------------------------------------------------------
# INITIAL BALANCE

# simple function to calculate the initial ASMD of non-adjusted populations
# (before the balance)
before <- function(scenario){
  
  dfs <- all_scenarios[[scenario]]
  
  results <- list()
  
  for (i in 1:length(dfs)){ # iteration over the list of datasets
    
    m.out0 <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
                      data = dfs[[i]],
                      method = NULL)  # no method 
    
    summary <- summary(m.out0)
    
    asmd <- summary$sum.all[-1,3] # extract ASMD's (ignore row `distance`)
    
    results[[i]] <- asmd # save result for each dataframe
  }
  
  return(results)
}

#-------------------------------------------------------------------------------
# ESTIMATION OF PS BY LOGISTIC REGRESSION

# The following functions return a list in which each element correspond to a
# dataset of the same scenario, which contains the vector of ASMD calculated
# for each dataset after the balance.

# Estimates PS by LR and applies the estimated PS using matching 
RL_matching <- function(scenario){ 
  
  dfs <- all_scenarios[[scenario]]
  
  results <- list()
  
  for (i in 1:length(dfs)){ # iterate over list of dataframes
    
    # 2 steps:
    
    # 1) `matchit` estimates PS, y applies matching
    
    m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6, # omitting complex variables
                      data = dfs[[i]],
                      method = "nearest", # Matching
                      replace = TRUE, 
                      caliper = 0.1, 
                      distance = "glm") # estimation method: logistic regresion
    
    # 2) `bal.tab` calculate ASMD's
    # we add the covariates that were not specified in the estimation formula
    balance <- bal.tab(m.out,
                       addl = ~ w7 + w8 + w9 + w10,
                       data = dfs[[i]])
    
    # extraction of ASDM's
    asmds <- abs(balance$Balance$Diff.Adj)[2:11] # Omitimos registro de `distance`
    results[[i]] <- asmds
  }
  return(results)
}

# Estimates PS by LR and applies the estimated PS using IPW 
RL_weight <- function(scenario) {
  
  dfs <- all_scenarios[[scenario]]
  results <- list()
  
  for (i in seq_along(dfs)) {
    
    # 3 steps:
    
    # 1) `matchit` estimates PS
    
    m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6, # omitting complex variables
                     data = dfs[[i]],
                     distance = "glm") # estimation method: logistic regression
    
    eps <- m.out$distance
    
    # 2) calculation of stabilized weights
    
    p_Tr <- mean(dfs[[i]]$Tr == 1)
    p_noTr <- 1 - p_Tr
    weights <- ifelse(dfs[[i]]$Tr == 1,
                      p_Tr / eps,
                      p_noTr / (1 - eps))
    
    # 2. `bal.tab` performs balance
    balance <- bal.tab(
      Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
      data = dfs[[i]],
      weights = weights,
      method = "weighting")
    
    asmds <- abs(balance$Balance$Diff.Adj)
    results[[i]] <- asmds
  }
  
  return(results)
}

#-------------------------------------------------------------------------------
# ESTIMATION OF PS BY MACHINE LEARNING



# HYPERPARAMETERS TESTING

# Returns estimated PS for a single dataset (df), using the specified ML method
# (method), adjusted with a specific value of the hyperparameter (hp)

# NOTE: In this simple implementation, only 1 hyperparameter for each ML method
# is tested! For CART, cp; for Random Forest, mtry; and for Gradient Boosting, 
# shrinkage.

ML_eps <- function(df, method, hp){ 
  
  if (method == "randomForest"){
    
    rf <- randomForest(factor(Tr) ~ w1 + w2 + w3 + w4 + w5 + w6,
                       data = df,
                       mtry = hp)
    eps <- predict(rf, type = "prob")[, 2]
    
  } else if (method == "gbm"){
    gbm <- gbm(
      formula = Tr ~ w1 + w2 + w3 + w4 + w5 + w6,
      data = df,
      distribution = "bernoulli",
      n.trees = 100,
      interaction.depth = 3,
      shrinkage = hp,
      n.minobsinnode = 10, 
      verbose= FALSE)
    eps <- predict(gbm, newdata = df, n.trees = 100, type = "response")
    
  } else if (method == "rpart"){
    cart <- rpart(
      formula = Tr ~ w1 + w2 + w3 + w4 + w5 + w6,
      data = df,
      method = "class",
      cp = rpart.control(hp))
    eps <- predict(cart, newdata = df, type = "prob")[, 2]
  }
  
  return(eps)
}

# Matching using the estimated PS (eps). Returns the ASMD after adjustment
matching <- function(df, eps){
  m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,      
                   data = df,
                   method = "nearest",
                   replace = TRUE,
                   caliper = 0.1,
                   distance = eps)
  
  balance <- bal.tab(m.out,
                     data = df)
  
  asmds <-  abs(balance$Balance$Diff.Adj)[2:11]
  
  return(asmds)
}

# IPW using the estimated PS (eps). Returns the ASMD after adjustment
ipw <- function(df, eps){
  
  eps <- pmin(pmax(eps, 0.01), 0.99)
  
  # Cálculo de los weights 
  p_Tr <- mean(df$Tr == 1)
  p_noTr <- 1 - p_Tr
  weights <- ifelse(df$Tr == 1,
                    p_Tr / eps,
                    p_noTr / (1 - eps))
  
  # `bal.tab` calcula el balance
  balance <- bal.tab(
    Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
    data = df,
    weights = weights,
    method = "weighting")
  
  asmds <- abs(balance$Balance$Diff.Adj)
  
  return(asmds)
} 

# integrates the all above functions: 
# 1) Estimates PS for all datasets using the ML method specified (method), adjusted
#    with a determined hyperparameter (hp)
# 2) Applies the estimated PS for each dataset using matching.

hp_testing_m <- function(hp, scenario, method){
  
  dfs <- all_scenarios[[scenario]]
  results <- data.frame()
  
  for (i in seq_along(dfs)) {
    # for each dataframe of the scenario
    eps <- ML_eps(dfs[[i]], method, hp) # Estimation of PS
    asmds <- matching(dfs[[i]], eps) # Matching using estimated PS
    results <- rbind(results, asmds) 
  }
  
  return(results)  
}

# Same, but in step 2) applies IPW
hp_testing_w <- function(hp, scenario, method){
  
  dfs <- all_scenarios[[scenario]]
  results <- data.frame()
  
  for (i in seq_along(dfs)) {
    eps <- ML_eps(dfs[[i]], method, hp)
    asmds <- ipw(dfs[[i]], eps)
    results <- rbind(results, asmds)
  }
  
  return(results)  
} 



# SELECTION OF BEST HYPERPARAMETER 

# helps to select the best hyperparameter for each scenario, based on the
# minimum overall balance (ASAM).
opt_hp <- function(hp_tests, scenario){
  
  # to calculate the ASAM, we use different covariates in each scenario,
  # since each one includes different confounding variables.
  # If the scenario is A, we use w1 to w6
  
  if (scenario == "A"){
    confounders <- c(1:6)
  } else if (scenario == "B"){
    confounders <- c(1:6, 10)
  } else if (scenario == "C"){
    confounders <- c(1:10)
  } else{
    confounders <- c(1, 3:7, 9, 10) 
  }
  
  v <- c()
  for (test in hp_tests){
    asam_mean <- mean(apply(test[confounders], 1, mean))
    v <- c(v, asam_mean)
  }
  return(which.min(v))
}


# similar to above function. Summarizes the results for a specific scenario.
# Returns a list with the average ASMD's and ASAM's of all datasets.
balance <- function(df, scenario){
  
  cov <- paste0("w", c(1:10))
  
  # To calculate the ASAM, we use different covariates in each scenario,
  # since each scenario includes different confounding variables.
 
  if (scenario == "A"){  # for scenario A, we use w1-w6
    confounders <- c(1:6)
  } else if (scenario == "B"){ # for scenario B, we use w1-w6 w1-w6 y w10
    confounders <- c(1:6, 10)
  } else if (scenario == "C"){ # for scenario C, we use w1-w10 
    confounders <- c(1:10)
  } else{
    confounders <- c(1, 3:7, 9, 10) # for D and E, we use w1, w3-w7, w9 y w10
  }
  
  asmd_mean <- apply(df, 2, mean)
  names(asmd_mean) <- cov
  
  asmd_sd <- apply(df, 2, sd)
  names(asmd_sd) <- cov
  
  asam_mean <- mean(asmd_mean[confounders])
  asam_sd <- sd(apply(df[confounders], 1, mean))
  
  result <- list("asmd_mean" = asmd_mean,
                 "asmd_sd" = asmd_sd,
                 "asam_mean" = asam_mean,
                 "asam_sd" = asam_sd)
  return(result)
}



#------------------------------------------------------------------------------
# ADDITIONALS 

# ESTIMATION OF PS FOR A SINGLE DATASET (instead of all datasets of an scenario)
# functions to estimate PS of a single dataset.
simple_RL <- function(df){
  m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6, 
                   data = df,
                   distance = "glm")
  eps <- m.out$distance
  return(eps)
}
simple_RF <- function(df, hp){
  rf <- randomForest(factor(Tr) ~ w1 + w2 + w3 + w4 + w5 + w6,
                     data = df,
                     mtry = hp)
  eps <- predict(rf, type = "prob")[, 2]
  return(eps)
}
simple_GBM <- function(df, hp){
  gbm <- gbm(
    formula = Tr ~ w1 + w2 + w3 + w4 + w5 + w6,
    data = df,
    distribution = "bernoulli",
    n.trees = 100,
    interaction.depth = 3,
    shrinkage = hp,
    n.minobsinnode = 10, 
    verbose= FALSE)
  eps <- predict(gbm, newdata = df, n.trees = 100, type = "response")
  return(eps)
}
simple_CART <- function(df, hp){
  cart <- rpart(
    formula = Tr ~ w1 + w2 + w3 + w4 + w5 + w6,
    data = df,
    method = "class",
    cp = rpart.control(hp))
  eps <- predict(cart, newdata = df, type = "prob")[, 2]
  return(eps)
}

# ESTIMATION OF ATE FOR A SINGLE DATASET.
# Estimates PS by indicated method (method) adjusted by hyperparamenter (hp)

# NOTE: Is designed to apply PS only with IPW

ate_est <- function(df, method, hp){
  
  if (method == "randomForest"){
    eps = simple_RF(df, hp)
  } else if (method == "gbm"){
    eps = simple_GBM(df, hp)
  } else if (method == "cart"){
    eps = simple_CART(df, hp)
  } else if (method == "logreg"){
    eps= simple_RL(df)
  }
  
  # IPW
  
  eps <- pmin(pmax(eps, 0.01), 0.99)
  
  # Cálculo de los weights 
  p_Tr <- mean(df$Tr == 1)
  p_noTr <- 1 - p_Tr
  weights <- ifelse(df$Tr == 1,
                    p_Tr / eps,
                    p_noTr / (1 - eps))
  
  d <- svydesign(ids = ~1, data = df, weights = ~weights)
  results <- svyglm(Y ~ Tr, design = d)
  
  ate_est <- coef(results)["Tr"]
  
  ate <- mean(df[df$Tr == 1]$Y - df[df$Tr == 0]$Y) 
  
  return(ate)
}


