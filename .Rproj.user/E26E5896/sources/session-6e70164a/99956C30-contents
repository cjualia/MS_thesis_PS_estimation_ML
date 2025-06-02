
# Generación de datos simulados
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
    tps <- (1 + exp(-(-3.5 + 1*w1 + 1*w2 + 0.1*w3 + 2*w4 + 2*w5 + 0.4*w6)))^ -1
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
  
  # Counterfactual Y1 corresponds to the outcome when treated,
  # therefore, the effect of the treatment is added (constant):
  effect <- -1
  Y1 <- Y0 + effect
  
  # Treatment assignment:
  unif1 <- runif(size,0,1) # Uniform continuous distribution (0,1)
  Tr <- ifelse(tps > unif1, 1, 0) # When PS > number, Tr=1 (treated); else, Tr=0 (non-treated)
  
  # continuous outcome Y
  Y    <- Tr*Y1 + (1-Tr)*Y0
  
  sim <- as.data.frame(cbind(w1, w2, w3 ,w4, w5, w6, w7, w8, w9, w10, Tr, Y, tps))
  
  sim <- sim[sim$tps >= 0.05 & sim$tps <= 0.95, ]
  
  return(sim)
}

rep_simulation <- function(nrep, size, scenario){
  all_sim <- list()
  for (r in seq(nrep)){
    df <- simulation(size, scenario)
    all_sim[[r]] <- df
  }
  return(all_sim)
}

# Balance de covariables inicial
before <- function(scenario){
  
  dfs <- all_scenarios[[scenario]]
  
  results <- list()
  
  for (i in 1:length(dfs)){ # Recorremos la lista de dataframes
    
    m.out0 <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,
                      data = dfs[[i]],
                      method = NULL) # No aplicamos ningún método)
    
    summary <- summary(m.out0)
    
    asmd <- summary$sum.all[-1,3] # Extraemos los ASMD's (ignoramos la línea `distance`)
    
    results[[i]] <- asmd # Guardamos el resultado para cada dataframe
  }
  
  return(results)
}

# Regresión Logística: Estimación de Propensity Scores + Matching
RL_matching <- function(scenario){ 
  
  dfs <- all_scenarios[[scenario]]
  
  results <- list()
  
  for (i in 1:length(dfs)){ # Recorremos la lista de dataframes
    
    # Se realiza en 2 pasos:
    
    # 1) `matchit` estima los PS, y realiza el matching
    
    m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6, # Omitimos las variables complejas
                      data = dfs[[i]],
                      method = "nearest", # Matching
                      replace = TRUE, 
                      caliper = 0.1, 
                      distance = "glm") # Método de estimación: regresión logística
    
    # 2) `bal.tab` calcula los ASMD's
    ## En este caso especificamos todas las variables porque nos interesa
    ## conocer el ASMD's de todas
    balance <- bal.tab(m.out,
                       addl = ~ w7 + w8 + w9 + w10,
                       data = dfs[[i]])
    
    # Extraemos los ASMD'S:
    asmds <- abs(balance$Balance$Diff.Adj)[2:11] # Omitimos registro de `distance`
    results[[i]] <- asmds
  }
  return(results)
}

# Regresión Logística: Estimación de Propensity Scores + IPW
RL_weight <- function(scenario) {
  
  dfs <- all_scenarios[[scenario]]
  results <- list()
  
  for (i in seq_along(dfs)) {
    
    # Se realiza en 3 pasos:
    
    # 1) `matchit` estima los PS
    
    m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6, # Omitimos las variables complejas
                     data = dfs[[i]],
                     distance = "glm") # Método de estimación: regresión logística
    
    eps <- m.out$distance
    
    # 2) Cálculo de weights
    
    p_Tr <- mean(dfs[[i]]$Tr == 1)
    p_noTr <- 1 - p_Tr
    weights <- ifelse(dfs[[i]]$Tr == 1,
                      p_Tr / eps,
                      p_noTr / (1 - eps))
    
    # 2. `bal.tab` calcula el balance
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

# Búsqueda de hiperparámetros óptimos (basado en el ASAM mínimo)

# `ML_eps` devuelve los PS estimados
ML_eps <- function(df, method, hp){ # Devuelve eps para 1 df con el hp especificado
  
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

# `matching` aplica el matching con los PS estimados y devuelve los ASMD's
matching <- function(df, eps){
  m.out <- matchit(Tr ~ w1 + w2 + w3 + w4 + w5 + w6 + w7 + w8 + w9 + w10,      
                   data = df,
                   method = "nearest",
                   replace = TRUE,
                   caliper = 0.1,
                   distance = eps)
  
  balance <- bal.tab(m.out,
                     addl = ~ w7 + w8 + w9 + w10,
                     data = df)
  
  asmds <-  abs(balance$Balance$Diff.Adj)[2:11]
  
  return(asmds)
}

# `ipw` aplica el IPW con los PS estimados y devuelve los ASMD's
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

# `hp_testing_m` integra las funciones anteriores para devolver los ASMD's
# resultantes de aplicar matching con los PS estimados por el método de ML especificado,
# ajustándolo con un determinado hiperparámetro (`hp`).
hp_testing_m <- function(hp, scenario, method){
  
  dfs <- all_scenarios[[scenario]]
  results <- data.frame()
  
  for (i in seq_along(dfs)) {
    # Para cada dataframe del escenario:
    eps <- ML_eps(dfs[[i]], method, hp) # Se estiman los eps por el ML especificado
    asmds <- matching(dfs[[i]], eps) # Se aplica matching, y se obtienen los ASMD's
    results <- rbind(results, asmds) 
  }
  
  return(results)  
}

# `hp_testing_w` hace lo mismo pero aplicando los PS estimados con IPW.
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

# `opt_hp` ayuda a seleccionar el hiperparámetro óptimo.
opt_hp <- function(hp_tests, scenario){
  
  # Para calcular el ASAM utilizamos distintas covariables en cada escenario; 
  # ya que cada uno presenta diferentes variables confusoras.
  # Si el escenario es A, utilizamos w1-w6
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


# `balance_results` muy similar a la función anterior. Simplemente resume los
# resultados para un determinado escenario. Devuelve una lista con el promedio
# y la desviación del ASAM, así como de los ASMD's.

balance <- function(df, scenario){
  
  cov <- paste0("w", c(1:10))
  
  # Para calcular el ASAM utilizamos distintas covariables en cada escenario; 
  # ya que cada uno presenta diferentes variables confusoras.
 
  if (scenario == "A"){  # Si el escenario es A, utilizamos w1-w6
    confounders <- c(1:6)
  } else if (scenario == "B"){ # Para B, utilizamos w1-w6 y w10
    confounders <- c(1:6, 10)
  } else if (scenario == "C"){ # Para C, utilizamos w1-w10 
    confounders <- c(1:10)
  } else{
    confounders <- c(1, 3:7, 9, 10) # Para D y E, utilizamos w1, w3-w7, w9 y w10
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

# Funciones de estimación de PS simples:
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

# `ate_estimation` devuelve la estimación del ATE para un dataframe. Estima los
# PS por el método indicado (`method`), y aplica IPW (no matching). 
# Devuelve la estimación del ATE calculada con el paquete `survey`.

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




