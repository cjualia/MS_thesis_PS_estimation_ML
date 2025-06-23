source("packages.R")

# Carga de datos

load("results/all_results.rda")
load("results/simple_estimations.rda")

#-------------------------------------------------------------------------------
# Distribución de PS

eps_distribution <- function(scenario) {
  Tr = rep(simple_estimations[[scenario]]$rf$Tr, times = 4)
  Tr <- factor(Tr, labels = c("No tratado", "Tratado"))
  
  eps = c(simple_estimations[[scenario]]$logreg$eps,
          simple_estimations[[scenario]]$rf$eps,
          simple_estimations[[scenario]]$gbm$eps,
          simple_estimations[[scenario]]$cart$eps)
  
  method = rep(c("RL", "RF", "GBM", "CART"),
               each = length(simple_estimations[[scenario]]$rf$eps))
  
  df1 <- data.frame(Tr = Tr, eps = eps, method = method, Tr = rep(Tr, times = 4))
  
  p1 <- ggplot(df1, aes(x = eps, fill = Tr)) +
    geom_density(alpha = 0.5) +
    facet_wrap(~ method, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = c("No tratado" = "#BF3EFF", "Tratado" = "#90EE90")) +
    labs(
      x = "Propensity score estimado",
      y = "Densidad",
      fill = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text = element_text(size = 10),
      legend.position = "bottom")
  
  tps <- simple_estimations[[scenario]]$rf$tps
  
  df2 <- data.frame(
    tps = rep(tps, times = 4),
    eps = c(
      simple_estimations[[scenario]]$logreg$eps,
      simple_estimations[[scenario]]$rf$eps,
      simple_estimations[[scenario]]$gbm$eps,
      simple_estimations[[scenario]]$cart$eps),
    model = rep(c("RL", "RF", "GBM", "CART"), each = length(tps)),
    Tr = Tr
  )
  
  p2 <- ggplot(df2, aes(x = tps, y = eps, color = Tr)) +
    geom_point(alpha = 0.3, size = 1.5) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
    facet_wrap(~ model, ncol = 2) +
    scale_color_manual(values = c("No tratado" = "#BF3EFF", "Tratado" = "#90EE90")) +
    labs(
      x = "Propensity Scores verdaderos",
      y = "Propensity Scores estimados",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      strip.text = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      legend.position = "bottom"
    )
  
  return(list("distribution" = p1, "estimation" = p2))
}

# Sólo para el caso base (A)
eps_distribution("A")$estimation
eps_distribution("B")$estimation
eps_distribution("C")$estimation

#-------------------------------------------------------------------------------
# Barplots ASAM

barplots <- function(scenario){
  
  if (scenario == "A"){
    results_m <- results_A_m 
    results_ipw <- results_A_w
  } else if (scenario == "B"){
    results_m <- results_B_m 
    results_ipw <- results_B_w
  } else if (scenario == "C"){
    results_m <- results_C_m 
    results_ipw <- results_C_w
  } else if (scenario == "D"){
    results_m <- results_D_m 
    results_ipw <- results_D_w
  } else if (scenario == "E"){
    results_m <- results_E_m 
    results_ipw <- results_E_w
  } 
  
  df <- data.frame(
    metodo = rep(c("RL", "RF", "GBM", "CART"), times = 2),
    tipo = rep(c("IPW", "Matching"), each = 4),
    asam = c(results_ipw$logr$asam_mean, results_ipw$rf$asam_mean, 
             results_ipw$gbm$asam_mean, results_ipw$cart$asam_mean,
             results_m$logr$asam_mean, results_m$rf$asam_mean, 
             results_m$gbm$asam_mean, results_m$cart$asam_mean),
    se = c(results_ipw$logr$asam_sd, results_ipw$rf$asam_sd, 
           results_ipw$gbm$asam_sd, results_ipw$cart$asam_sd,
           results_m$logr$asam_sd, results_m$rf$asam_sd, 
           results_m$gbm$asam_sd, results_m$cart$asam_sd)/sqrt(100))
  
  plot <- ggplot(df, aes(x = metodo, y = asam, fill = tipo)) +
    
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    
    geom_errorbar(aes(ymin = asam - se, ymax = asam + se),
                  position = position_dodge(width = 0.7),
                  width = 0.2) +
    
    labs(x = NULL, y = "ASAM", fill = NULL,
         title = paste("Escenario", scenario)) +
    
    scale_fill_manual(values = c("IPW" = "#98F5FF", "Matching" = "#5F9EA0")) +
    
    coord_cartesian(ylim = c(0, 0.25)) +  # <<--- Aquí se fija el eje Y
    
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 9),
          legend.text = element_text(size = 8),
          legend.title = element_blank(),
          plot.title = element_text(size = 10, hjust = 0.5))
  
  return(plot)
}


barplots_results <- list("A" <- barplots("A"),
                         "B" <- barplots("B"),
                         "C" <- barplots("C"),
                         "D" <- barplots("D"),
                         "E" <- barplots("E"))

#-------------------------------------------------------------------------------

# Loveplots

loveplot <- function(scenario, method){
  if (scenario == "A"){
    
    confounders <- c(1:6)
    if (method == "matching"){results <- results_A_m} 
    else if (method == "ipw") {results <- results_A_w}
    
  } else if (scenario == "B"){
    confounders <- c(1:6, 10)
    
    if (method == "matching"){results <- results_B_m} 
    else if (method == "ipw") {results <- results_B_w}
    
  } else if (scenario == "C"){
    confounders <- c(1:10)
    if (method == "matching"){results <- results_C_m} 
    else if (method == "ipw") {results <- results_C_w}
    
  } else if (scenario == "D"){
    confounders <- c(1, 3:7, 9, 10)
    if (method == "matching"){results <- results_D_m} 
    else if (method == "ipw") {results <- results_D_w}
    
  } else if (scenario == "E"){
    confounders <- c(1, 3:7, 9, 10)
    if (method == "matching"){results <- results_E_m} 
    else if (method == "ipw") {results <- results_E_w}
  }
  
  asmd_means <- c(results$logr$asmd_mean[confounders], 
                  results$rf$asmd_mean[confounders],
                  results$gbm$asmd_mean[confounders], 
                  results$cart$asmd_mean[confounders])
  
  asmd_se <- c(results$logr$asmd_sd[confounders], 
                results$rf$asmd_sd[confounders],
                results$gbm$asmd_sd[confounders], 
                results$cart$asmd_sd[confounders])/sqrt(100)
  
  model <- rep(c("RL","RF", "GBM", "CART"), each = length(confounders))
  
  cov <- rep(paste0("w", confounders), times = 4)
  
  df <- data.frame(as.factor(model),as.factor(cov),asmd_means, asmd_se)
  
  plot <- ggplot(df, aes(x = asmd_means, y = reorder(cov, asmd_means), color = model)) +
    geom_point(position = position_dodge(width = 0.3), size = 2) +
    geom_errorbarh(
      aes(xmin = asmd_means - asmd_se, xmax = asmd_means + asmd_se),
      height = 0.2,
      position = position_dodge(width = 0.5)) +
    
    scale_color_manual(values = c("RL" = "#0000FF", "RF" = "#FF3E96", 
                                  "GBM" = "#00FF00","CART" = "#00F5FF")) + 
    
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "black") +
    coord_cartesian(xlim = c(0, 0.4)) +
    labs(x = "ASMD", y = "Confusores", color = "Método") +
    #facet_wrap(~ escenario) +
    theme_minimal() +
    theme(legend.position = "bottom",
          axis.title = element_text(size = 10),
          legend.text = element_text(size = 10),
          legend.title = element_blank())
  
  return(plot)
} 

loveplot_results <- list("A_m" = loveplot_A_m <- loveplot("A", "matching"),
                         "B_m" = loveplot_B_m <- loveplot("B", "matching"),
                         "C_m" =loveplot_C_m <- loveplot("C", "matching"),
                         "D_m" =loveplot_D_m <- loveplot("D", "matching"),
                         "E_m" =loveplot_E_m <- loveplot("E", "matching"),
                         
                         "A_w" =loveplot_A_w <- loveplot("A", "ipw"),
                         "B_w" =loveplot_B_w <- loveplot("B", "ipw"),
                         "C_w" =loveplot_C_w <- loveplot("C", "ipw"),
                         "D_w" =loveplot_D_w <- loveplot("D", "ipw"),
                         "E_w" =loveplot_E_w <- loveplot("E", "ipw"))

#-------------------------------------------------------------------------------
Fig_barplots <- loveplot_results[1:5][[1]] + loveplot_results[6:10][[1]] + loveplot_results[1:5][[2]] + 
  loveplot_results[6:10][[2]] + loveplot_results[1:5][[3]] + loveplot_results[6:10][[3]] +
  plot_layout(ncol = 2)

Fig_balance <- barplots_results[[1]] + barplots_results[[2]] + barplots_results[[3]] + plot_layout(ncol = 3)
