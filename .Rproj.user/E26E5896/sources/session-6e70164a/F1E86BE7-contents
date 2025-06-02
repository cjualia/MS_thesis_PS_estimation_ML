

source("packages.R")
source("functions.R")

# Carga de todos los datos simulados:

load("data/all_scenarios.rda")

#-------------------------------------------------------------------------------

# 1. Balance de covariables inicial

## La función `before` devuelve una lista en la que cada elemento guarda el vector
## de ASMD's para covariable (en total, 10) de cada dataframe (en total, 100). 
## Convertimos la lista en una matriz 10 x 100 (dataframes x covariables) con `do.call`.

ASMD_A <- do.call(cbind, before("A"))
ASMD_B <- do.call(cbind, before("B"))
ASMD_C <- do.call(cbind, before("C"))
ASMD_D <- do.call(cbind, before("D"))
ASMD_E <- do.call(cbind, before("E"))

## Hallamos el promedio y la desviación de ASMD para cada covariable:

mean_ASMD_A <- apply(ASMD_A, 1, mean) ; sd_ASMD_A <- apply(ASMD_A, 1, sd)
mean_ASMD_B <- apply(ASMD_B, 1, mean) ; sd_ASMD_B <- apply(ASMD_B, 1, sd)
mean_ASMD_C <- apply(ASMD_C, 1, mean) ; sd_ASMD_C <- apply(ASMD_C, 1, sd)
mean_ASMD_D <- apply(ASMD_D, 1, mean) ; sd_ASMD_D <- apply(ASMD_D, 1, sd)
mean_ASMD_E <- apply(ASMD_E, 1, mean) ; sd_ASMD_E <- apply(ASMD_E, 1, sd)

## De forma similar calculamos el ASAM para cada dataframe (media de ASMD's). 
## Después, el promedio y la desviación entre todos los dataframes.

## IMPORTANTE: Sólo incluimos en el ASAM las variables confusoras, porque es
## para las cuales nos interesa evaluar el balance.
## Las variables confusoras son distintas para cada escenario.
## --> A: w1-w6
## --> B: w1-w6 y w10
## --> C: w1-w10
## --> D, E: w1, w3-w7 y w9-w10

ASAM_A <- apply(ASMD_A[c(1:6), ], 2, mean) # Confusores: w1-6
mean_ASAM_A <- mean(ASAM_A) ; sd_ASAM_A <- sd(ASAM_A)

ASAM_B <- apply(ASMD_B[c(1:6,10), ], 2, mean) # Confusores: w1-w6 y w10
mean_ASAM_B <- mean(ASAM_B) ; sd_ASAM_B <- sd(ASAM_B)

ASAM_C <- apply(ASMD_C, 2, mean) # Confusores: w1-w10
mean_ASAM_C <- mean(ASAM_C) ; sd_ASAM_C <- sd(ASAM_C)

ASAM_D <- apply(ASMD_D[c(1, 3:7, 9, 10), ], 2, mean) # Confusores: w1, w3-w7 y w9-w10
mean_ASAM_D <- mean(ASAM_D) ; sd_ASAM_D <- sd(ASAM_D)

ASAM_E <- apply(ASMD_E[c(1, 3:7, 9, 10), ], 2, mean) # Confusores: w1, w3-w7 y w9-w10
mean_ASAM_E <- mean(ASAM_E) ; sd_ASAM_E <- sd(ASAM_E)

## Guardamos esta información en una lista:

inicial <- list("A" = list("asam_mean" = mean_ASAM_A,
                           "asam_sd" = sd_ASAM_A,
                           "asmd_mean" = mean_ASMD_A,
                           "asmd_sd" = sd_ASMD_A),
                "B" = list("asam_mean" = mean_ASAM_B,
                           "asam_sd" = sd_ASAM_B,
                           "asmd_mean" = mean_ASMD_B,
                           "asmd_sd" = sd_ASMD_B),
                "C" = list("asam_mean" = mean_ASAM_C,
                           "asam_sd" = sd_ASAM_C,
                           "asmd_mean" = mean_ASMD_C,
                           "asmd_sd" = sd_ASMD_C),
                "D"= list("asam_mean" = mean_ASAM_D,
                          "asam_sd" = sd_ASAM_D,
                          "asmd_mean" = mean_ASMD_D,
                          "asmd_sd" = sd_ASMD_D),
                "E" = list("asam_mean" = mean_ASAM_E,
                       "asam_sd" = sd_ASAM_E,
                       "asmd_mean" = mean_ASMD_E,
                       "asmd_sd" = sd_ASMD_E))

save(inicial, file = "results/inicial.rda")

#-------------------------------------------------------------------------------

# 2. Estimación de Propensity Scores con Regresión Logística

## IMPORTANTE: Ocultamos a los modelos de estimación las variables complejas, aunque
## sean confusoras o afecten al tratamiento. Esto es importante para cumplir con uno 
## de los objetivos principales de este estudio, que es comprobar si técnicas de ML 
## son capaces de mejorar la estimación de los PS gracias a su capacidad para detectar
## relaciones complejas/no lineales e interacciones, que a menudo se omiten en los
## modelos tradicionales (regresión logística). 

## + Aplicando Matching

RL_m_A <- RL_matching("A")
RL_m_B <- RL_matching("B")
RL_m_C <- RL_matching("C")
RL_m_D <- RL_matching("D")
RL_m_E <- RL_matching("E")

## + Aplicando IPW

RL_w_A <- RL_weight("A")
RL_w_B <- RL_weight("B")
RL_w_C <- RL_weight("C")
RL_w_D <- RL_weight("D")
RL_w_E <- RL_weight("E")


## Preparamos los resultados de cada escenario, calculando el promedio y desviación de 
## los ASMD's obtenidos en todos los dataframes, así como del ASAM. El procedimiento
## es similar al anterior (punto 1).

## + Resultados del Matching

### Creación de matrices
ASMD_A_RL_m <- do.call(cbind, RL_m_A)
ASMD_B_RL_m <- do.call(cbind, RL_m_B)
ASMD_C_RL_m <- do.call(cbind, RL_m_C)
ASMD_D_RL_m <- do.call(cbind, RL_m_D)
ASMD_E_RL_m <- do.call(cbind, RL_m_E)

### Promedio y desviación de ASMD para cada covariable:
mean_ASMD_A_RL_m <- apply(ASMD_A_RL_m, 1, mean) ; sd_ASMD_A_RL_m <- apply(ASMD_A_RL_m, 1, sd)
mean_ASMD_B_RL_m <- apply(ASMD_B_RL_m, 1, mean) ; sd_ASMD_B_RL_m <- apply(ASMD_B_RL_m, 1, sd)
mean_ASMD_C_RL_m <- apply(ASMD_C_RL_m, 1, mean) ; sd_ASMD_C_RL_m <- apply(ASMD_C_RL_m, 1, sd)
mean_ASMD_D_RL_m <- apply(ASMD_D_RL_m, 1, mean) ; sd_ASMD_D_RL_m <- apply(ASMD_D_RL_m, 1, sd)
mean_ASMD_E_RL_m <- apply(ASMD_E_RL_m, 1, mean) ; sd_ASMD_E_RL_m <- apply(ASMD_E_RL_m, 1, sd)

### Promedio y desviación del ASAM:
ASAM_A_RL_m <- apply(ASMD_A_RL_m[c(1:6), ], 2, mean) 
mean_ASAM_A_RL_m <- mean(ASAM_A_RL_m) ; sd_ASAM_A_RL_m <- sd(ASAM_A_RL_m)

ASAM_B_RL_m <- apply(ASMD_B_RL_m[c(1:6,10), ], 2, mean) 
mean_ASAM_B_RL_m <- mean(ASAM_B_RL_m) ; sd_ASAM_B_RL_m <- sd(ASAM_B_RL_m)

ASAM_C_RL_m <- apply(ASMD_C_RL_m, 2, mean) 
mean_ASAM_C_RL_m <- mean(ASAM_C_RL_m) ; sd_ASAM_C_RL_m <- sd(ASAM_C_RL_m)

ASAM_D_RL_m <- apply(ASMD_D_RL_m[c(1, 3:7, 9, 10), ], 2, mean) 
mean_ASAM_D_RL_m <- mean(ASAM_D_RL_m) ; sd_ASAM_D_RL_m <- sd(ASAM_D_RL_m)

ASAM_E_RL_m <- apply(ASMD_E_RL_m[c(1, 3:7, 9, 10), ], 2, mean) 
mean_ASAM_E_RL_m <- mean(ASAM_E_RL_m) ; sd_ASAM_E_RL_m <- sd(ASAM_E_RL_m)


## + Resultados del IPW

### Creación de matrices
ASMD_A_RL_w <- do.call(cbind, RL_w_A)
ASMD_B_RL_w <- do.call(cbind, RL_w_B)
ASMD_C_RL_w <- do.call(cbind, RL_w_C)
ASMD_D_RL_w <- do.call(cbind, RL_w_D)
ASMD_E_RL_w <- do.call(cbind, RL_w_E)

### Promedio y desviación de ASMD para cada covariable:
mean_ASMD_A_RL_w <- apply(ASMD_A_RL_w, 1, mean) ; sd_ASMD_A_RL_w <- apply(ASMD_A_RL_w, 1, sd)
mean_ASMD_B_RL_w <- apply(ASMD_B_RL_w, 1, mean) ; sd_ASMD_B_RL_w <- apply(ASMD_B_RL_w, 1, sd)
mean_ASMD_C_RL_w <- apply(ASMD_C_RL_w, 1, mean) ; sd_ASMD_C_RL_w <- apply(ASMD_C_RL_w, 1, sd)
mean_ASMD_D_RL_w <- apply(ASMD_D_RL_w, 1, mean) ; sd_ASMD_D_RL_w <- apply(ASMD_D_RL_w, 1, sd)
mean_ASMD_E_RL_w <- apply(ASMD_E_RL_w, 1, mean) ; sd_ASMD_E_RL_w <- apply(ASMD_E_RL_w, 1, sd)

### Promedio y desviación del ASAM:
ASAM_A_RL_w <- apply(ASMD_A_RL_w[c(1:6), ], 2, mean) 
mean_ASAM_A_RL_w <- mean(ASAM_A_RL_w) ; sd_ASAM_A_RL_w <- sd(ASAM_A_RL_w)

ASAM_B_RL_w <- apply(ASMD_B_RL_w[c(1:6,10), ], 2, mean) 
mean_ASAM_B_RL_w <- mean(ASAM_B_RL_w) ; sd_ASAM_B_RL_w <- sd(ASAM_B_RL_w)

ASAM_C_RL_w <- apply(ASMD_C_RL_w, 2, mean) 
mean_ASAM_C_RL_w <- mean(ASAM_C_RL_w) ; sd_ASAM_C_RL_w <- sd(ASAM_C_RL_w)

ASAM_D_RL_w <- apply(ASMD_D_RL_w[c(1, 3:7, 9, 10), ], 2, mean) 
mean_ASAM_D_RL_w <- mean(ASAM_D_RL_w) ; sd_ASAM_D_RL_w <- sd(ASAM_D_RL_w)

ASAM_E_RL_w <- apply(ASMD_E_RL_w[c(1, 3:7, 9, 10), ], 2, mean) 
mean_ASAM_E_RL_w <- mean(ASAM_E_RL_w) ; sd_ASAM_E_RL_w <- sd(ASAM_E_RL_w)

## Guardamos esta información en listas:

# Resultados del matching
RL_m <- list("A" = list("asam_mean" = mean_ASAM_A_RL_m,
                           "asam_sd" = sd_ASAM_A_RL_m,
                           "asmd_mean" = mean_ASMD_A_RL_m,
                           "asmd_sd" = sd_ASMD_A_RL_m),
                "B" = list("asam_mean" = mean_ASAM_B_RL_m,
                           "asam_sd" = sd_ASAM_B_RL_m,
                           "asmd_mean" = mean_ASMD_B_RL_m,
                           "asmd_sd" = sd_ASMD_B_RL_m),
                "C" = list("asam_mean" = mean_ASAM_C_RL_m,
                           "asam_sd" = sd_ASAM_C_RL_m,
                           "asmd_mean" = mean_ASMD_C_RL_m,
                           "asmd_sd" = sd_ASMD_C_RL_m),
                "D"= list("asam_mean" = mean_ASAM_D_RL_m,
                          "asam_sd" = sd_ASAM_D_RL_m,
                          "asmd_mean" = mean_ASMD_D_RL_m,
                          "asmd_sd" = sd_ASMD_D_RL_m),
                "E" = list("asam_mean" = mean_ASAM_E_RL_m,
                           "asam_sd" = sd_ASAM_E_RL_m,
                           "asmd_mean" = mean_ASMD_E_RL_m,
                           "asmd_sd" = sd_ASMD_E_RL_m))

# Resultados del IPW

RL_w <- list("A" = list("asam_mean" = mean_ASAM_A_RL_w,
                        "asam_sd" = sd_ASAM_A_RL_w,
                        "asmd_mean" = mean_ASMD_A_RL_w,
                        "asmd_sd" = sd_ASMD_A_RL_w),
             "B" = list("asam_mean" = mean_ASAM_B_RL_w,
                        "asam_sd" = sd_ASAM_B_RL_w,
                        "asmd_mean" = mean_ASMD_B_RL_w,
                        "asmd_sd" = sd_ASMD_B_RL_w),
             "C" = list("asam_mean" = mean_ASAM_C_RL_w,
                        "asam_sd" = sd_ASAM_C_RL_w,
                        "asmd_mean" = mean_ASMD_C_RL_w,
                        "asmd_sd" = sd_ASMD_C_RL_w),
             "D"= list("asam_mean" = mean_ASAM_D_RL_w,
                       "asam_sd" = sd_ASAM_D_RL_w,
                       "asmd_mean" = mean_ASMD_D_RL_w,
                       "asmd_sd" = sd_ASMD_D_RL_w),
             "E" = list("asam_mean" = mean_ASAM_E_RL_w,
                        "asam_sd" = sd_ASAM_E_RL_w,
                        "asmd_mean" = mean_ASMD_E_RL_w,
                        "asmd_sd" = sd_ASMD_E_RL_w))

save(RL_m, RL_w, file = "results/RL_results.rda")

#-------------------------------------------------------------------------------

# 2. Tuning RANDOM FOREST

### La función `hp_testing_m` realiza los siguientes pasos:
### 1) Estimación de PS según el método de ML especificado ("randomForest", "gbm", "rpart")
###    ajustando por un hiperparámetro especificado (hp). En el caso de RF, es 
###   `mtry`; para GBM, `shrinkage`, y para CART, `cp`.
### 2) Matching según los PS estimados.
### 3) Devuelve los ASMD's

## 2.1. Testar diferentes hiperparámetros

### Aplicamos dicha función con `lapply()` para testar varios hiperparámetros. El
### resultado es una lista en la que cada elemento guarda los resultados obtenidos 
### con un hiperparámetro distinto. Los resultados, son dataframes 10 x 100 con 
### los datos de ASMD para cada covariable (10 en total), en cada dataframe (100 en total).

### Para los Random Forest, testaremos `mtry` = 2, 3, 4, 5, 6

### Resultados para Matching

hp_test_A_m <- lapply(c(2:6), function(x){ hp_testing_m(x, "A", "randomForest")})
hp_test_B_m <- lapply(c(2:6), function(x){ hp_testing_m(x, "B", "randomForest")})
hp_test_C_m <- lapply(c(2:6), function(x){ hp_testing_m(x, "C", "randomForest")})
hp_test_D_m <- lapply(c(2:6), function(x){ hp_testing_m(x, "D", "randomForest")})
hp_test_E_m <- lapply(c(2:6), function(x){ hp_testing_m(x, "E", "randomForest")})

### Resultados para IPW

hp_test_A_w <- lapply(c(2:6), function(x){ hp_testing_w(x, "A", "randomForest")})
hp_test_B_w <- lapply(c(2:6), function(x){ hp_testing_w(x, "B", "randomForest")})
hp_test_C_w <- lapply(c(2:6), function(x){ hp_testing_w(x, "C", "randomForest")})
hp_test_D_w <- lapply(c(2:6), function(x){ hp_testing_w(x, "D", "randomForest")})
hp_test_E_w <- lapply(c(2:6), function(x){ hp_testing_w(x, "E", "randomForest")})


## 2.2. Selección del hiperparámetro óptimo

### La función `opt_hp` está diseñada para indicar qué hiperparámetro de los testados
### anteriormente da menor ASAM.

### Resultados para Matching

opt_hp(hp_test_A_m, "A") # mtry = 2 para escenario A / matching
opt_hp(hp_test_B_m, "B") # mtry = 2 para escenario B / matching
opt_hp(hp_test_C_m, "C") # mtry = 2 para escenario C / matching
opt_hp(hp_test_D_m, "D") # mtry = 2 para escenario D / matching
opt_hp(hp_test_E_m, "E") # mtry = 2 para escenario E / matching

### Resultados para IPW

opt_hp(hp_test_A_w, "A") # mtry = 3 para escenario A / IPW
opt_hp(hp_test_B_w, "B") # mtry = 3 para escenario B / IPW
opt_hp(hp_test_C_w, "C") # mtry = 2 para escenario C / IPW
opt_hp(hp_test_D_w, "D") # mtry = 2 para escenario D / IPW
opt_hp(hp_test_E_w, "E") # mtry = 2 para escenario E / IPW

#--------------------------------------------------------------------------------

# 3. Tuning GRADIENT BOOSTING MACHINE

## 3.1. Testar diferentes hiperparámetros

### Para los Gradient Boosting Machines, testaremos `shrinkage` = 0.1, 0.01 y 0.001.

### Resultados para Matching

gbm_hp_test_A_m <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_m(x, "A", "gbm")})
gbm_hp_test_B_m <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_m(x, "B", "gbm")})
gbm_hp_test_C_m <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_m(x, "C", "gbm")})
gbm_hp_test_D_m <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_m(x, "D", "gbm")})
gbm_hp_test_E_m <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_m(x, "E", "gbm")})

### Resultados para IPW

gbm_hp_test_A_w <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_w(x, "A", "gbm")})
gbm_hp_test_B_w <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_w(x, "B", "gbm")})
gbm_hp_test_C_w <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_w(x, "C", "gbm")})
gbm_hp_test_D_w <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_w(x, "D", "gbm")})
gbm_hp_test_E_w <- lapply(c(0.1, 0.01, 0.001), function(x){ hp_testing_w(x, "E", "gbm")})

## 3.2. Selección del hiperparámetro óptimo

### Resultados para Matching

opt_hp(gbm_hp_test_A_m, "A") # shrinkage = 0.1 para escenario A / matching
opt_hp(gbm_hp_test_B_m, "B") # shrinkage = 0.1 para escenario B / matching
opt_hp(gbm_hp_test_C_m, "C") # shrinkage = 0.001 para escenario C / matching
opt_hp(gbm_hp_test_D_m, "D") # shrinkage = 0.001 para escenario D / matching
opt_hp(gbm_hp_test_E_m, "E") # shrinkage = 0.001 para escenario E / matching

### Resultados para IPW

opt_hp(gbm_hp_test_A_w, "A") # shrinkage = 0.1 para escenario A / IPW
opt_hp(gbm_hp_test_B_w, "B") # shrinkage = 0.1 para escenario B / IPW
opt_hp(gbm_hp_test_C_w, "C") # shrinkage = 0.1 para escenario C / IPW
opt_hp(gbm_hp_test_D_w, "D") # shrinkage = 0.1 para escenario D / IPW
opt_hp(gbm_hp_test_E_w, "E") # shrinkage = 0.1 para escenario E / IPW


#-------------------------------------------------------------------------------

# 4. Tuning CLASSIFICATION & REGRESSION TREES

## 4.1. Testar diferentes hiperparámetros

### Para los Classification and Regression Trees, testaremos `cp` = 0.1, 0.01, 0.001

### Resultados para Matching

cart_hp_test_A_m <- lapply(c(0.01, 0.001), function(x){ hp_testing_m(x, "A", "rpart")})
cart_hp_test_B_m <- lapply(c(0.01, 0.001), function(x){ hp_testing_m(x, "B", "rpart")})
cart_hp_test_C_m <- lapply(c(0.01, 0.001), function(x){ hp_testing_m(x, "C", "rpart")})
cart_hp_test_D_m <- lapply(c(0.01, 0.001), function(x){ hp_testing_m(x, "D", "rpart")})
cart_hp_test_E_m <- lapply(c(0.01, 0.001), function(x){ hp_testing_m(x, "E", "rpart")})

### Resultados para IPW

cart_hp_test_A_w <- lapply(c(0.01, 0.001), function(x){ hp_testing_w(x, "A", "rpart")})
cart_hp_test_B_w <- lapply(c(0.01, 0.001), function(x){ hp_testing_w(x, "B", "rpart")})
cart_hp_test_C_w <- lapply(c(0.01, 0.001), function(x){ hp_testing_w(x, "C", "rpart")})
cart_hp_test_D_w <- lapply(c(0.01, 0.001), function(x){ hp_testing_w(x, "D", "rpart")})
cart_hp_test_E_w <- lapply(c(0.01, 0.001), function(x){ hp_testing_w(x, "E", "rpart")})

## 4.2. Selección del hiperparámetro óptimo

### Resultados para Matching

opt_hp(cart_hp_test_A_m, "A") # cp = 0.01 para escenario A / matching
opt_hp(cart_hp_test_B_m, "B") # cp = 0.001 para escenario B / matching
opt_hp(cart_hp_test_C_m, "C") # cp = 0.001 para escenario C / matching
opt_hp(cart_hp_test_D_m, "D") # cp = 0.001 para escenario D / matching
opt_hp(cart_hp_test_E_m, "E") # cp = 0.001 para escenario E / matching

### Resultados para IPW

opt_hp(cart_hp_test_A_w, "A") # cp = 0.001 para escenario A / IPW
opt_hp(cart_hp_test_B_w, "B") # cp = 0.001 para escenario B / IPW
opt_hp(cart_hp_test_C_w, "C") # cp = 0.001 para escenario C / IPW
opt_hp(cart_hp_test_E_w, "E") # cp = 0.001 para escenario D / IPW
opt_hp(cart_hp_test_D_w, "D") # cp = 0.001 para escenario E / IPW


## Guardamos todos los resultados:

save(rf_hp_test_A_m, rf_hp_test_B_m, rf_hp_test_C_m, rf_hp_test_D_m, rf_hp_test_E_m,
     rf_hp_test_A_w, rf_hp_test_B_w, rf_hp_test_C_w, rf_hp_test_D_w, rf_hp_test_E_w,
     gbm_hp_test_A_m, gbm_hp_test_B_m, gbm_hp_test_C_m, gbm_hp_test_D_m, gbm_hp_test_E_m,
     gbm_hp_test_A_w, gbm_hp_test_B_w, gbm_hp_test_C_w, gbm_hp_test_D_w, gbm_hp_test_E_w,
     cart_hp_test_A_m, cart_hp_test_B_m, cart_hp_test_C_m, cart_hp_test_D_m, cart_hp_test_E_m,
     cart_hp_test_A_m, cart_hp_test_B_m, cart_hp_test_C_m, cart_hp_test_D_m, cart_hp_test_E_m,
     file = "results/hp_optimization.rda")

#-------------------------------------------------------------------------------

# NOTA: Se pueden consultar los valores de ASAM promedio obtenidos para cada hp con el
# siguiente código; en este caso, para el escenario E utilizando IPW. Se imprimen los
# ASAM en el orden en que se han especificado los hp:

for (df in cart_hp_test_E_w){
  asam <- mean(apply(df, 2, mean))
  print(asam)
}

# Es decir, si hemos testado para mtry= 2,3,4,5 y 6; el primer valor que se imprime 
# corresponde a mtry=2.

#-------------------------------------------------------------------------------

# 5. Análisis de cada escenario:

## En primer lugar, extraemos los datos con los que trabajaremos. 

rf_A_m <- rf_hp_test_A_m[[1]] ; rf_A_w <- rf_hp_test_A_w[[2]]
rf_B_m <- rf_hp_test_B_m[[1]] ; rf_B_w <- rf_hp_test_B_w[[2]]
rf_C_m <- rf_hp_test_C_m[[1]] ; rf_C_w <- rf_hp_test_C_w[[1]]
rf_D_m <- rf_hp_test_D_m[[1]] ; rf_D_w <- rf_hp_test_D_w[[1]]
rf_E_m <- rf_hp_test_E_m[[1]] ; rf_E_w <- rf_hp_test_E_w[[1]]

gbm_A_m <- gbm_hp_test_A_m[[1]] ; gbm_A_w <- gbm_hp_test_A_w[[1]]
gbm_B_m <- gbm_hp_test_B_m[[1]] ; gbm_B_w <- gbm_hp_test_B_w[[1]]
gbm_C_m <- gbm_hp_test_C_m[[3]] ; gbm_C_w <- gbm_hp_test_C_w[[1]]
gbm_D_m <- gbm_hp_test_D_m[[3]] ; gbm_D_w <- gbm_hp_test_D_w[[1]]
gbm_E_m <- gbm_hp_test_E_m[[3]] ; gbm_E_w <- gbm_hp_test_E_w[[1]]

cart_A_m <- gbm_hp_test_A_m[[1]] ; cart_A_w <- cart_hp_test_A_w[[2]]
cart_B_m <- gbm_hp_test_B_m[[2]] ; cart_B_w <- cart_hp_test_B_w[[2]]
cart_C_m <- gbm_hp_test_C_m[[2]] ; cart_C_w <- cart_hp_test_C_w[[2]]
cart_D_m <- gbm_hp_test_D_m[[2]] ; cart_D_w <- cart_hp_test_D_w[[2]]
cart_E_m <- gbm_hp_test_E_m[[2]] ; cart_E_w <- cart_hp_test_E_w[[2]]


# Reunimos los resultados por escenario y por método de aplicación:

## Escenario A / Matching                     ## Escenario A / IPW

A_m <- list(rf_A_m, gbm_A_m, cart_A_m) ; A_w <- list(rf_A_w, gbm_A_w, cart_A_w)
names(A_m) <- c("rf", "gbm", "cart")   ; names(A_w) <- c("rf", "gbm", "cart")

## Escenario B / Matching                     ## Escenario B / IPW

B_m <- list(rf_B_m, gbm_B_m, cart_B_m) ; B_w <- list(rf_B_w, gbm_B_w, cart_B_w)
names(B_m) <- c("rf", "gbm", "cart")   ; names(B_w) <- c("rf", "gbm", "cart")

## Escenario C / Matching                     ## Escenario C / IPW

C_m <- list(rf_C_m, gbm_C_m, cart_C_m) ; C_w <- list(rf_C_w, gbm_C_w, cart_C_w)
names(C_m) <- c("rf", "gbm", "cart")   ; names(C_w) <- c("rf", "gbm", "cart")

## Escenario D / Matching                     ## Escenario D / IPW

D_m <- list(rf_D_m, gbm_D_m, cart_D_m) ; D_w <- list(rf_D_w, gbm_D_w, cart_D_w)
names(D_m) <- c("rf", "gbm", "cart")   ; names(D_w) <- c("rf", "gbm", "cart")

## Escenario E / Matching                     ## Escenario E / IPW

E_m <- list(rf_E_m, gbm_E_m, cart_E_m) ; E_w <- list(rf_E_w, gbm_E_w, cart_E_w)
names(E_m) <- c("rf", "gbm", "cart")   ; names(E_w) <- c("rf", "gbm", "cart")


## Aplicamos la función `balance_results` a cada lista:

## Escenario A / Matching                                  ## Escenario A / IPW
results_A_m <- lapply(A_m, function(x){balance(x, "A")}) ; results_A_w <- lapply(A_w, function(x){balance(x, "A")})

## Escenario B / Matching                                  ## Escenario B / IPW
results_B_m <- lapply(B_m, function(x){balance(x, "B")}) ; results_B_w <- lapply(B_w, function(x){balance(x, "B")})

## Escenario C / Matching                                  ## Escenario C / IPW
results_C_m <- lapply(C_m, function(x){balance(x, "C")}) ; results_C_w <- lapply(C_w, function(x){balance(x, "C")})

## Escenario D / Matching                                  ## Escenario D / IPW
results_D_m <- lapply(D_m, function(x){balance(x, "D")}) ; results_D_w <- lapply(D_w, function(x){balance(x, "D")})

## Escenario E / Matching                                  ## Escenario E / IPW
results_E_m <- lapply(E_m, function(x){balance(x, "E")}) ; results_E_w <- lapply(E_w, function(x){balance(x, "E")})


## Añadimos a estos resultados los de la regresión logística:

# Matching

results_A_m[["logr"]] <- list("asam_mean" = RL_m$A$asam_mean,
                              "asam_sd" = RL_m$A$asam_sd,
                              "asmd_mean" = RL_m$A$asmd_mean,
                              "asmd_sd" = RL_m$A$asmd_sd)

results_B_m[["logr"]] <- list("asam_mean" = RL_m$B$asam_mean,
                              "asam_sd" = RL_m$B$asam_sd,
                              "asmd_mean" = RL_m$B$asmd_mean,
                              "asmd_sd" = RL_m$B$asmd_sd)

results_C_m[["logr"]] <- list("asam_mean" = RL_m$C$asam_mean,
                              "asam_sd" = RL_m$C$asam_sd,
                              "asmd_mean" = RL_m$C$asmd_mean,
                              "asmd_sd" = RL_m$C$asmd_sd)

results_D_m[["logr"]] <- list("asam_mean" = RL_m$D$asam_mean,
                              "asam_sd" = RL_m$D$asam_sd,
                              "asmd_mean" = RL_m$D$asmd_mean,
                              "asmd_sd" = RL_m$D$asmd_sd)

results_E_m[["logr"]] <- list("asam_mean" = RL_m$E$asam_mean,
                              "asam_sd" = RL_m$E$asam_sd,
                              "asmd_mean" = RL_m$E$asmd_mean,
                              "asmd_sd" = RL_m$E$asmd_sd)

# IPW

results_A_w[["logr"]] <- list("asam_mean" = RL_w$A$asam_mean,
                              "asam_sd" = RL_w$A$asam_sd,
                              "asmd_mean" = RL_w$A$asmd_mean,
                              "asmd_sd" = RL_w$A$asmd_sd)

results_B_w[["logr"]] <- list("asam_mean" = RL_w$B$asam_mean,
                              "asam_sd" = RL_w$B$asam_sd,
                              "asmd_mean" = RL_w$B$asmd_mean,
                              "asmd_sd" = RL_w$B$asmd_sd)

results_C_w[["logr"]] <- list("asam_mean" = RL_w$C$asam_mean,
                              "asam_sd" = RL_w$C$asam_sd,
                              "asmd_mean" = RL_w$C$asmd_mean,
                              "asmd_sd" = RL_w$C$asmd_sd)

results_D_w[["logr"]] <- list("asam_mean" = RL_w$D$asam_mean,
                              "asam_sd" = RL_w$D$asam_sd,
                              "asmd_mean" = RL_w$D$asmd_mean,
                              "asmd_sd" = RL_w$D$asmd_sd)

results_E_w[["logr"]] <- list("asam_mean" = RL_w$E$asam_mean,
                              "asam_sd" = RL_w$E$asam_sd,
                              "asmd_mean" = RL_w$E$asmd_mean,
                              "asmd_sd" = RL_w$E$asmd_sd)


# Guardamos los resultados:

save(results_A_m, results_A_w, results_B_m, results_B_w, results_C_m, results_C_w,
     results_D_m, results_D_w, results_E_m, results_E_w, file = "results/all_results.rda")

#-------------------------------------------------------------------------------
# 6. Distribución de PS estimados 

## Trabajaremos con un único dataframe por escenario. Seleccionamos el primero.

dfA <- all_scenarios$A[[1]]
dfB <- all_scenarios$B[[1]]
dfC <- all_scenarios$C[[1]]

## Realizamos la estimación de los PS por cada técnica, para cada uno de estos 
## dataframes. 
## Para las técnicas de ML, utilizamos como hiperparámetros los seleccionados anteriormente:

# Regresión logística

eps_RL_A <- simple_RL(dfA) ; eps_RL_B <- simple_RL(dfB) ; eps_RL_C <- simple_RL(dfC)

# Random Forest

eps_RF_A <- simple_RF(dfA, 3) ; eps_RF_B <- simple_RF(dfB, 3) ;  eps_RF_C <- simple_RF(dfC, 3)

# GBM

eps_GBM_A <- simple_GBM(dfA, 0.001) ; eps_GBM_B <- simple_GBM(dfB, 0.001) ;  eps_GBM_C <- simple_GBM(dfC, 0.001)

# CART

eps_CART_A <- simple_CART(dfA, 0.1) ; eps_CART_B <- simple_CART(dfB, 0.1) ;  eps_CART_C <- simple_CART(dfC,0.1)

simple_estimations <- list("A" = list("logreg" = data.frame("Tr" = dfA$Tr, "tps" = dfA$tps, "eps" = eps_RL_A),
                                      "rf" = data.frame("Tr" = dfA$Tr,"tps" = dfA$tps, "eps" = eps_RF_A),
                                      "gbm" = data.frame("Tr" = dfA$Tr,"tps" = dfA$tps, "eps" = eps_GBM_A),
                                      "cart" = data.frame("Tr" = dfA$Tr,"tps" = dfA$tps, "eps" = eps_CART_A)),
                           
                           "B" = list("logreg" = data.frame("Tr" = dfB$Tr,"tps" = dfB$tps, "eps" = eps_RL_B),
                                      "rf" = data.frame("Tr" = dfB$Tr,"tps" = dfB$tps, "eps" = eps_RF_B),
                                      "gbm" = data.frame("Tr" = dfB$Tr,"tps" = dfB$tps, "eps" = eps_GBM_B),
                                      "cart" = data.frame("Tr" = dfB$Tr,"tps" = dfB$tps, "eps" = eps_CART_B)),
                           
                           "C" = list("logreg" = data.frame("Tr" = dfC$Tr,"tps" = dfC$tps, "eps" = eps_RL_C),
                                      "rf" = data.frame("Tr" = dfC$Tr,"tps" = dfC$tps, "eps" = eps_RF_C),
                                      "gbm" = data.frame("Tr" = dfC$Tr,"tps" = dfC$tps, "eps" = eps_GBM_C),
                                      "cart" = data.frame("Tr" = dfC$Tr,"tps" = dfC$tps, "eps" = eps_CART_C)))

save(simple_estimations, file = "results/simple_estimations.rda")

#-------------------------------------------------------------------------------

# 7. Estimación del ATE en cada dataframe de cada escenario

ATE_est_RL_A <- lapply(all_scenarios$A, function(x){ate_est(x, "logreg")})
ATE_est_RL_B <- lapply(all_scenarios$B, function(x){ate_est(x, "logreg")})
ATE_est_RL_B <- lapply(all_scenarios$C, function(x){ate_est(x, "logreg")})


ATE_est_rf_A <- lapply(all_scenarios$A, function(x){ate_est(x, "randomForest", 3)})
ATE_est_rf_B <- lapply(all_scenarios$B, function(x){ate_est(x, "randomForest", 3)})
ATE_est_rf_C <- lapply(all_scenarios$C, function(x){ate_est(x, "randomForest", 2)})

ATE_est_gbm_A <- lapply(all_scenarios$A, function(x){ate_est(x, "gbm", 0.1)})
ATE_est_gbm_B <- lapply(all_scenarios$B, function(x){ate_est(x, "gbm", 0.1)})
ATE_est_gbm_C <- lapply(all_scenarios$C, function(x){ate_est(x, "gbm", 0.1)})

ATE_est_cart_A <- lapply(all_scenarios$A, function(x){ate_est(x, "cart", 0.001)})
ATE_est_cart_B <- lapply(all_scenarios$B, function(x){ate_est(x, "cart", 0.001)})
ATE_est_cart_C <- lapply(all_scenarios$C, function(x){ate_est(x, "cart", 0.001)})

# Sesgo en la estimación del ATE

bias_RL_A <- abs(mean(unlist(ATE_est_RL_A)) - (-1))
bias_RL_B <-abs(mean(unlist(ATE_est_RL_B)) - (-1))
bias_RL_C <-abs(mean(unlist(ATE_est_RL_C)) - (-1))

bias_RF_A <-abs(mean(unlist(ATE_est_rf_A)) - (-1))
bias_RF_B <-abs(mean(unlist(ATE_est_rf_B)) - (-1))
bias_RF_C <-abs(mean(unlist(ATE_est_rf_C)) - (-1))

bias_GBM_A <- abs(mean(unlist(ATE_est_gbm_A)) - (-1))
bias_GBM_B <-abs(mean(unlist(ATE_est_gbm_B)) - (-1))
bias_GBM_C <-abs(mean(unlist(ATE_est_gbm_C)) - (-1))

bias_CART_A <-abs(mean(unlist(ATE_est_cart_A)) - (-1))
bias_CART_B <-abs(mean(unlist(ATE_est_cart_B)) - (-1))
bias_CART_C <-abs(mean(unlist(ATE_est_cart_C)) - (-1))

