

source("packages.R")
source("functions.R")

# `all_scenarios` is a dataframe that stores all the dataframes used
# for the analysis. Each column represents a different scenario (A–E).

all_scenarios <- as.data.frame(sapply(c("A","B","C","D", "E"), 
                                      function(x){rep_simulation(100,1000,x)}))

# Data saved in the `data` folder.
save(all_scenarios, file = "data/all_scenarios.rda")


# Ejemplo 1: Se accede a todos los dataframes del escenario "A" con:
# all_scenarios$A

# Ejemplo 2: Se accede al primer dataframe del escenario "A" con:
# all_scenarios$A[[1]]

