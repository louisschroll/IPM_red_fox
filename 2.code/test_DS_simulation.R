#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/test_DS_simulation.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-28
#'
#' Script description:
#' Test the accuracy and good functioning of the functions sim_DS_data use to 
#' simulate distance sampling data and run_DS_model that run a distance-sampling 
#' model in jags.
#'
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(jagsUI)       

# Load R functions
path_to_Rfunc <- "2.code/R_func"
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

## Test simulation and model ----
# Simulate distance-sampling data
DS_data <- sim_DS_data(nsites = 50)                  # Line transect (default)

# Run DS model and summarize posterior output
out1 <- run_DS_model(DS_data = DS_data$data)
print(out1, 2)
sum(DS_data$N.true)

# Check convergence
MCMCvis::MCMCtrace(out1, params = c("sigma", "beta0", "N_gic"), pdf = F)

## Test the model with repeated simulations ----
nsimu <- 100
df_simu <- tibble()

for (i in 1:nsimu) {
  DS_data <- sim_DS_data(nsites = 50)
  out1 <- run_DS_model(DS_data$data)
  df_simu <- df_simu %>% bind_rows(
    tibble(num_simu = i, 
           N_gic = out1$sims.list$N_gic, 
           alpha0 = out1$sims.list$alpha0, 
           beta0 = out1$sims.list$beta0, 
           psi = out1$sims.list$psi)
  )
  
}

df_simu_summary <- df_simu %>% select(num_simu, N_gic) %>% 
  group_by(num_simu) %>% 
  summarise(mean = mean(N_gic),
            med = median(N_gic))


boxplot(df_simu_summary$mean)
abline(h = 200, col = "red")


## Look the effect of data augmentation, i.e. number of pseudo-individuals (nz)
# Number of pseudo-individuals to test
nz_to_test <- c(50, 100, 200, 300, 500)

# Number of simulation for each nz
nsimu <- 1000

# list of nz values
nz_values <- rep(c(50, 100, 200, 300, 500), each = nsimu)

# Initialize df_simu
df_simu <- tibble()

for (i in 1:(nsimu * length(nz_to_test))) {
  # Simulate the data and run the model
  DS_data <- sim_DS_data(nsites = 50)
  out1 <- run_DS_model(DS_data$data, nz = nz_values[i])
  # Add new estimates (mean) to df_simu
  df_simu <- df_simu %>% bind_rows(
    tibble(num_simu = i, 
           nz = nz_values[i],
           N_gic = out1$mean$N_gic, 
           sigma = out1$mean$sigma, 
           beta0 = out1$mean$beta0, 
           psi = out1$mean$psi)
  )
}

plot1 <- df_simu %>% 
  ggplot(aes(x = nz, y = N_gic, color = as.factor(nz))) +
  geom_violin() +
  geom_boxplot() +
  theme_minimal() +
  geom_hline(yintercept = 200) +
  labs(title = "Population size estimates with multiple simulations",
       y = "Population size estimate",
       x = "Number of pseudo-individuals (DA)") +
  ylim(y = c(0, 400)) +
  theme(legend.position = "none")

plot2 <- df_simu %>% 
  ggplot(aes(x = nz, y = sigma, color = as.factor(nz))) +
  geom_violin() +
  geom_boxplot() +
  theme_minimal() +
  geom_hline(yintercept = 0.15) +
  labs(title = "Detection parameter estimates with multiple simulations",
       y = expression("Detection parameter" ~ sigma),
       x = "Number of pseudo-individuals (DA)") +
  #ylim(y = c(0, 400)) +
  theme(legend.position = "none")

library(patchwork)
plot1 + plot2
