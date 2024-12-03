#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/script_simulation_DS.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-22
#'
#' Script description:
#'
#'
#' -----------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(jagsUI)

library(nimble)
library(nimbleDistance)

# Load functions
path_to_Rfunc <- "2.code/R_func"
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

## Simulation of the data ------------------------------------------------------
# Set simulation parameters
n_years <- 10
n_age_class <- 5
N0 <- 200

n_sites <- 50
transect_len <- 6
size_hunting_area <- 250
dist_max <- 0.6
mean.sigma <- 0.15
# Simulate population dynamic
N <- sim_pop_dyn(n_years = n_years,
                 n_age_class = n_age_class,
                 N0 = N0)

N %>%
  colSums() %>% 
  as_tibble() %>% 
  mutate(year = 1:nrow(.)) %>% 
  ggplot(aes(x = year, y = value)) +
  geom_line() +
  geom_point() +
  theme_minimal()

t(N) %>% as_tibble() %>% 
  mutate(year = 1:nrow(.),
         Ntot = rowSums(.)) %>% 
  rename_with(~ paste0("age_", seq_along(.)), starts_with("V")) %>% 
  pivot_longer(-year,
               names_to = "age_class") %>% 
  ggplot(aes(x = year, y = value, group = as.factor(age_class), colour = as.factor(age_class))) +
  geom_line() +
  theme_minimal()

# Simulate distance-sampling data based on matrix N
DS_data <- sim_DS_data(Ntot = colSums(N), 
                       nsites = n_sites, 
                       mean.sigma = mean.sigma,
                       dist_max = dist_max,
                       transect_len = transect_len,
                       size_study_area = size_study_area)

# Simulate age-at-harvest matrix
harvest_rate <- 0.2
harvest_data <- sim_age_data(N, harvest_rate = harvest_rate)

## Analysis of DS data independently for each year -----------------------------
N_estimates <- N_estimates_upper <- N_estimates_lower <- N_estimates2 <- c()

for (i in 1:n_years){
  out1 <- run_DS_model(DS_data %>% filter(year == i) %>% select(-year), 
                       nsites = n_sites, 
                       transect_len = transect_len,
                       nz = 600)
  N_estimates <- c(N_estimates, out1$mean$N_gic)
  N_estimates2 <- c(N_estimates2, out1$mean$N_gic2)
  N_estimates_upper <- c(N_estimates_upper, out1$q2.5$N_gic)
  N_estimates_lower <- c(N_estimates_lower, out1$q97.5$N_gic)
}

tibble(year = 1:ncol(N),
       N_real = colSums(N),
       N_estimates = N_estimates,
       N_estimates_upper = N_estimates_upper,
       N_estimates_lower = N_estimates_lower) %>% 
  pivot_longer(-c(year, N_estimates_upper, N_estimates_lower),
               names_to = "Pop_size") %>%
  mutate(N_estimates_upper = ifelse(Pop_size == "N_estimates", N_estimates_upper, NA),
         N_estimates_lower = ifelse(Pop_size == "N_estimates", N_estimates_lower, NA)) %>% 
  ggplot(aes(x = year, y = value, group = as.factor(Pop_size), colour = as.factor(Pop_size))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(N_estimates_lower), 
                  ymax = c(N_estimates_upper)), 
              linetype=2, alpha=0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 400))


# ------

simData_list <- list(DS_data1 = DS_data, DS_data2 = DS_data %>% filter(site < 41))

ISDM_out <- runIntDSModel(DS_data_list = simData_list, 
                          harvest_data = C, 
                          W = dist_max, 
                          size_hunting_area = size_study_area)





