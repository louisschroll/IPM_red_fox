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
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

library(tidyverse)
library(jagsUI)

path_to_Rfunc <- "2.code/R_func"
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

Tmax <- 10
n_years <- Tmax + 1
N <- sim_pop_dyn(Tmax = Tmax)

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
  pivot_longer(-year,
               names_to = "age_class") %>% 
  ggplot(aes(x = year, y = value, group = as.factor(age_class), colour = as.factor(age_class))) +
  geom_line() +
  theme_minimal()

DS_data <- sim_DS_data(Ntot = colSums(N), nsites = 50, transect_len = 6)
N_estimates <- N_estimates_upper <- N_estimates_lower <- N_estimates2 <- c()

for (i in 1:(Tmax + 1)){
  out1 <- run_DS_model(DS_data %>% filter(year == i) %>% select(-year), 
                       nsites = 50, 
                       transect_len = 6,
                       nz = 600)
  N_estimates <- c(N_estimates, out1$mean$N_gic)
  N_estimates2 <- c(N_estimates2, out1$mean$N_gic2)
  N_estimates_upper <- c(N_estimates_upper, out1$q2.5$N_gic)
  N_estimates_lower <- c(N_estimates_lower, out1$q97.5$N_gic)
}

t(N) %>% as_tibble() %>% 
  mutate(year = 1:nrow(.),
         Ntot = rowSums(.),
         N_estimates = N_estimates,
         N_estimates_upper = N_estimates_upper,
         N_estimates_lower = N_estimates_lower) %>% 
  pivot_longer(-c(year, N_estimates_upper, N_estimates_lower),
               names_to = "age_class") %>%
  mutate(N_estimates_upper = ifelse(age_class == "N_estimates", N_estimates_upper, NA),
         N_estimates_lower = ifelse(age_class == "N_estimates", N_estimates_lower, NA)) %>% 
  ggplot(aes(x = year, y = value, group = as.factor(age_class), colour = as.factor(age_class))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(N_estimates_lower), 
                  ymax = c(N_estimates_upper)), 
                  linetype=2, alpha=0.1) +
  theme_minimal()


# ------

simData_list <- list(DS_data1 = DS_data, DS_data2 = DS_data %>% filter(site<41))

dyn_output <- run_dynamic_DS_model(DS_data)


t(N) %>% as_tibble() %>% 
  mutate(year = 1:nrow(.),
         Ntot = rowSums(.),
         N_estimates = N_estimates,
         N_estimates_upper = N_estimates_upper,
         N_estimates_lower = N_estimates_lower,
         N_est2 = dyn_output$mean$Ntot) %>% 
  pivot_longer(-c(year, N_estimates_upper, N_estimates_lower),
               names_to = "age_class") %>%
  mutate(N_estimates_upper = ifelse(age_class == "N_estimates", N_estimates_upper, NA),
         N_estimates_lower = ifelse(age_class == "N_estimates", N_estimates_lower, NA)) %>% 
  ggplot(aes(x = year, y = value, group = as.factor(age_class), colour = as.factor(age_class))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(N_estimates_lower), 
                  ymax = c(N_estimates_upper)), 
              linetype=2, alpha=0.1) +
  theme_minimal()




