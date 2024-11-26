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

dyn_output <- run_dynamic_DS_model(DS_data)


t(N) %>% as_tibble() %>% 
  mutate(year = 1:nrow(.),
         Ntot = rowSums(.),
         N_estimates = N_estimates,
         N_estimates_upper = N_estimates_upper,
         N_estimates_lower = N_estimates_lower,
         N_est2 = N_estimates2) %>% 
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



# Data bundle
jags.data <- list(
  Nhat = out4$mean$totalN,
  var.Nhat = out4$sd$totalN ^ 2,
  T = length(out4$mean$totalN)
)

# Write JAGS model file
cat(
  file = "model3.txt",
  "model {
    # Priors and linear models
    mu.lam ~ dunif(0, 10) # Prior for mean growth rate
    sig.lam ~ dunif(0, 10) # Prior for sd of growth rate
    sig2.lam <- pow(sig.lam, 2)
    tau.lam <- pow(sig.lam, -2)
    sig.ystar ~ dunif(0, 10000) # Prior for sd of observation process
    sig2.ystar <- pow(sig.ystar, 2)
    tau.ystar <- pow(sig.ystar, -2)
    
    # Likelihood
    # Model for the initial population size: uniform priors
    N[1] ~ dunif(0, 500)
    # Process model over time: our model of population dynamics
    for (t in 1:(T - 1)) {
      lambda[t] ~ dnorm(mu.lam, tau.lam)
      N[t + 1] <- N[t] * lambda[t]
    }
    
    # Observation process for the abundance estimates with their SEs
    for (t in 1:T) {
      Nhat[t] ~ dnorm(ystar[t], tau.se[t])
      tau.se[t] <- 1 / var.Nhat[t] # Assumed known and given by variance of Nhat
      ystar[t] ~ dnorm(N[t], tau.ystar)
    }
  }
")

# Initial values
inits <- function() {
  list(sig.lam = runif(1, 0, 1))
}
# Parameters monitored
parameters <- c("lambda",
                "mu.lam",
                "sig2.ystar",
                "sig2.lam",
                "sig.ystar",
                "sig.lam",
                "N",
                "ystar")
# MCMC settings
ni <- 50000
nb <- 25000
nc <- 3
nt <- 25
na <- 10000
# Call JAGS from R (ART <1 min), check convergence and summarize posteriors
out5 <- jags(
  jags.data,
  inits,
  parameters,
  "model3.txt",
  n.iter = ni,
  n.burnin = nb,
  n.chains = nc,
  n.thin = nt,
  n.adapt = na,
  parallel = TRUE
)
traceplot(out5) # Not shown
print(out5, 3) # Not shown

