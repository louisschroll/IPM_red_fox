cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(jagsUI)
library(popdemo)
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
transect_len <- 2
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
  theme_minimal() +
  coord_cartesian(ylim = c(0, 500))

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
                       size_study_area = size_hunting_area)

# Simulate age-at-harvest matrix
harvest_rate <- 0.2
harvest_data <- sim_age_data(N, harvest_rate = harvest_rate)

# Analysis of DS data independently for each year -----------------------------
N_estimate_DS <- N_estimates_upper <- N_estimates_lower <- N_estimates2 <- c()

for (i in 1:n_years){
  out1 <- run_DS_model(DS_data %>% filter(year == i) %>% select(-year),
                       nsites = n_sites,
                       transect_len = transect_len,
                       nz = 600)
  N_estimate_DS <- c(N_estimate_DS, out1$mean$N_gic)
  N_estimates2 <- c(N_estimates2, out1$mean$N_gic2)
  N_estimates_upper <- c(N_estimates_upper, out1$q2.5$N_gic)
  N_estimates_lower <- c(N_estimates_lower, out1$q97.5$N_gic)
}

N_estimate_tibble <- tibble(year = 1:ncol(N),
       N_real = colSums(N),
       N_DS = N_estimate_DS,
       q97.5 = N_estimates_upper,
       q2.5 = N_estimates_lower) %>%
  pivot_longer(-c(year, q2.5, q97.5),
               names_to = "id", values_to = "mean") %>%
  mutate(q97.5 = ifelse(id == "N_DS", q97.5, NA),
         q2.5 = ifelse(id == "N_DS", q2.5, NA)) 

N_estimate_tibble %>%
  ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(q2.5),
                  ymax = c(q97.5)),
              linetype=2, alpha=0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 400))



# --------------
# Write IPM Code
IPM.Code <- nimble::nimbleCode({
  ####################
  # POPULATION MODEL #
  ####################
  
  # Initial population size t = 1
  for (a in 1:n_age_class){
    N_dec[a] ~ dunif(1, 500)
    N[a, 1] <- round(N_dec[a])
    # N[a, 1] ~ dcat(rep(1/500, 500))
  }
  
  # Population dynamic t>1
  for (t in 1:(n_years - 1)) {
    # Age class 0 (index = 1): local reproduction
    N[1, t + 1] ~ dpois(surv_rate * sum(fec_rate * N[2:n_age_class, t]))
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors
    for (a in 1:(n_age_class - 2)) {
      N[a + 1, t + 1] ~ dbinom(surv_rate, N[a, t])
    }
    
    # Age class 4+ (index = n_age_class = 5): age class 4 and 5+ survivors
    N[n_age_class, t + 1] ~ dbinom(surv_rate, N[n_age_class - 1, t] + N[n_age_class, t])
  }
  
  # Derived quantities
  for (t in 1:n_years){
    N_tot[t] <- sum(N[1:n_age_class, t])
  }
  
  ###########
  # PRIORS  #
  ###########
  fec_rate <- 1.5 # dunif(0, 20) # Recruitment
  surv_rate ~ dunif(0, 1)      # Survival
  
  ################
  # COUNT MODULE #
  ################
  ## Parameters:
  # N_tot = population size at a given time
  # h = age- and time-dependent probability of dying from hunting
  
  ## Data:
  # DS_estimate = estimation of population size using distance-sampling
  
  ## Likelihood
  for (t in 1:n_years){
    # succprob[t] <- kappa / (kappa + N_tot[t])
    # DS_estimate[t] ~ dnegbin(prob = succprob[t], size = kappa)
    DS_estimate[t] ~ dpois(N_tot[t])
  }
  
  # Priors
  # kappa ~ dunif(min = 0.01, max = 100)
  
  #########################
  # AGE-AT-HARVEST MODULE #
  #########################
  ## Parameters:
  # N = number of individuals in a given age class at a given time
  # h = age- and time-dependent probability of dying from hunting
  
  ## Data:
  # C = age-at-harvest matrix
  
  ## Priors
  harvest_rate ~ dunif(0, 1)
  
  ## Likelihood
  for(t in 1:n_years){
    for(a in 1:n_age_class){
      C[a, t] ~ dbin(harvest_rate, N[a, t])
    }
  }
})


# IPM_prepareInputData
## Reformat data into lists for Nimble
input_data <- list(
  nim.data = list(
    DS_estimate = round(N_estimate_DS),
    C = harvest_data
  ),
  
  nim.constants = list(
    n_years = n_years,
    n_age_class = nrow(harvest_data)
  ))

# IPM_simulateInits
initVals <- list(
  N = N,
  N_tot = colSums(N),
  N_dec = N[1:n_age_class, 1],
  harvest_rate = 0.2,
  # kappa = 3,
  succprob = 3 / (3 + round(N_estimate_DS)),
  surv_rate = runif(1, min = 0, max = 1)
)

# IPM_setupModel
## Set parameters to monitor
params <- c("N_tot", "N", 
            # "kappa", 
            "surv_rate")

niter = 40000
nthin = 1
nburn = 5000
nchains = 3

model_setup <- list(
  modelParams = params,
  initVals = initVals,
  mcmcParams = list(
    niter = niter,
    nthin = nthin,
    nburn = nburn,
    nchains = nchains
  )
)



model <- nimbleModel(IPM.Code, 
                     constants = input_data$nim.constants, 
                     data = input_data$nim.data,
                     inits = initVals)
model$initializeInfo()
simulate(model)
model$calculate("logProb_N")


# IPM_runModel
IPM_out <- nimbleMCMC(code = IPM.Code,
                       data = input_data$nim.data, 
                       constants = input_data$nim.constants,
                       inits = model_setup$initVals, 
                       monitors = model_setup$modelParams,
                       nchains = model_setup$mcmcParams$nchains, 
                       niter = model_setup$mcmcParams$niter, 
                       nburnin = model_setup$mcmcParams$nburn, 
                       thin = model_setup$mcmcParams$nthin, 
                       samplesAsCodaMCMC = TRUE)

MCMCvis::MCMCtrace(IPM_out)

IPM_tibble <- map(IPM_out, . %>% 
                     as_tibble()) %>% 
  bind_rows() 

N_estimate_IPM <- IPM_tibble %>% 
  select(starts_with("N_tot")) %>% 
  janitor::clean_names() %>%
  pivot_longer(everything(), names_to = "id") %>% 
  mutate(year = as.numeric(str_extract(id, "\\d+")),
         id = "N_IPM") %>% 
  group_by(year, id) %>% 
  summarise(mean = mean(value),
            median = median(value),
            q2.5 = quantile(value, probs = 0.025),
            q97.5 = quantile(value, probs = 0.975))


N_estimate_IPM %>% 
  bind_rows(tibble(mean = c(colSums(N)),
                   year = 1:n_years,
                   id = "N_real")) %>% 
  bind_rows(N_estimate_tibble %>% filter(id == "N_DS")) %>% 
  ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(q2.5), 
                  ymax = c(q97.5)), 
              linetype = 2, alpha = 0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 500))

