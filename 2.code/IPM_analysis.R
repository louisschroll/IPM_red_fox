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
transect_len <- 6
size_hunting_area <- 250
dist_max <- 0.6
mean.sigma <- 0.15
# Simulate population dynamic
pop_dyn_list <- sim_pop_dyn(n_years = n_years,
                 n_age_class = n_age_class,
                 N0 = N0)

N <- pop_dyn_list$N

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

# Simulate reproduction data
reprod_data <- sim_reprod_data(pop_dyn_list)


# Analysis of DS data across all years -----------------------------
DS_out <- DS_runModel(
  data_DS = DS_data,
  nsites = n_sites,
  transect_len = transect_len,
  nz = 1000
)

MCMCvis::MCMCtrace(DS_out)

N_estimate_tibble <- map(DS_out, as_tibble) %>% 
  bind_rows() %>% 
  select(starts_with("N_gic[")) %>% 
  janitor::clean_names() %>% 
  pivot_longer(everything(), names_to = "year", values_to = "N") %>% 
  mutate(year = as.numeric(str_remove(year, "n_gic_"))) %>% 
  group_by(year) %>% 
  summarise(mean = mean(N),
            q2.5 = quantile(N, probs = 0.025),
            q97.5 = quantile(N, probs = 0.975)) %>% 
  bind_rows(tibble(year = 1:n_years,
                   mean = colSums(N)), .id = "id")

N_estimate_tibble %>% 
  ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(q2.5),
                  ymax = c(q97.5)),
              linetype=2, alpha=0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 500)) +
  labs(x = "Year", y = "N")

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
  
  
  # First year (reproduction not modeled separately)
  B[1:n_age_class, 1] <- 0
  L[1:n_age_class, 1] <- 0
  R[1:n_age_class, 1] <- 0
  
  # Population dynamic t>1
  for (t in 1:(n_years - 1)) {
    # Age class 0 (index = 1): local reproduction
    for (a in 1:n_age_class) {
      # Breeding Population Size: Number of females that reproduce
      B[a, t + 1] ~ dbin(pregnancy_rate, round(N[a, t] / 2))
      
      # Litter Size (in utero): Number of pups produced by females of age class a
      L[a, t + 1] ~ dpois(B[a, t + 1] * rho * 0.5)
      
      # Number Recruits: Number of pups surviving to emerge from the den
      R[a, t + 1] ~ dbin(surv_rate, L[a, t + 1])
    }
    N[1, t + 1] <- sum(R[1:n_age_class, t + 1])
    
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
  surv_rate ~ dunif(0, 1)      # Survival
  pregnancy_rate ~ dunif(0, 1)
  rho ~ dunif(0, 20)
  
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
  
  #########################
  # PLACENTAL SCAR MODULE #
  #########################
  
  ### Parameters:
  # rho = expected number of placental scars (fetuses)
  # pregnancy_rate = pregnancy rate
  
  ## Data:
  # P1 = individual placental scar counts
  # P1_age = individual ages associated with P1
  # P1_year = year associated with P1
  
  # P2 = individual presence/absence of placental scars
  # P2_age = individual ages associated with P2
  # P2_year = year associated with P2
  
  
  ### Likelihood (litter size)

  for(x in 1:n_obs_litter){
    obs_litter_size[x] ~ dpois(rho)
  }

  ### Likelihood (pregnancy rate)

  for(x in 1:n_obs_female){
    isBreeding[x] ~ dbern(pregnancy_rate)
  }
})


# IPM_prepareInputData
## Reformat data into lists for Nimble
input_data <- list(
  nim.data = list(
    DS_estimate = round(N_estimate_tibble$mean),
    C = harvest_data,
    obs_litter_size = reprod_data$nbFoetus,
    isBreeding = reprod_data$isBreeding
  ),
  
  nim.constants = list(
    n_years = n_years,
    n_age_class = nrow(harvest_data),
    n_obs_litter = length(reprod_data$nbFoetus),
    n_obs_female = length(reprod_data$isBreeding)
  ))

# IPM_simulateInits
initVals <- list(
  N = N,
  N_tot = colSums(N),
  N_dec = N[1:n_age_class, 1],
  B = pop_dyn_list$B,
  L = pop_dyn_list$L,
  R = pop_dyn_list$NbRecruits,
  harvest_rate = 0.2,
  # kappa = 3,
  # succprob = 3 / (3 + round(N_estimate_DS)),
  surv_rate = runif(1, min = 0, max = 1),
  pregnancy_rate = runif(1, min = 0, max = 1),
  rho = rpois(1, 4)
)

# IPM_setupModel
## Set parameters to monitor
params <- c("N_tot", "N", 
            # "kappa", 
            "surv_rate",
            "pregnancy_rate",
            "rho")

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
  bind_rows(N_estimate_tibble %>% filter(id == 1)) %>% 
  ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(q2.5), 
                  ymax = c(q97.5)), 
              linetype = 2, alpha = 0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300))

