
# Write IPM Code
IPM.Code <- nimble::nimbleCode({
  ####################
  # POPULATION MODEL #
  ####################
  
  # Initial population size t = 1
  for (a in 1:n_age_class){
    N_dec[a] ~ dunif(1, 500)
    N[a, 1] <- round(N_dec[a])
  }
  
  # Population dynamic t>1
  for (t in 1:(n_years - 1)) {
    # Age class 0 (index = 1): local reproduction
    N[1, t + 1] ~ dpois(surv_rate * sum(fec_rate * N[2:n_age_class, t]))
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors
    for (a in 1:(n_age_class - 2)) {
      N[a + 1, t + 1] ~ dbinom(N[a, t], surv_rate)
    }
    
    # Age class 4+ (index = n_age_class = 5): age class 4 and 5+ survivors
    N[n_age_class, t + 1] ~ dbinom(N[n_age_class - 1, t] + N[n_age_class, t], surv_rate)
  }
  
  # Derived quantities
  for (t in 1:n_years){
    N_tot[t] <- sum(N[1:n_age_class, t])
  }
  
  ###########
  # PRIORS  #
  ###########
  fec_rate <- 1.6 # dunif(0, 20) # Recruitment
  surv_rate <- 0.56 # ~ dunif(0, 1)      # Survival
  
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
    succprob[t] <- kappa / (kappa + N_tot[t])
    DS_estimate[t] ~ dnegbin(prob = succprob[t], size = kappa)
    # DS_estimate[t] ~ dpois(N_tot[t])
  }
  
  # Priors
  kappa ~ dunif(min = 0.01, max = 10)
  
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
      C[a, t] ~ dpois(harvest_rate * N[a, t]) #dbin(harvest_rate, N[a, t])
    }
  }
})


# IPM_prepareInputData
## Reformat data into lists for Nimble
input_data <- list(
  nim.data = list(
    DS_estimate = round(N_estimates),
    C = harvest_data
  ),
  
  nim.constants = list(
    n_years = n_years,
    n_age_class = nrow(harvest_data)
  ))

# IPM_simulateInits
initVals <- list(
  N = N,
  N_tot = N_estimates,
  N_dec = N[1:n_age_class, 1],
  harvest_rate = 0.2,
  kappa = 3,
  succprob = 3 / (3 + round(N_estimates))
)

# IPM_setupModel
## Set parameters to monitor
params <- c("N_tot", "N")

niter = 10000
nthin = 1
nburn = 1000
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

N_estimate <- IPM_tibble %>% 
  select(starts_with("N_tot")) %>% 
  janitor::clean_names() %>%
  pivot_longer(everything(), names_to = "id") %>% 
  mutate(year = as.numeric(str_extract(id, "\\d+")),
         id = str_extract(id, "N_tot")) %>% 
  group_by(year, id) %>% 
  summarise(mean = mean(value),
            median = median(value),
            q2.5 = quantile(value, probs = 0.025),
            q97.5 = quantile(value, probs = 0.975))


N_estimate %>% 
  bind_rows(tibble(mean = c(colSums(N)),
                   year = 1:n_years,
                   id = "N_real")) %>% 
  ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
  geom_line() +
  geom_ribbon(aes(ymin = c(q2.5), 
                  ymax = c(q97.5)), 
              linetype = 2, alpha = 0.1) +
  theme_minimal() +
  coord_cartesian(ylim = c(0, 300))

