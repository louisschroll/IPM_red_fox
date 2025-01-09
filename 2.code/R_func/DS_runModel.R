

DS_runModel <- function(data_DS,
                         dist_max = 0.6,
                         transect_len = 2,
                         size_study_area = 250,
                         nsites = 50,
                         nz = 200) {
  B <- dist_max
  years <- unique(data_DS$year)
  nyears <- length(years)
  
  # Data augmentation: add "pseudo-individuals"
  nind <- nrow(data_DS)
  y <- c(data_DS$N_obs, rep(0, nz))
  
  # Augment 'site' and 'year' with random values for pseudo-individuals
  site <- c(data_DS$site, sample(1:nsites, nz, replace = TRUE))
  year <- c(data_DS$year, sample(years, nz, replace = TRUE))
  
  # Replace NA in d with B + 1 for pseudo-individuals
  d <- c(ifelse(is.na(data_DS$d), B + 1, data_DS$d), runif(nz, 0, B))
  
  # Bundle and summarize data set
  nimble.constants <- list(
    nyears = nyears,
    nsites = nsites,
    B = B,
    transect_len = transect_len,
    size_study_area = size_study_area,
    nind = nind,
    nz = nz, 
    year = year
  )
  
  nimble.data <- list(y = y, d = d, site = site)
  
  # Nimble model for line transect HDS (NOT point transects!)
  nimble_code <- nimbleCode({
    # Prior distributions
    for (t in 1:nyears){
      beta0[t] ~ dunif(-2, 10)   # Intercept of lambda for each year
    }
    alpha0 ~ dunif(-10, 10)               # Intercept of log(sigma)
    log(sigma) <- alpha0
    
    for (t in 1:nyears) {
      # psi is a derived parameter under DA for stratified populations
      psi[t] <- sum(lambda[1:nsites, t]) / (nind[t] + nz[t])
      
      # 'Likelihood'
      for (i in 1:(nind[t] + nz[t])) {
        # i is index for individuals
        z[i, t] ~ dbern(psi[t])                               # Data augmentation variables
        d[i, t] ~ dunif(0, B + 1)                          # Allow range up to B + 1
        is_observed[i, t] <- (d[i, t] <= B)                  # Indicator for valid observations
        p[i, t] <- is_observed[i, t] * exp(-d[i, t] * d[i, t] / (2 * sigma * sigma))
        y[i, t] ~ dbern(z[i, t] * p[i, t])                       # Basic Bernoulli random variable
        site[i, t] ~ dcat(site.probs[1:nsites, t])
      }
      
      # Linear models for abundance and detection
      for (s in 1:nsites) {
        # s is index for sites
        N[s, t] ~ dpois(lambda[s, t])                  # Realized abundance at site s in year t
        log(lambda[s, t]) <- beta0[t]                 # Different lambda for each year
        site.probs[s, t] <- lambda[s, t] / sum(lambda[1:nsites, t])
      }
      
      # Derived parameters for each year
      Ntotal[t] <- sum(N[1:nsites, t])
      area[t] <- nsites * 2 * B * transect_len
      D[t] <- Ntotal[t] / area[t]
      N_gic[t] <- D[t] * size_study_area
    }
  })
  
  # Inits
  inits <- function() {
    alpha0 <- 0
    sigma <- exp(alpha0)
    beta0 <- rep(0, nyears)
    lambda <- matrix(exp(beta0), nrow = nsites, ncol = nyears)
    
    list(
      beta0 = beta0,
      alpha0 = alpha0,
      sigma = sigma,
      z = y,
      N = matrix(5, nrow = nsites, ncol = nyears)
    )
  }
  
  # Parameters to save
  params <- c("sigma", "beta0", "psi", "Ntotal", "N_gic", "D")
  
  # Build and compile the Nimble model
  nimble_model <- nimbleModel(
    code = nimble_code,
    data = nimble.data,
    constants = nimble.constants,
    inits = inits()
  )
  nimble_model$initializeInfo()
  print(nimble_model$calculate())
  
  compiled_model <- compileNimble(nimble_model)
  
  # Configure and build the MCMC
  mcmc_conf <- configureMCMC(nimble_model, monitors = params)
  mcmc <- buildMCMC(mcmc_conf)
  compiled_mcmc <- compileNimble(mcmc, project = nimble_model)
  
  # Run the MCMC
  mcmc_output <- runMCMC(
    compiled_mcmc,
    niter = 12000,
    nburnin = 2000,
    nchains = 3,
    thin = 1,
    setSeed = TRUE
  )
  
  return(mcmc_output)
}

# Analysis of DS data across all years -----------------------------
DS_out <- run_DS_model(
  data_DS = DS_data,
  nsites = n_sites,
  transect_len = transect_len,
  nz = 200
)

MCMCvis::MCMCtrace(DS_out)

N_estimate_DS <- map(DS_out, as_tibble) %>% bind_rows() %>% select(starts_with("N_gic["))
# N_estimates_upper <- map(DS_out, as_tibble) %>% bind_rows() %>% pull(N_gic_all) %>% quantile(probs = 0.025)
# N_estimates_lower <- map(DS_out, as_tibble) %>% bind_rows() %>% pull(N_gic_all) %>% quantile(probs = 0.975)

N_estimate_tibble <- tibble(
  N_DS = N_estimate_DS,
  q97.5 = N_estimates_upper,
  q2.5 = N_estimates_lower
)

N_estimate_tibble %>%
  ggplot(aes(x = 1, y = N_DS)) +
  geom_point() +
  geom_errorbar(aes(ymin = q2.5, ymax = q97.5), width = 0.1) +
  theme_minimal() +
  labs(x = "Year", y = "Estimated Population Size")
