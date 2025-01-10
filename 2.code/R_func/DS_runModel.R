DS_runModel <- function(data_DS,
                        dist_max = 0.6,
                        transect_len = 2,
                        size_study_area = 250,
                        nsites = 50,
                        nz = 300) {
  B <- dist_max
  years <- unique(data_DS$year)
  n_years <- length(years)
  
  # Data augmentation: Create matrices for variables
  nind_per_year <- sapply(years, function(yr) sum(data_DS$year == yr))
  max_nind <- max(nind_per_year) + nz
  
  y <- matrix(0, nrow = max_nind, ncol = n_years)
  d <- matrix(B + 1, nrow = max_nind, ncol = n_years)  # B+1 for pseudo-individuals
  site <- matrix(NA, nrow = max_nind, ncol = n_years)
  
  for (t in 1:n_years) {
    observed <- data_DS$year == t
    num_obs <- sum(observed)
    
    y[1:num_obs, t] <- data_DS$N_obs[observed]
    d[1:num_obs, t] <- ifelse(is.na(data_DS$d[observed]), B + 1, data_DS$d[observed])
    site[1:num_obs, t] <- data_DS$site[observed]
    
    # Fill pseudo-individuals
    y[(num_obs + 1):max_nind, t] <- 0
    d[(num_obs + 1):max_nind, t] <- runif(max_nind - num_obs, 0, B)
    site[(num_obs + 1):max_nind, t] <- sample(1:nsites, max_nind - num_obs, replace = TRUE)
  }
  
  # Bundle and summarize data set
  nimble.constants <- list(
    nyears = n_years,
    nsites = nsites,
    B = B,
    transect_len = transect_len,
    size_study_area = size_study_area,
    max_nind = max_nind
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
      psi[t] <- sum(lambda[1:nsites, t]) / max_nind
      
      # 'Likelihood'
      for (i in 1:max_nind) {
        # i is index for individuals
        z[i, t] ~ dbern(psi[t])                           # Data augmentation variables
        d[i, t] ~ dunif(0, B + 1)                         # Allow range up to B + 1
        is_observed[i, t] <- (d[i, t] <= B)               # Indicator for valid observations
        p[i, t] <- is_observed[i, t] * exp(-d[i, t] * d[i, t] / (2 * sigma * sigma))
        y[i, t] ~ dbern(z[i, t] * p[i, t])                # Basic Bernoulli random variable
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
    beta0 <- rep(0, n_years)
    lambda <- matrix(exp(beta0), nrow = nsites, ncol = n_years)
    is_observed <- (d <= B) 
    p <- is_observed * exp(-d * d / (2 * sigma * sigma))
    
    list(
      beta0 = beta0,
      alpha0 = alpha0,
      sigma = sigma,
      z = y,
      N = matrix(5, nrow = nsites, ncol = n_years),
      p = p,
      is_observed = is_observed
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

# # Analysis of DS data across all years -----------------------------
# DS_out2 <- DS_runModel(
#   data_DS = DS_data,
#   nsites = n_sites,
#   transect_len = transect_len,
#   nz = 200
# )
# 
# MCMCvis::MCMCtrace(DS_out2)
# 
# N_estimate_tibble2 <- map(DS_out2, as_tibble) %>% 
#   bind_rows() %>% 
#   select(starts_with("N_gic[")) %>% 
#   janitor::clean_names() %>% 
#   pivot_longer(everything(), names_to = "year", values_to = "N") %>% 
#   mutate(year = as.numeric(str_remove(year, "n_gic_"))) %>% 
#   group_by(year) %>% 
#   summarise(mean = mean(N),
#             q2.5 = quantile(N, probs = 0.025),
#             q97.5 = quantile(N, probs = 0.975)) %>% 
#   bind_rows(tibble(year = 1:n_years,
#                    mean = colSums(N)), .id = "id")
# 
# N_estimate_tibble2 %>% 
#   ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
#   geom_line() +
#   geom_ribbon(aes(ymin = c(q2.5),
#                   ymax = c(q97.5)),
#               linetype=2, alpha=0.1) +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, 500)) +
#   labs(x = "Year", y = "N")
