#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-28
#'
#' Script description:
#' @data_DS: distance-sampling data
#' @dist_max: distance maximale at which an observation can be made
#' @transect_len: length of the transects (identical for all transects)
#' @size_study_area: size of the entire study area considered 
#' @nsites: number of prospected sites
#' @nz: number of generated pseudo-individuals for data augmentation
#' -------------------------------------------------------------------------------


run_DS_model <- function(data_DS,
                         dist_max = 0.6,
                         transect_len = 2,
                         size_study_area = 250,
                         nsites = 50,
                         nz = 200){
  B <- dist_max
  # data_DS <- DS_data %>% filter(year == 1)
  # Data augmentation: add "pseudo-individuals"
  nind <- nrow(data_DS)
  y <- c(data_DS$N_obs, rep(0, nz))     # Augmented detection indicator y
  
  # Modification: Augment 'site' with random integers for pseudo-individuals
  site <- c(data_DS$site, sample(1:nsites, nz, replace = TRUE))
  
  # Replace NA in d with B + 1 for pseudo-individuals
  d <-  c(ifelse(is.na(data_DS$d), B + 1, data_DS$d), runif(nz, 0, B))     # Use B + 1 to indicate unobservable
  
  # Bundle and summarize data set
  nimble.constants <- list(
    nsites = nsites,
    B = B,
    transect_len = transect_len,
    size_study_area = size_study_area,
    nind = nind,
    nz = nz
  )
  
  nimble.data <- list(y = y, d = d, site = site)
  
  # Nimble model for line transect HDS (NOT point transects!)
  nimble_code <- nimbleCode({
    # Prior distributions
    beta0 ~ dunif(-2, 10)   # Intercept of lambda
    alpha0 ~ dunif(-10, 10)  # Intercept of log(sigma)
    log(sigma) <- alpha0
    
    # psi is a derived parameter under DA for stratified populations
    psi <- sum(lambda[1:nsites]) / (nind + nz)
    
    # 'Likelihood'
    for (i in 1:(nind + nz)) {
      # i is index for individuals
      z[i] ~ dbern(psi)                               # Data augmentation variables
      d[i] ~ dunif(0, B + 1)                          # Allow range up to B + 1
      is_observed[i] <- (d[i] <= B)                  # Indicator for valid observations
      p[i] <- is_observed[i] * exp(-d[i] * d[i] / (2 * sigma * sigma))
      y[i] ~ dbern(z[i] * p[i])                       # Basic Bernoulli random variable
      site[i] ~ dcat(site.probs[1:nsites])            # Population distribution among sites
    }
    
    # Linear models for abundance and detection
    for (s in 1:nsites) {
      # s is index for sites
      N[s] ~ dpois(lambda[s])                         # Realized abundance at site s
      log(lambda[s]) <- beta0                         # same lambda at all site
      site.probs[s] <- lambda[s] / sum(lambda[1:nsites])
    }
    
    # Derived parameter: total population size across all sites
    Ntotal <- sum(z[1:(nind + nz)])
    area <- nsites * 2 * B * transect_len
    D <- Ntotal / area
    N_gic <- D * size_study_area
  })
  
  # Inits
  inits <- function() {
    alpha0 <- 0
    sigma <- exp(alpha0)
    beta0 <- 0
    lambda <- rep(exp(beta0), nsites)
    
    list(
      beta0 = beta0,
      alpha0 = alpha0,
      sigma = sigma,
      z = y,
      N = rep(5, nsites)
    )
  }
  
  # Parameters to save
  params <- c("sigma", "beta0", "psi", "Ntotal", "D", "N_gic")
  
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



# Analysis of DS data independently for each year -----------------------------
# N_estimate_DS <- N_estimates_upper <- N_estimates_lower <- N_estimates2 <- c()
# 
# for (i in 1:n_years){
#   print(i/n_years)
#   DS_out <- run_DS_model(data_DS = DS_data %>% filter(year == i) %>% select(-year),
#                          nsites = n_sites,
#                          transect_len = transect_len,
#                          nz = 200)
#   N_estimate_DS <- c(N_estimate_DS, map(DS_out, as_tibble) %>% bind_rows() %>% pull(N_gic) %>% mean())
#   N_estimates_upper <- c(N_estimates_upper, map(DS_out, as_tibble) %>% bind_rows() %>% pull(N_gic) %>% quantile(probs = 0.025))
#   N_estimates_lower <- c(N_estimates_lower, map(DS_out, as_tibble) %>% bind_rows() %>% pull(N_gic) %>% quantile(probs = 0.975))
# }
# 
# N_estimate_tibble <- tibble(year = 1:n_years,
#                             N_real = colSums(N),
#                             N_DS = N_estimate_DS,
#                             q97.5 = N_estimates_upper,
#                             q2.5 = N_estimates_lower) %>%
#   pivot_longer(-c(year, q2.5, q97.5),
#                names_to = "id", values_to = "mean") %>%
#   mutate(q97.5 = ifelse(id == "N_DS", q97.5, NA),
#          q2.5 = ifelse(id == "N_DS", q2.5, NA)) 
# 
# N_estimate_tibble %>%
#   ggplot(aes(x = year, y = mean, group = as.factor(id), colour = as.factor(id))) +
#   geom_line() +
#   geom_ribbon(aes(ymin = c(q2.5),
#                   ymax = c(q97.5)),
#               linetype=2, alpha=0.1) +
#   theme_minimal() +
#   coord_cartesian(ylim = c(0, 500))
# 


