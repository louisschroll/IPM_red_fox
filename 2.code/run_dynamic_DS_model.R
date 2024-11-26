
run_dynamic_DS_model <- function(DS_data,
                                 dist_max = 0.6,
                                 transect_len = 2,
                                 size_study_area = 250,
                                 nz = 200,
                                 age_max = 5){
  B <- dist_max
  nsites <- max(DS_data$site, na.rm = T)
  nyears <- length(unique(DS_data$year))
  npseudo_ind <- nz * nyears
  # Data augmentation: add "pseudo-individuals"
  nind <- nrow(DS_data)
  y <- c(DS_data$y, rep(0, npseudo_ind))        # Augmented detection indicator y
  site <- c(DS_data$site, rep(NA, npseudo_ind)) # Augmented site indicator,
  # unknown (i.e., NA) for augmented inds.
  d <- c(DS_data$d, rep(NA, npseudo_ind))    # Augmented distance data (with NAs)
  
  year <- c(DS_data$year, rep(unique(DS_data$year), each = nz))
  
  # mean.sj <- 0.37
  # mean.sa <- 0.56
  # mean.f1 <- 1.4
  # mean.fa <- 1.6
  # surv_rate <- c(mean.sj, rep(mean.sa, age_max))
  # fec_rate <- c(mean.f1, rep(mean.fa, 4))
  
  # Bundle and summarize data set
  win.data <- list(
    nsites = nsites,
    nyears = nyears,
    B = B,
    transect_len = transect_len,
    size_study_area = size_study_area,
    nind = nind,
    npseudo_ind = npseudo_ind,
    y = y,
    d = d,
    site = site,
    year = year
    # age_max = age_max,
    # surv_rate = surv_rate,
    # fec_rate = fec_rate
  )
  
  
  # JAGS model for line transect HDS (NOT point transects!)
  cat("model{
    ## Distance-sampling
    # Prior distributions
    for (t in 1:nyears) {  
      beta0[t] ~ dunif(0,10)   # Intercept of lambda
    }
    alpha0 ~ dunif(-10,10)  # Intercept of log(sigma)
    log(sigma) <- alpha0
    
    # Likelihood
    for (t in 1:nyears) {                          # t is index for years
      for (s in 1:nsites) {                        # s is index for sites
        
        # Model for abundance
        area_site[t, s] <- 2 * B * transect_len
        lambda[t, s] <- beta0[t]
        N_site[t, s] ~ dpois(lambda[t, s])               # Realized abundance at site s in year t
        site.probs[t, s] <- lambda[t, s] / sum(lambda[t, 1:nsites])  # Site probability
  
        # Detection model
        # log(sigma[t, s]) <- alpha0                  # Detection parameter (sigma) can vary by year and site
      }
    }
  
    # psi is a derived parameter under DA for stratified populations
    psi <- sum(lambda[1:nyears, 1:nsites]) / (nind + npseudo_ind)
      
    for (i in 1:(nind + npseudo_ind)) {                      # i is index for individuals
      z[i] ~ dbern(psi)                             # Data augmentation variables
      # year[i] ~ dcat(year.probs[1:nyears])          # Year assignment for each individual ??
      site[i] ~ dcat(site.probs[year[i], 1:nsites]) # Site assignment based on year-specific site probabilities
      d[i] ~ dunif(0, B)                            # Distance uniformly distributed
      p[i] <- exp(-d[i] * d[i] / (2 * sigma * sigma))  # Detection function sigma[year[i], site[i]]
      y[i] ~ dbern(z[i] * p[i])                     # Bernoulli random variable
    }
  
    # Derived parameter for total population size across all sites
    # Ntotal_ds <- sum(z[1:(nind + npseudo_ind)])
    size_area_sampled <- nsites * 2 * B * transect_len
    # D <- Ntotal_ds / size_area_sampled
    # N_gic <- D * size_study_area
    for (t in 1:nyears){
      Ntot[t] <- sum(N_site[t, 1:nsites])
      N_gic[t] <- Ntot[t] / size_area_sampled * size_study_area
    }
  }
", fill = TRUE, file = "model1.txt")
  
  
  # Inits
  inits <- function() {
    list(alpha0 = 0,
         beta0 = rep(0, nyears),
         z = y)
        # N = matrix(10, nrow = age_max, ncol = nyears+1)) #rpois(n = age_max * (nyears+1), lambda = c(120, 50, 40, 30, 20))
  }
  
  # Parameters to save
  params <- c("sigma", "psi", "Ntot", "D", "N_gic")
  
  # MCMC settings
  ni <- 12000
  nb <- 2000
  nt <- 1
  nc <- 3
  
  # Run the model in JAGS
  out1 <- jags(
    win.data,
    inits,
    params,
    "model1.txt",
    n.thin = nt,
    n.chains = nc,
    n.burnin = nb,
    n.iter = ni,
    parallel = TRUE
  )
  
  return(out1)
}



