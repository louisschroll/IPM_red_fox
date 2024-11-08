#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-28
#'
#' Script description:
#' @DS_data: distance-sampling data
#' @dist_max: distance maximale at which an observation can be made
#' @transect_len: length of the transects (identical for all transects)
#' @size_study_area: size of the entire study area considered 
#' @nsites: number of prospected sites
#' @nz: number of generated pseudo-individuals for data augmentation
#' -------------------------------------------------------------------------------


run_DS_model <- function(DS_data,
                         dist_max = 0.6,
                         transect_len = 2,
                         size_study_area = 250,
                         nsites = 50,
                         nz = 200){
  data <- DS_data
  B <- dist_max

  # Data augmentation: add "pseudo-individuals"
  nind <- nrow(data)
  y <- c(data[, 2], rep(0, nz))     # Augmented detection indicator y
  site <- c(data[, 1], rep(NA, nz)) # Augmented site indicator,
  # unknown (i.e., NA) for augmented inds.
  d <- c(data[, 5], rep(NA, nz))    # Augmented distance data (with NAs)
  
  # Bundle and summarize data set
  win.data <- list(
    nsites = nsites,
    B = B,
    transect_len = transect_len,
    size_study_area = size_study_area,
    nind = nind,
    nz = nz,
    y = y,
    d = d,
    site = site
  )
  

  # JAGS model for line transect HDS (NOT point transects!)
  cat(
    "
model{
  # Prior distributions
  beta0 ~ dunif(0,10)   # Intercept of lambda
  alpha0 ~ dunif(-10,10)  # Intercept of log(sigma)
  log(sigma) <- alpha0

  # psi is a derived parameter under DA for stratified populations
  psi <- sum(lambda[1:nsites]) / (nind + nz)

  # 'Likelihood'
  for(i in 1:(nind + nz)){                          # i is index for individuals
    z[i] ~ dbern(psi)                               # Data augmentation variables
    d[i] ~ dunif(0, B)                              # distance uniformly distributed
    p[i] <- exp(-d[i] * d[i] / (2 * sigma * sigma)) # Detection function
    y[i] ~ dbern(z[i] * p[i])                       # Basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites])            # Population distribution among sites
  }

  # Linear models for abundance and for detection
  for(s in 1:nsites){                               # s is index for sites
    # Model for abundance
    N[s] ~ dpois(lambda[s])                         # Realized abundance at site s
    log(lambda[s]) <- beta0                         # same lambda at all site
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])

    # detection
    # log(sigma[s]) <- alpha0
  }

  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind + nz)])
  area <- nsites * 2 * B * transect_len
  D <- Ntotal / area
  N_gic <- D * size_study_area
  
  # N_gic <- Ntotal / (nsites * 2 * B * transect_len) * size_study_area
}
", fill = TRUE, file = "model1.txt")
  
  
  # Inits
  zst <- c(rep(1, sum(y)), rep(0, nz)) # Initial values for DA variables
  inits <- function() {
    list(beta0 = 0,
         alpha0 = 0,
         z = zst)
  }
  
  # Parameters to save
  params <- c("sigma", "beta0", "psi", "Ntotal", "D", "N_gic")
  
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
