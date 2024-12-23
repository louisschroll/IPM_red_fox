#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/application_ds_model.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-28
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------

# cat("\014")              # clear the console
# rm(list = ls())          # remove all variables of the work space


run_DS_model <- function(DS_data){
  data <- DS_data$data
  nsites <- DS_data$nsites
  B <- DS_data$dist_max
  transect_len <- DS_data$transect_len
  size_study_area <- DS_data$size_study_area
  
  # Data augmentation: add a bunch of "pseudo-individuals"
  nz <- 200                        # number of pseudo-individuals
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
  

  # BUGS model for line transect HDS (NOT point transects!)
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
  for(i in 1:(nind + nz)){                 # i is index for individuals
    z[i] ~ dbern(psi)                      # Data augmentation variables
    d[i] ~ dunif(0, B)                     # distance uniformly distributed
    p[i] <- exp(-d[i] * d[i] / (2 * sigma * sigma)) # Detection function
    y[i] ~ dbern(z[i] * p[i])                    # Basic Bernoulli random variable
    site[i] ~ dcat(site.probs[1:nsites])   # Population distribution among sites
  }

  # Linear models for abundance and for detection
  for(s in 1:nsites){                      # s is index for sites
    # Model for abundance
    N[s] ~ dpois(lambda[s])                # Realized abundance at site s
    log(lambda[s]) <- beta0                # same lambda at all site
    site.probs[s] <- lambda[s] / sum(lambda[1:nsites])

    # detection
    # log(sigma[s]) <- alpha0
  }

  # Derived parameter: total population size across all sites
  Ntotal <- sum(z[1:(nind + nz)])
  area <- nsites * 2 * B * transect_len
  D <- Ntotal / area
  N_gic <- D * size_study_area
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
  params <- c("alpha0", "beta0", "psi", "Ntotal", "D", "N_gic")
  
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


