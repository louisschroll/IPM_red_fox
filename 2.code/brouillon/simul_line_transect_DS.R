
library(jagsUI)
library(tidyverse)

simulate_DS_data <- function(N = 200, sigma = 30){
  #' Simulate non-hierarchical line transect data under CDS
  #' Function arguments:
  #'    N: number of individuals along transect with distance u(-100, 100)
  #'    sigma: scale parameter of half-normal detection function
  #' Function subjects N individuals to sampling, and then retains the value
  #' of x=distance only for individuals that are captured
  
  xall <- runif(N, -100, 100) # Distances of all N individuals
  g <- function(x, sig) exp(-x^2/(2*sig^2))
  p <- g(xall, sig=sigma) # detection probability
  y <- rbinom(N, 1, p) # some inds. are detected and their distance measured
  x <- xall[y==1]      # this has direction (right or left transect side)
  x <- abs(x)          # now it doesn't have direction
 
  return(list(nind = length(x), N = N, sigma = sigma, xall = xall, x = x))
}

# Obtain a data set for analysis
#set.seed(205)               # If you want to get same results
tmp <- simulate_DS_data(sigma = 30)

plot_simdata <- function(sim_data) {
  op <- par(mfrow = c(1, 2))
  
  with(sim_data, {
    # Plot the detection function
    curve(
      exp(-x ^ 2 / (2 * sigma ^ 2)), 0, 100,
      xlab = "Distance (x)",
      ylab = "Detection prob.",
      main = "Detection function",
      lwd = 2,
      ylim = c(0, 1))
    text(80, 0.9, paste("sigma:", sigma))
    
    # Histogram of all individual present
    hist(
      abs(xall),
      nclass = 10,
      xlab = "Distance (x)",
      col = "grey",
      main = "True (grey) \nand observed distances (blue)")
    # Histogram of observed individual
    hist(x, col = "blue", add = TRUE)
  })
  
  par(op)
}


plot_simdata(tmp)

run_DS_model <- function(DS_data){
  # Analysis of continuous data using data augmentation (DA)
  nind <- DS_data$nind
  nz <- 200 # Augment observed data with nz = 200 zeroes
  y <- c(rep(1, nind), rep(0, nz)) # Augmented inds. have y=0 by definition
  x <- c(DS_data$x, rep(NA, nz)) # Value of distance are missing for the augmented
  B = 100 #max(x, na.rm = TRUE) * 1.2 # B must superior than maximal distance
  
  # Bundle and summarize data set
  win.data <- list(
    nind = nind,
    nz = nz,
    x = x,
    y = y,
    B = B
  )
  
  # Save text file with BUGS model
  cat("
  model {
    # Priors
    sigma ~ dunif(0,1000)  # Half-normal scale
    psi ~ dunif(0,1)       # DA parameter

    # Likelihood
    for(i in 1:(nind+nz)){
      # Process model
      z[i] ~ dbern(psi)   # DA variables
      x[i] ~ dunif(0, B)  # Distribution of distances
      # Observation model
      logp[i] <- -((x[i] * x[i]) / (2 * sigma * sigma)) # Half-normal detection fct.
      p[i] <- exp(logp[i])
      mu[i] <- z[i] * p[i]
      y[i] ~ dbern(mu[i]) # Simple Bernoulli measurement error process
    }
    # Derived quantities
    N <- sum(z[1:(nind + nz)]) # Population size
    # D <- N / 60                # Density, with A = 60 km^2 when B = 500
  }
  ", fill=TRUE, file="model1.txt")
  
  # Inits
  zst <- y
  inits <- function(){ list (psi=runif(1), z=zst, sigma=runif(1,40,200)) }
  
  # Params to save
  params <- c("N", "sigma")
  
  # Experience the raw power of BUGS and summarize marginal posteriors
  
  out1 <- jags(win.data, 
               inits, 
               params, 
               "model1.txt", 
               n.thin = 1,
               n.chains = 3, 
               n.burnin = 1000, 
               n.iter = 11000,
               DIC = FALSE,
               parallel = TRUE) 
  return(out1)
}

out1 <- run_DS_model(DS_data = tmp)
print(out1, 3)
MCMCvis::MCMCtrace(out1, pdf = FALSE)


N_estimate <- sig_estimate <- c()
for (t in 1:500){
  tmp <- simulate_DS_data(N = 200, sigma = 30)
  out1 <- run_DS_model(DS_data = tmp)
  print(out1, 3)
  N_estimate <- c(N_estimate, unlist(out1$samples[,1]) %>% mean())
  sig_estimate <- c(sig_estimate, unlist(out1$samples[,2]) %>% mean())
}

boxplot(N_estimate)
boxplot(sig_estimate)
