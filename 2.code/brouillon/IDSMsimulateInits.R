#' Simulate values for initializing integrated model
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param survVarT logical. If TRUE, survival is simulated including annual variation.
#'
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

simulateInits <- function(nim.data,
                          nim.constants,
                          survVarT,
                          initVals.seed = 42) {
  set.seed(initVals.seed)
  
  # Limits and constants
  n_areas <- nim.constants$n_areas
  n_age_class <- nim.constants$n_age_class
  n_years <- nim.constants$n_years
  n_sites <- nim.constants$n_sites
  
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  
  # Vital rates
  ## Area-specific survival parameters
  mean.survival <- runif(1, 0.40, 0.45)
  h.sigma.S <- runif(1, 0.05, 0.2)

  mu.S <- rnorm(n_areas,
                qlogis(mean.survival),
                sd = h.sigma.S)
  
  sigmaT.S <- runif(1, 0.05, 0.2)
  sigmaR.S <- runif(1, 0.05, 0.2)
  
  Mu.S <- rep(NA, n_areas)
  S <-  matrix(NA, nrow = n_areas, ncol = n_years - 1)

  epsT.S <- rep(0, n_years - 1)
  #epsT.S <- rnorm(n_years-1, 0, sigmaT.S)
  epsR.S <- matrix(0, nrow = n_areas, ncol = n_years-1)
  #epsR.S <- matrix(rnorm(n_areas*n_years, 0, sigmaR.S), nrow = n_areas, ncol = n_years)
  
  for (x in 1:n_areas) {
    Mu.S[x] <- plogis(mu.S[x])
    S[x, 1:(n_years-1)] <- plogis(qlogis(Mu.S[x]) + epsT.S[1:(n_years-1)] + epsR.S[x, 1:(n_years-1)])
  }
  
  ## Area-specific reproductive parameters
  h.Mu.R  <- runif(1, 1.5, 3)
  h.sigma.R <- runif(1, 0.05, 0.2)
  
  h.Mu.betaR.R <- runif(1, 0.01, 0.1)
  h.sigma.betaR.R <- runif(1, 0.05, 0.1)
  
  Mu.R <- rlnorm(n_areas, meanlog = log(h.Mu.R), sdlog = h.sigma.R)
  
  sigmaT.R <- runif(1, 0.05, 0.2)
  sigmaR.R <- runif(1, 0.05, 0.2)
  
  R_year <- matrix(NA, nrow = n_areas, ncol = n_years)
  
  epsT.R <- rep(0, n_years)
  #epsT.R <- rnorm(N_year, 0, sigmaT.R)
  epsR.R <- matrix(0, nrow = n_areas, ncol = n_years)
  #epsR.R <- matrix(rnorm(n_areas*n_years, 0, sigmaR.R), nrow = n_areas, ncol = n_years)
  
  for (x in 1:n_areas) {
    R_year[x, 1:n_years] <- exp(log(Mu.R[x]) + betaR.R[x] * RodentOcc[x, 1:n_years] + epsT.R[1:n_years] + epsR.R[x, 1:n_years])
  }
  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  log.mean.sigma <- runif(1, 3.5, 5.5)
  sd.det.area <- runif(1, 0.05, 0.2)
  
  #log.mean.sigma.area <- rnorm(n_areas, log.mean.sigma, sd = sd.det.area)
  log.mean.sigma.area <- rep(log.mean.sigma, n_areas)
  
  sd.det.year <- runif(1, 0.05, 0.2)
  sd.det.residual <- runif(1, 0.05, 0.2)
  
  sigma <- esw <- p <- matrix(NA, nrow = n_areas, ncol = n_years)
  
  epsT.det <- rep(0, n_years)
  #epsT.det <- rnorm(n_years, 0, sd = sd.det.year)
  epsR.det <- matrix(0, nrow = n_areas, ncol = n_years)
  #epsR.det <- matrix(rnorm(n_areas*n_years, 0, sd.det.residual), nrow = n_areas, ncol = n_years)
  
  for (x in 1:n_areas) {
    sigma[x, 1:n_years] <- exp(log.mean.sigma.area[x] + epsT.det[1:n_years] + epsR.det[x, 1:n_years])
    esw[x, 1:n_years] <- sqrt(pi * sigma[x, 1:n_years] ^ 2 / 2)
    p[x, 1:n_years] <- min(esw[x, 1:n_years], W) / W
  }
  
  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  Mu.D1 <- rep(NA, n_areas)
  sigma.D <- runif(n_areas, 0.1, 2)
  
  N_exp <- Density <- array(0, dim = c(n_areas, n_age_class, max(n_sites), n_years))
  
  for (x in 1:n_areas) {
    D_x_sum <- nim.data$N_a_line_year[x, 2, , ] / (L[x, , ] * W * 2)
    D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
    Mu.D1[x] <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))
    
    for (j in 1:n_sites[x]) {
      Density[x, 2, j, 1] <- Mu.D1[x]
      Density[x, 1, j, 1] <- Density[x, 2, j, 1] * R_year[x, 1] # Juveniles
      
      N_exp[x, 1, j, 1] <- Density[x, 1, j, 1] * L[x, j, 1] * W * 2
      N_exp[x, 2, j, 1] <- Density[x, 2, j, 1] * L[x, j, 1] * W * 2
    }
  }
  
  ## Population projection over time
  for (x in 1:n_areas) {
    for (j in 1:n_sites[x]) {
      for (t in 2:n_years) {
        Density[x, 2, j, t] <- sum(Density[x, 1:n_age_class, j, t-1]) * S[x, t-1] # Adults
        Density[x, 1, j, t] <- Density[x, 2, j, t] * R_year[x, t] # Juveniles
        
        N_exp[x, 1:n_age_class, j, t] <- Density[x, 1:n_age_class, j, t] * L[x, j, t] * W * 2
      }
    }
  }
  
  ## Area-specific population size and density
  N_tot_exp <- matrix(NA, nrow = n_areas, ncol = n_years)
  
  for (x in 1:n_areas) {
    for (t in 1:n_years) {
      N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:n_sites[x], t] + N_exp[x, 2, 1:n_sites[x], t])    ## Summing up expected number of birds in covered area;
    }
  }
  
  ## Area-, year-, and age-class specific density (for monitoring)
  meanDens <- array(NA, dim = c(n_areas, n_age_class, n_years))
  
  for (x in 1:n_areas) {
    for (a in 1:n_age_class) {
      for (t in 1:n_years) {
        meanDens[x, a, t] <- mean(Density[x, a, 1:n_sites[x], t])
      }
    }
  }
  
  meanDens2 <- array(NA, dim = c(n_areas, n_age_class, n_years))
  
  # Apply over the first three dimensions (areas, age classes, and years)
  meanDens2 <- apply(expand.grid(1:n_areas, 1:n_age_class, 1:n_years), 1, function(indices) {
    x <- indices[1]  # area
    a <- indices[2]  # age class
    t <- indices[3]  # year
    mean(Density[x, a, 1:n_sites[x], t])
  }) %>%
    array(dim = c(n_areas, n_age_class, n_years))
  
  
  # Assembly #
  InitVals <- list(
    Mu.D1 = Mu.D1,
    sigma.D = sigma.D,
    eps.D1 = matrix(0, nrow = nim.constants$n_areas, ncol = max(n_sites)),
    
    Mu.R = Mu.R,
    h.Mu.betaR.R = h.Mu.betaR.R,
    h.sigma.betaR.R = h.sigma.betaR.R,
    h.Mu.R = h.Mu.R,
    h.sigma.R = h.sigma.R,
    sigmaT.R = sigmaT.R,
    sigmaR.R = sigmaR.R,
    epsT.R = epsT.R,
    epsR.R = epsR.R,
    epsA.R =  log(Mu.R) - log(h.Mu.R),
    R_year = R_year,
    
    log.mean.sigma.area = log.mean.sigma.area,
    log.mean.sigma = log.mean.sigma,
    sd.det.area = sd.det.area,
    sd.det.year = sd.det.year,
    sd.det.residual = sd.det.residual,
    epsT.det = epsT.det,
    epsR.det = epsR.det,
    eps.det.area = log.mean.sigma.area - log.mean.sigma,
    sigma = sigma,
    sigma2 = sigma ^ 2,
    esw = esw,
    p = p,
    
    mean.survival = mean.survival,
    h.sigma.S = h.sigma.S,
    mu.S = mu.S,
    Mu.S = Mu.S,
    sigmaT.S = sigmaT.S,
    sigmaR.S = sigmaR.S,
    epsT.S = epsT.S,
    epsR.S = epsR.S,
    epsA.S = mu.S - logit(mean.survival),
    S = S,
    
    Density = Density,
    meanDens = meanDens,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp
  )
  
  return(InitVals)
}