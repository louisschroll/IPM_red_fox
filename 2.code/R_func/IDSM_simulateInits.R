#' Simulate values for initializing integrated model
#'
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#'
#' @return A list containing one complete set of initial values for the model.
#' @export
#'
#' @examples

IDSM_simulateInits <- function(nim.data,
                          nim.constants) {
 # set.seed(initVals.seed)
  
  # Limits and constants
  n_age_class <- nim.constants$n_age_class
  n_years <- nim.constants$n_years
  n_sites <- nim.constants$n_sites
  
  size_hunting_area <- nim.constants$size_hunting_area
  L <- nim.data$L
  W <- nim.constants$W
  pi <- 3.141593
  
  # Vital rates
  ## survival parameters
  mean.survival <- runif(1, 0.25, 0.45)

  ## reproductive parameters
  mean.recruitment <- runif(1, 1.5, 3)
  
  # Detection parameters #
  #----------------------#
  
  ## Area-specific detection parameters
  log.mean.sigma <- runif(1, 3.5, 5.5)
  sd.det.area <- runif(1, 0.05, 0.2)
  
  sigma <- esw <- p <- c()
  
  sigma[1:n_years] <- exp(log.mean.sigma)
  esw[1:n_years] <- sqrt(pi * sigma[1:n_years] ^ 2 / 2)
  p[1:n_years] <- min(esw[1:n_years], W) / W
  
  # Population model #
  #------------------#
  
  ## Initial densities / population sizes
  sigma.D <- runif(1, 0.1, 2)
  
  N_exp <- Density <- array(0, dim = c(n_age_class, n_sites, n_years))

  D_x_sum <- nim.data$DS_count[ , ] / (L[ , ] * W * 2)
  D_data <- D_x_sum[which(!is.na(D_x_sum) & D_x_sum > 0)]
  Mu.D1 <- runif(1, quantile(D_data, 0.25), quantile(D_data, 0.75))
  
  for (j in 1:n_sites) {
    Density[2:n_age_class, j, 1] <- runif(n_age_class-1, 0, 1)
    Density[1, j, 1] <- sum(Density[2:n_age_class, j, 1]) / 2 * mean.recruitment
    N_exp[1:n_age_class, j, 1] <- Density[1:n_age_class, j, 1] * L[j, 1] * W * 2
  }
  
  ## Population projection over time
  for (j in 1:n_sites) {
    for (t in 2:n_years) {
      for(a in 2:(n_age_class-1)){
        Density[a, j, t] <- Density[a - 1, j, t - 1] * mean.survival
      }		
      Density[n_age_class, j, t] <- (Density[n_age_class-1, j, t-1] + Density[n_age_class, j, t-1]) * mean.survival
      Density[1, j, t] <- sum(Density[2:n_age_class, j, t]) / 2 * mean.recruitment
      N_exp[1:n_age_class, j, t] <- Density[1:n_age_class, j, t] * L[j, t] * W * 2
    }
  }

  
  ## Area-specific population size and density
  N_tot_exp <- matrix(NA, nrow = n_sites, ncol = n_years)
  
    for (t in 1:n_years) {
      N_tot_exp[t] <- sum(N_exp[1, 1:n_sites, t] + N_exp[2, 1:n_sites, t])    ## Summing up expected number of birds in covered area;
    }
  
  
  ## Area-, year-, and age-class specific density (for monitoring)
  meanDens <- N_tot_gic <- array(NA, dim = c(n_age_class, n_years))
  
    for (a in 1:n_age_class) {
      for (t in 1:n_years) {
        meanDens[a, t] <- mean(Density[a, 1:n_sites, t])
        N_tot_gic[a, t] <- round(meanDens[a, t] * size_hunting_area)
      }
    }
  
  # Age-at-harvest parameters
  harvest_rate <- runif(1, 0, 1)
  C <- matrix(rbinom(n = n_age_class * n_years, 
                     size = N_tot_gic, 
                     prob = harvest_rate),
              ncol = n_years)
  
  
  # Assembly #
  InitVals <- list(
    mean.recruitment = mean.recruitment,
    
    log.mean.sigma = log.mean.sigma,
    sigma = sigma,
    sigma2 = sigma ^ 2,
    esw = esw,
    p = p,
    
    mean.survival = mean.survival,
    
    Density = Density,
    meanDens = meanDens,
    N_exp = N_exp,
    N_tot_exp = N_tot_exp,
    N_tot_gic = N_tot_gic,
    
    harvest_rate = harvest_rate,
    C = C
  )
  
  return(InitVals)
}