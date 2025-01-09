#' HEADER ------------------------------------------------------------------------
#'
#' Script name:
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-21
#'
#' Script description:
#' @n_years: number of years over which we let the population evolve
#' @n_age_class: number of age class
#' @N0: population size at the beginning of the study (t = 0)
#' -------------------------------------------------------------------------------


sim_pop_dyn <- function(n_years = 10,
                        N0 = 200,
                        n_age_class = 5) {
  burnin_year <- 20
  n_years_tot <- n_years + burnin_year
  # Define mean of the demographic parameters
  mean.sj <- 0.56 # 0.37
  mean.sa <- 0.56
  mean.f1 <- 0 #1.4
  mean.fa <- 1.5
  surv_rate <- c(mean.sj, rep(mean.sa, n_age_class-1))
  fec_rate <- rep(mean.fa, n_age_class-1)
  pregnancy_rate <- 0.4
  rho <- 4
  # Define population matrix and initial stage-specific population sizes
  N <- B <- L <- NbRecruits <- matrix(NA, nrow = n_age_class, ncol = n_years_tot)
  # colnames(N) <- paste0("Y", 1:n_years_tot)
  # R <- mean.fa
  # S <- mean.sa
  # M <- matrix(c(R, R, R, R, R,
  #               S, 0, 0, 0, 0,
  #               0, S, 0, 0, 0,
  #               0, 0, S, 0, 0,
  #               0, 0, 0, S, S), ncol = n_age_class) %>% t()
  # 
  # eignevalue <- eigs(M, "ss")
  # stable_stage_prop <- eignevalue / sum(eignevalue)
  N[1:n_age_class, 1] <- round(N0 * c(0.4, 0.3, 0.1, 0.1, 0.1))

  # Project population
  for (t in 1:(n_years_tot - 1)) {
    # Age class 0 (index = 1): local reproduction
    for(a in 1:n_age_class){
      
      # Breeding Population Size: Number of females that reproduce
      B[a, t+1] <- rbinom(1, prob = pregnancy_rate, size = round(N[a, t] / 2))
      
      # Litter Size (in utero): Number of pups produced by females of age class a
      L[a, t+1] <- rpois(1, B[a, t+1] * rho)
      
      # Number Recruits: Number of pups surviving to emerge from the den
      NbRecruits[a, t+1] <- rbinom(1, prob = surv_rate[1], size = L[a, t+1])
    } 
    N[1, t + 1] <- sum(NbRecruits[1:n_age_class, t+1])
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors
    for (a in 1:(n_age_class - 2)) {
      N[a + 1, t + 1] <- rbinom(1, N[a, t], surv_rate[a])
    }
    
    # Age class 4+ (index = n_age_class = 5): age class 4 and 5+ survivors
    N[n_age_class, t + 1] <- rbinom(1, N[n_age_class - 1, t] + N[n_age_class, t], surv_rate[n_age_class])
    
    if (sum(N[, t + 1]) == 0)
      break # Stop calculation if pop. extinct
    
  }
  return(list(N = N[, (burnin_year+1):(n_years_tot)], 
              B = B[, (burnin_year+1):(n_years_tot)], 
              L = L[, (burnin_year+1):(n_years_tot)], 
              NbRecruits = NbRecruits[, (burnin_year+1):(n_years_tot)]))
}
sim_pop_dyn()
