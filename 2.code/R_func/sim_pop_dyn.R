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
  # Define mean of the demographic parameters
  mean.sj <- 0.54 # 0.37
  mean.sa <- 0.54
  mean.f1 <- 0 #1.4
  mean.fa <- 1.578
  surv_rate <- c(mean.sj, rep(mean.sa, n_age_class-1))
  fec_rate <- rep(mean.fa, n_age_class-1)
  
  # Define population matrix and initial stage-specific population sizes
  N <- matrix(NA, nrow = n_age_class, ncol = n_years)
  colnames(N) <- paste0("Y", 1:n_years)
  R <- mean.fa
  S <- mean.sa
  M <- matrix(c(0, R, R, R, R,
                S, 0, 0, 0, 0,
                0, S, 0, 0, 0,
                0, 0, S, 0, 0,
                0, 0, 0, S, S), ncol = n_age_class) %>% t()
  
  eignevalue <- eigs(M, "ss")
  stable_stage_prop <- eignevalue / sum(eignevalue)
  N[1:n_age_class, 1] <- round(N0 * stable_stage_prop)

  # Project population
  for (t in 1:(n_years - 1)) {
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors
    for (a in 1:(n_age_class - 2)) {
      N[a + 1, t + 1] <- N[a, t] * surv_rate[a] #rbinom(1, N[a, t], surv_rate[a])
    }
    
    # Age class 4+ (index = n_age_class = 5): age class 4 and 5+ survivors
    N[n_age_class, t + 1] <- (N[n_age_class - 1, t] + N[n_age_class, t])*surv_rate[n_age_class] #rbinom(1, N[n_age_class - 1, t] + N[n_age_class, t], surv_rate[n_age_class])
    
    # Age class 0 (index = 1): local reproduction
    N[1, t + 1] <- surv_rate[1] * sum(fec_rate * N[2:n_age_class, t + 1]) #rpois(1, surv_rate[1] * sum(fec_rate * N[2:n_age_class, t + 1]))
    
    if (sum(N[, t + 1]) == 0)
      break # Stop calculation if pop. extinct
    
  }
  return(round(N))
}
