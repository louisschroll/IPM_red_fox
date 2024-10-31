#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-21
#'
#' Script description:
#' @Tmax: number of years over which we let the population evolve
#' @age_max: greater age class
#' @N0: population size at the beginning of the study (t = 0)
#' -------------------------------------------------------------------------------


sim_pop_dyn <- function(Tmax = 10,
                        N0 = 200,
                        age_max = 5){
  
  # Define mean of the demographic parameters
  mean.sj <- 0.37
  mean.sa <- 0.56
  mean.f1 <- 1.4
  mean.fa <- 1.6
  surv_rate <- c(mean.sj, rep(mean.sa, age_max))
  fec_rate <- c(mean.f1, rep(mean.fa, 4))
  
  # Define population matrix and initial stage-specific population sizes
  N <- matrix(NA, nrow = 5, ncol = Tmax + 1)
  colnames(N) <- paste0("Y", 1:(Tmax+1))
  N[1:age_max, 1] <- round(N0 * c(0.5, 0.2, 0.13, 0.1, 0.07))
  
  # Project population
  for (t in 1:Tmax){
    # Age class 0 (index = 1): local reproduction
    N[1, t+1] <- rpois(1, surv_rate[1] * sum(fec_rate * N[1:age_max, t]))
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
    for(a in 1:(age_max-2)){
      N[a+1, t+1] <- rbinom(1, N[a, t], surv_rate[a])
    }		
    
    # Age class 4+ (index = age_max = 5): age class 4 and 5+ survivors
    N[age_max, t+1] <- rbinom(1, N[age_max-1, t] + N[age_max, t], surv_rate[age_max])
    
    if (sum(N[,t+1]) == 0) break # Stop calculation if pop. extinct
    
  }
  return(N)
}



