#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-29
#'
#' Script description:
#' Simulate age-structure data
#' @N: Population matrix with one column for each year and one line for each 
#' age class
#' @harvest_rate: mean proportion of animal killed, identical for each age-class
#' @output: similar to N, with the number of killed animal for each age-class
#' 
#' Example:
#' N <- sim_pop_dyn(Tmax = 10)
#' sim_age_data(N)
#' 
#' -------------------------------------------------------------------------------


sim_age_data <- function(N, harvest_rate = 0.2){
  # matrix(rpois(n = length(c(N)), 
  #              lambda = c(N) * harvest_rate),
  #        ncol = ncol(N))
  matrix(rbinom(n = length(c(N)), 
                size = c(N),
                prob = harvest_rate),
         ncol = ncol(N))
  
}

