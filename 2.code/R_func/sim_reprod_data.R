


sim_reprod_data <- function(pop_dyn_list){
  rho <- 4
  pregnancy_rate <- 0.4
  n_obs_female <- 500
  
  isBreeding <- rbinom(n_obs_female, size = 1, pregnancy_rate)
  
  nbFoetus <- rpois(sum(isBreeding), rho)
  
  return(list(isBreeding = isBreeding, nbFoetus = nbFoetus))
}





