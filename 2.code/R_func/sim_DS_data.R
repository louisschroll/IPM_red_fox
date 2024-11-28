#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/sim_DS_data.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-21
#'
#' Script description:
#' @sd.sigma.year: standard deviation of random year variation in detection probability
#' @sd.sigma.site: standard deviation of random line variation in detection probability
#'
#' -------------------------------------------------------------------------------


sim_DS_data <- function (type = c("line", "point"),
                         nsites = 50,
                         Ntot = c(200, 230, 190),
                         mean.sigma = 0.15,
                         sd.sigma.year = 0,
                         sd.sigma.site = 0,
                         dist_max = 0.6,
                         transect_len = 2,
                         size_study_area = 250,
                         discard0 = FALSE){
  type <- match.arg(type)
  n_years <- length(Ntot)
  size_area_sampled <- 2 * dist_max * transect_len
  
  lambda <- Ntot * size_area_sampled / size_study_area
  
  N <- matrix(NA, nrow = nsites, ncol = n_years)
  for (i in 1:n_years){
    N[1:nsites, i] <- rpois(nsites, lambda[i])
  }
  
  N_vect <- c(N)
  
  # Calculate year- and site-specific sigma parameter
  epsilonT <- rnorm(n_years, mean = 0, sd = sd.sigma.year)
  epsilonJ <- rnorm(nsites, mean = 0, sd = sd.sigma.site)
  sigma <- c(exp(log(mean.sigma) + outer(epsilonJ, epsilonT, FUN = '+')))
  
  # Simulate distance-sampling data
  tibble(year = rep(rep(1:n_years, each = nsites), N_vect),
         site = rep(rep(1:nsites, n_years), N_vect),
         sigma = rep(sigma, N_vect)) %>% 
    mutate(d = round(runif(n = length(site), min = 0, max = dist_max), 2),
           p = exp(-d * d / (2 * (sigma ^ 2))),
           N_obs = rbinom(nrow(.), 1, p)) %>% 
    filter(N_obs == 1) %>% 
    complete(year, site) %>% 
    mutate(N_obs = ifelse(is.na(N_obs), 0, N_obs),
           effort = transect_len) %>% 
    select(-c(sigma, p)) 
}


