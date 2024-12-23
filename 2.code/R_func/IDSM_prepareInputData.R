

IDSM_prepareInputData <- function(DS_data, 
                                  harvest_data,
                                  W, 
                                  size_hunting_area,
                                  beta_param){

  n_sites <- max(DS_data$site)
  n_years <- max(DS_data$year)

  DS_count_tbl <- DS_data %>%
    group_by(year, site) %>%
    summarise(count = sum(N_obs, na.rm = T),
              L = mean(effort)) %>%
    ungroup() 
  
  DS_count <- DS_count_tbl %>%
    pull(count) %>%
    matrix(nrow = n_sites, ncol = n_years)
  
  L <- DS_count_tbl %>%
    pull(L) %>%
    matrix(nrow = n_sites, ncol = n_years)

  # retrieve distance data
  # n_obs_dist: numeric vector, number of observed distances at each study area
  # d: numeric matrix, nrow = n_areas and ncol = max(n_obs_dist), observed distances
  # Year_obs: matrix, same dim than d, year of observed distances in d
  d <- DS_data %>%
    filter(!is.na(d)) %>%
    pull(d)
  
  year_obs <- DS_data %>%
    filter(!is.na(d)) %>%
    pull(year)
  
  # # Prior for initial density
  # N0 <- 150
  # R <- 1.6
  # S <- 0.56
  # M <- matrix(c(0, R, R, R, R,
  #               S, 0, 0, 0, 0,
  #               0, S, 0, 0, 0,
  #               0, 0, S, 0, 0,
  #               0, 0, 0, S, S), ncol = n_age_class) %>% t()
  # 
  # eignevalue <- eigs(M, "ss")
  # stable_stage_prop <- eignevalue / sum(eignevalue)
  # intit_density <- round(N0 * stable_stage_prop) / size_hunting_area
  # mu <- intit_density
  # sigma <- 0.1
  # 
  # a <- ((1 - mu) / (sigma * sigma) - (1 / mu)) * mu ^ 2
  # b <- a * ( (1/mu) - 1)
  # 
  # beta_param <- matrix(c(a, b), ncol = 2)

  ## Reformat data into lists for Nimble
  input_data <- list(
    nim.data = list(
      d = d,
      L = L,
      DS_count = DS_count,
      C = harvest_data
    ),
    
    nim.constants = list(
      n_sites = n_sites,
      n_years = n_years,
      n_obs_dist = length(d),
      n_age_class = nrow(harvest_data),
      #year_obs = year_obs,
      W = W,
      size_hunting_area = size_hunting_area,
      beta_param = beta_param,
      pi = 3.141593
    ))

  return(input_data)
}



