prepareInputData <- function(DS_data_list, 
                             harvest_data,
                             W, 
                             size_hunting_area){
  # 
  # n_areas <- 1
  # n_sites <- max(DS_data$site)
  # n_years <- max(DS_data$year)
  # 
  # DS_count <- DS_data %>% 
  #   group_by(year, site) %>% 
  #   summarise(count = sum(N_obs, na.rm = T)) %>% 
  #   ungroup() %>% 
  #   pull(count) %>% 
  #   array(dim = c(n_areas, n_sites, n_years))
  # 
  # DS_dist <- DS_data %>% 
  #   filter(!is.na(d)) %>% 
  #   pull(d)
  # 
  # ## Reformat data into lists for Nimble
  # input_data <- list(
  #   nim.data = list(
  #     d = DS_dist,
  #     L = DS_data$effort,
  #     DS_count = DS_count),
  #   
  #   nim.constants = list(
  #     n_areas = n_areas,
  #     n_sites = n_sites,
  #     n_years = n_years))
  
  n_areas <- length(simData_list)
  n_sites <- unlist(map(simData_list, . %>% pull(site) %>% max()))
  n_years <- unname(unlist(map(simData_list, . %>% pull(year) %>% max()))[1]) 
  
  
  sim_data_count <- simData_list %>% 
    bind_rows(.id = "area") %>% 
    complete(area, year, site) %>%
    group_by(area, year, site) %>%
    summarise(count = sum(N_obs, na.rm = T),
              L = min(effort)) %>%
    ungroup() 
  
  DS_count <- sim_data_count %>%  
    pull(count) %>%
    array(dim = c(n_areas, n_sites, n_years))
  
  L_array <- sim_data_count %>% 
    pull(L) %>%
    array(dim = c(n_areas, n_sites, n_years))
  
  # retrieve distance data
  # n_obs_dist: numeric vector, number of observed distances at each study area
  # d: numeric matrix, nrow = n_areas and ncol = max(n_obs_dist), observed distances
  # Year_obs: matrix, same dim than d, year of observed distances in d
  n_obs_dist <- unlist(map(simData_list, . %>% filter(!is.na(d)) %>% nrow()))
  
  dist_tbl <- simData_list %>%
    map(. %>% filter(!is.na(d)) %>%
          mutate(row_nb = row_number())) %>%
    bind_rows(.id = "area") %>%
    complete(area, row_nb) # note: NA are added for dimension purpose but are not used in the analysis
  
  d <- dist_tbl %>%
    pull(d) %>%
    matrix(nrow = n_areas)
  
  Year_obs <- dist_tbl %>%
    pull(year) %>%
    matrix(nrow = n_areas)
  
  
  ## Reformat data into lists for Nimble
  input_data <- list(
    nim.data = list(
      d = d,
      L = L_array,
      DS_count = DS_count,
      C = harvest_data
    ),
    
    nim.constants = list(
      n_areas = n_areas,
      n_sites = n_sites,
      n_years = n_years,
      n_obs_dist = n_obs_dist,
      n_age_class = nrow(harvest_data),
      W = W,
      size_hunting_area = size_hunting_area,
      pi = 3.141593
    ))
  
  return(input_data)
}



