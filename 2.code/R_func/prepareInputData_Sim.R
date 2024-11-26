prepareInputData_Sim <- function(SimData){
  
  n_areas <- 1
  n_sites <- max(DS_data$site)
  n_years <- max(DS_data$year)
  
  DS_count <- DS_data %>% 
    as_tibble() %>% 
    select(-c(u, v)) %>% 
    group_by(year, site) %>% 
    summarise(count = sum(y, na.rm = T)) %>% 
    ungroup() %>% 
    pull(count) %>% 
    array(dim = c(n_areas, n_sites, n_years))
    
    L <- SimData$DS.data$L
    y <- SimData$DS.data$d
    Year_obs <- SimData$DS.data$d_year
    sumR_obs <- SimData$Rep.data$sumR_obs
    sumAd_obs <- SimData$Rep.data$sumAd_obs
    sumR_obs_year <- SimData$Rep.data$sumR_obs_year
    N_sites <- SimData$SimParams$Jmax
    N_obs <- length(SimData$DS.data$d)
    N_sumR_obs <- SimData$Rep.data$N_sumR_obs
  
  
  ## Reformat data into lists for Nimble
  input_data <- list(
    nim.data = list(
      d = d,
      L = L,
      DS_count = DS_count
    ),
    
    nim.constants = list(
      n_areas = 1,
      n_years = 
      
    )
  )
  
  ## Return data
  return(input_data)
}