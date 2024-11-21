# ---------
library(Distance)

library(Distance)

ds_data_year <- DS_data %>% filter(year == 4) %>% select(-c(year, u, v)) %>% mutate(distance = d)

# Example region table
region.table <- data.frame(
  Region.Label = "StudyArea",
  Area = 250
)

# Example sample table
sample.table <- data.frame(
  Sample.Label = unique(ds_data_year$site),  # Unique sample sites
  Region.Label = "StudyArea",
  Effort = 7  # Assuming each site had an effort of 1 unit
)

# Example observation table
obs.table <- data.frame(
  object = seq_len(nrow(ds_data_year)), # Unique observation IDs
  Region.Label = "StudyArea",
  Sample.Label = ds_data_year$site
)

# Fit detection function
detection_model <- ds(data = ds_data_year, 
                      transect = "line", 
                      key = "hn",
                      region_table = region.table,
                      sample_table = sample.table,
                      obs_table = obs.table)
summary(detection_model)
gof_ds(detection_model)

N_estimates2 <- c()
for (i in 1:(Tmax + 1)){
  ds_data_year <- DS_data %>% filter(year == i) %>% select(-c(year, u, v)) %>% mutate(distance = d)
  
  # Example region table
  region.table <- data.frame(
    Region.Label = "StudyArea",
    Area = 250
  )
  
  # Example sample table
  sample.table <- data.frame(
    Sample.Label = unique(ds_data_year$site),  # Unique sample sites
    Region.Label = "StudyArea",
    Effort = 4  # Assuming each site had an effort of 1 unit
  )
  
  # Example observation table
  obs.table <- data.frame(
    object = seq_len(nrow(ds_data_year)), # Unique observation IDs
    Region.Label = "StudyArea",
    Sample.Label = ds_data_year$site
  )
  
  # Fit detection function
  detection_model <- ds(data = ds_data_year, 
                        transect = "line", 
                        key = "hn",
                        region_table = region.table,
                        sample_table = sample.table,
                        obs_table = obs.table)
  N_estimates2 <- c(N_estimates2, detection_model$dht$individuals$N$Estimate)
  # N_estimates2 <- c(N_estimates2, out1$mean$N_gic2)
  # N_estimates_upper <- c(N_estimates_upper, out1$q2.5$N_gic)
  # N_estimates_lower <- c(N_estimates_lower, out1$q97.5$N_gic)
}