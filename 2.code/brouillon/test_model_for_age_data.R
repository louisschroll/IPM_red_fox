#' HEADER ------------------------------------------------------------------------
#'
#' Script name:  ~/IPM_red_fox/2.code/brouillon/test_model_for_age_data.R
#' Author:       Louis Schroll
#' Email:        louis.schroll@ens-lyon.fr
#' Date:         2024-10-29
#'
#' Script description:
#'
#'
#' -------------------------------------------------------------------------------

cat("\014")              # clear the console
rm(list = ls())          # remove all variables of the work space

# Load packages
library(tidyverse)
library(jagsUI)       

# Load R functions
path_to_Rfunc <- "2.code/R_func"
sapply(paste0(path_to_Rfunc, "/", list.files(path_to_Rfunc)), source)

Tmax <- 10
N <- sim_pop_dyn(Tmax = Tmax)
C <- sim_age_data(N)

run_harvest_model <- function(C){
  age_max <- nrow(C)
  # Bundle data
  jags.data <- list(
    C = C,
    age_max = nrow(C),
    Tmax = ncol(C),
    h = 0.2
  )
  
  # Write model in jags
  cat("model{
      ################################
      #### PRIORS AND CONSTRAINTS ####
      ################################
      
      #h ~ dunif(0, 1)
      for(t in 1:Tmax) {
      lambda[1, t] ~ dgamma(10, 0.1)
      lambda[2, t] ~ dgamma(5, 0.1)
      lambda[3, t] ~ dgamma(2, 0.1)
      lambda[4, t] ~ dgamma(2, 0.1)
      lambda[5, t] ~ dgamma(2, 0.1)
        for(a in 1:age_max){
           
          N[a, t] ~ dpois(lambda[a, t])
          # h[a, t] ~ dunif(0, 1)
        }
      }
      
  
      ######################################
      #### AGE-AT-HARVEST MODULE ####
      ######################################
      
      ### Parameters:
      # N = number of individuals in a given age class at a given time
      # h = age- and time-dependent probability of dying from winter hunting
      
      ### Data:
      # C = age-at-harvest matrix

      ### Likelihood
      
      for(t in 1:Tmax){
        for(a in 1:age_max){
          C[a, t] ~ dbin(h, N[a, t]) #h[a, t]
        }
      }
      
      for (t in 1:Tmax){
        Ntot[t] <- sum(N[1:age_max, t])
      }
}", fill = TRUE, file = "model_age.txt")
  
  
  # Inits
  inits <- function() {
    list(#h = runif(1, 0, 1),
      lambda = matrix(rgamma(n = age_max * (Tmax+1), shape = 2, rate = 0.1),
                      nrow = age_max))
  }
  
  # Parameters to save
  params <- c("h", "Ntot", "N")
  
  # MCMC settings
  ni <- 12000
  nb <- 2000
  nt <- 1
  nc <- 3
  
  
  # Run the model in JAGS
  jags.output <- jags(
    jags.data,
    inits,
    params,
    "model_age.txt",
    n.thin = nt,
    n.chains = nc,
    n.burnin = nb,
    n.iter = ni,
    parallel = TRUE
  )
  return(jags.output)
}

jags.output <- run_harvest_model(C = C)
print(jags.output, 2)

# Check convergence
MCMCvis::MCMCtrace(jags.output, params = c("Ntot"), pdf = F)

dd <- tibble(N_true = colSums(N),
       N_estim = jags.output$mean$Ntot) %>% 
  mutate(residuals = N_true - N_estim,
         year = 1:(Tmax+1))

dd

sum(abs(dd$residuals))

dd %>% pivot_longer(-c(year, residuals)) %>%
  ggplot(aes(x = year, y = value, colour = name)) +
  geom_point() +
  geom_line() +
  labs(x = "Year", 
       y = "Residual (N_true - N_estim)", 
       title = "True and estimated value of N by year") +
  theme_minimal() +
  ylim(c(0, 320))

ggplot(dd, aes(x = year, y = residuals)) +
  geom_line() +
  geom_point() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "red") +
  labs(x = "Year", 
       y = "Residual (N_true - N_estim)", 
       title = "Residuals Over Time") +
  theme_minimal()

ggplot(dd, aes(x = residuals)) +
  geom_histogram(binwidth = 5,
                 fill = "skyblue",
                 color = "black") +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red") +
  labs(x = "Residual (N_true - N_estim)", 
       y = "Frequency", 
       title = "Histogram of Residuals") +
  theme_minimal()

ggplot(dd, aes(x = N_estim, y = residuals)) +
  geom_point() +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             color = "red") +
  labs(x = "Estimated Ntot (N_estim)", 
       y = "Residual (N_true - N_estim)", 
       title = "Residuals vs Estimated Ntot") +
  theme_minimal()

dd

sum(abs(dd$residuals))


N_age_estim <- jags.output$mean$N
N
n_simu <- 100

sum_diff_matrix <- N-N_age_estim
for (i in 1:n_simu){
  N <- sim_pop_dyn(Tmax = Tmax)
  C <- sim_age_data(N)
  jags.output <- run_harvest_model(C = C)
  
  N_age_estim <- jags.output$mean$N
  diff_matrix <- N-N_age_estim
  sum_diff_matrix <- sum_diff_matrix + diff_matrix
}


# Convert the matrix to a data frame
N_U_df <- as.data.frame(sum_diff_matrix)
colnames(N_U_df) <- paste0("Y", 1:ncol(N_U_df)) 
N_U_df$row <- 1:nrow(N_U_df)                    

# Reshape to long format
N_U_long <- pivot_longer(N_U_df, cols = starts_with("Y"), 
                         names_to = "Year", values_to = "Difference")

ggplot(N_U_long, aes(x = Year, y = factor(row), fill = Difference)) +
  geom_tile(color = "white") +                       # Add gridlines
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  labs(x = "Year", y = "Age Group", fill = "N - U") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

