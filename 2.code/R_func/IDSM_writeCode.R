

#' Write integrated distance sampling model code
#'
#' @return an R call object specifying the model structure for integrated
#' distance sampling model.
#' @export
#'
#' @examples

IDSM_writeCode <- function() {
  IDSM.code <- nimble::nimbleCode({
    # n_areas = number of areas
    # n_sites = number of sites in area x
    # n_age_class = number of age classes
    # n_years = number of years
    # N_sumR_obs = number of data points in juvenile:adult ratio counts
    
    # N_exp[a, j, t] = Number of age class a individuals in site j of area x in year t
    # Density[a, j, t] = Density of age class a individuals in site j of area x in year t
    # L[j, t] = length of transect line in site j of area x in year t
    # W = truncation distance for line transect surveys
    
    # S[t] = annual survival from year t to t+1 in area x
    # R_year[t] = recruitment rate in year t in area x
    # p[t] = average distance sampling detection rate in area x in year t
    # sigma[t] = average distance sampling detection decay rate in area x in year t
    
    ####################
    # POPULATION MODEL #
    ####################
    
    #-----------------------------------------#
    # Initial population size/density (t = 1) #
    #-----------------------------------------#
    for (j in 1:n_sites) {
      # Age class 0 (index = 1): reproduction
      #Density[1, j, 1] <- sum(Density[2:n_age_class, j, 1]) / 2 * mean.recruitment #R_year[1]
      
      # Age classes 1 to 4 (indeces = 2, 3, 4, 5)
      for(a in 1:n_age_class){
        Density[a, j, 1] ~ dunif(0, 1) #dbeta(beta_param[a, 1], beta_param[a, 2])
      }
      
      ## Adult and juvenile numbers
      N_exp[1:n_age_class, j, 1] <- Density[1:n_age_class, j, 1] * L[j, 1] * W * 2
    }
    
    #-------------------------------#
    # Population dynamics for t > 1 #
    #-------------------------------#
    for (j in 1:n_sites) {
      for (t in 2:n_years) {
        # Age class 0 (index = 1): reproduction
        Density[1, j, t] <- sum(Density[2:n_age_class, j, t]) / 2 * mean.recruitment  #R_year[t]
        
        # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors 
        for(a in 2:(n_age_class-1)){
          Density[a, j, t] <- Density[a - 1, j, t - 1] * mean.survival  #S[t - 1]
        }		
        
        # Age class 4+ (index = n_age_class = 5): age class 4 and 5+ survivors
        Density[n_age_class, j, t] <- (Density[n_age_class-1, j, t-1] + Density[n_age_class, j, t-1]) * mean.survival  #S[t-1]
        
        ## Adult and juvenile numbers
        N_exp[1:n_age_class, j, t] <- Density[1:n_age_class, j, t] * L[j, t] * W * 2
      }
    }
    
    #--------------------#
    # Derived parameters #
    #--------------------#
    
    ## Area- and year-specific total densities
    for (t in 1:n_years) {
      N_tot_exp[t] <- sum(N_exp[1:n_age_class, 1:n_sites, t])
    }
    
    ## Area-, year-, and age-class specific density (for monitoring)
    for (a in 1:n_age_class) {
      for (t in 1:n_years) {
        meanDens[a, t] <- mean(Density[a, 1:n_sites, t])
        N_tot_gic_age[a, t] <- round(meanDens[a, t] * size_hunting_area)
      } # t
    } # a
    
    for(t in 1:n_years){
      N_tot_gic[t] <- sum(N_tot_gic_age[1:n_age_class, t])
    }
    
    ############################
    # Distance sampling module #
    ############################
    
    ## Age-specific line transect counts
    # DS_count[j, t] = number of individuals detected in site j in year t
    for (j in 1:n_sites) {
      for (t in 1:n_years) {
        DS_count[j, t] ~ dbin(p, sum(N_exp[1:n_age_class, j, t])) #dpois(p * sum(N_exp[1:n_age_class, j, t]))
      }
    }
    
    ## Line transect observation distances (likelihood using nimbleDistance::dHN)
    # N_obs = number of observations in detection distance data in area x
    # d[i] = i'th entry in detection distance data for area x
    for (i in 1:n_obs_dist) {
      d[i] ~ dHN(sigma = sigma, # sigma[year_obs[i]],
                 Xmax = W,
                 point = 0)
    }
    
    # PARAMETER MODELS/CONSTRAINTS
    ## Distance sampling detection parameters
    # Detection decay
    log(sigma) <- 0.15 #log.mean.sigma
    sigma2 <- sigma * sigma
    
    # Effective strip width
    esw <- sqrt(pi * sigma2 / 2)
    
    # Average detection rate
    p <- min(esw, W) / W
    
    # ## Annual recruitment rates
    # R_year[1:n_years] <- mean.recruitment
    # 
    # ## Annual survival probabilities
    # logit(Mu.S) <- mean.survival
    # S[1:(n_years - 1)] <- Mu.S
    
    ###########
    # PRIORS  #
    ###########
    mean.recruitment <- 1.578 # dunif(0, 20) # Recruitment
    mean.survival <- 0.54 # ~ dunif(0, 1)      # Survival
    log.mean.sigma ~ dunif(-10, 1) # Detection
    
    ##################################
    # Module for age-at-harvest data #
    ##################################
    
    ## Parameters:
    # N = number of individuals in a given age class at a given time
    # h = age- and time-dependent probability of dying from hunting
    
    ## Data:
    # C = age-at-harvest matrix
    
    ## Priors
    harvest_rate <- 0.2 # ~ dunif(0, 1)

    ## Likelihood

    for(t in 1:n_years){
      for(a in 1:n_age_class){
        C[a, t] ~ dbin(harvest_rate, N_tot_gic_age[a, t])
      }
    }
  })
  return(IDSM.code)
}