

#' Write integrated distance sampling model code
#'
#' @param survVarT logical. If TRUE, writes code for a model including random
#' year variation in survival probability. If FALSE, assumed constant survival
#' probability across time.
#' @param telemetryData logical. If TRUE, uses information from telemetry data
#' from Lierne. If FALSE, only line transect data is used.
#' @return an R call object specifying the model structure for integrated
#' distance sampling model.
#' @export
#'
#' @examples

writeIDSMCode <- function(survVarT) {
  IDSM.code <- nimble::nimbleCode({
    # n_areas = number of areas
    # n_sites[x] = number of sites in area x
    # age_max = number of age classes
    # n_years = number of years
    # N_sumR_obs[x] = number of data points in juvenile:adult ratio counts
    
    # N_exp[x, a, j, t] = Number of age class a individuals in site j of area x in year t
    # Density[x, a, j, t] = Density of age class a individuals in site j of area x in year t
    # L[x, j, t] = length of transect line in site j of area x in year t
    # W = truncation distance for line transect surveys
    
    # Mu.D1[x] = average initial density in area x
    
    # S[x, t] = annual survival from year t to t+1 in area x
    # R_year[x, t] = recruitment rate in year t in area x
    # p[x, t] = average distance sampling detection rate in area x in year t
    # sigma[x, t] = average distance sampling detection decay rate in area x in year t
    
    # eps.D1[x, j] = random site effect on initial density area x (site j)

    ####################
    # POPULATION MODEL #
    ####################
    for (x in 1:n_areas) {
      #-----------------------------------------#
      # Initial population size/density (t = 1) #
      #-----------------------------------------#
      for (j in 1:n_sites[x]) {
        ## Adult densities
        Density[x, 2, j, 1] <- exp(log(Mu.D1[x]) + eps.D1[x, j])
        
        ## Juvenile densities
        Density[x, 1, j, 1] <- Density[x, 2, j, 1] * R_year[x, 1]
        
        ## Adult and juvenile numbers
        N_exp[x, 1:age_max, j, 1] <- Density[x, 1:age_max, j, 1] * L[x, j, 1] * W * 2
      }
      
      #-------------------------------#
      # Population dynamics for t > 1 #
      #-------------------------------#
      for (j in 1:n_sites[x]) {
        for (t in 2:n_years) {
          ## Adult densities
          Density[x, 2, j, t] <- sum(Density[x, 1:age_max, j, t - 1]) * S[x, t - 1]
          
          ## Juvenile densities
          Density[x, 1, j, t] <- (Density[x, 2, j, t] / 2) * R_year[x, t]
          
          ## Adult and juvenile numbers
          N_exp[x, 1:age_max, j, t] <- Density[x, 1:age_max, j, t] * L[x, j, t] * W * 2
        }
      }
      
      #--------------------#
      # Derived parameters #
      #--------------------#
      
      ## Area- and year-specific total densities
      for (t in 1:n_years) {
        N_tot_exp[x, t] <- sum(N_exp[x, 1, 1:n_sites[x], t] + N_exp[x, 2, 1:n_sites[x], t])
      }
      
      ## Area-, year-, and age-class specific density (for monitoring)
      for (a in 1:age_max) {
        for (t in 1:n_years) {
          meanDens[x, a, t] <- mean(Density[x, a, 1:n_sites[x], t])
        } # t
      } # a
    } # x
    
    
    ############################
    # Distance sampling module #
    ############################
    for (x in 1:n_areas) {
      ## Age-specific line transect counts
      # DS_count[x, j, t] = number of individuals detected in site j of area x in year t
      for (j in 1:n_sites[x]) {
        for (t in 1:n_years) {
          DS_count[x, j, t] ~ dpois(p[x, t] * sum(N_exp[x, 1:age_max, j, t]))
        }
      }
      
      ## Line transect observation distances (likelihood using nimbleDistance::dHN)
      # N_obs[x] = number of observations in detection distance data in area x
      # d[x, i] = i'th entry in detection distance data for area x
      for (i in 1:N_obs[x]) {
        d[x, i] ~ dHN(sigma = sigma[x, Year_obs[x, i]],
                      Xmax = W,
                      point = 0)
      }
    } # x
    
    ################################
    # PARAMETER MODELS/CONSTRAINTS #
    ################################
    
    for (x in 1:n_areas) {
      ## Distance sampling detection parameters
      
      for (t in 1:n_years) {
        # Detection decay
        log(sigma[x, t]) <- log.mean.sigma.area[x] + epsT.dd[t] + epsR.dd[x, t]
        sigma2[x, t] <- sigma[x, t] * sigma[x, t]
        
        # Effective strip width
        esw[x, t] <- sqrt(pi * sigma2[x, t] / 2)
        
        # Average detection rate
        p[x, t] <- min(esw[x, t], W) / W
      }
      
      ## Annual recruitment rates
      R_year[x, 1:n_years] <- exp(log(Mu.R[x]) + epsT.R[1:n_years] + epsR.R[x, 1:n_years])
      
      ## Annual survival probabilities
      logit(Mu.S[x]) <- mu.S[x]
      
      if (survVarT) {
        logit(S[x, 1:(n_years - 1)]) <- logit(Mu.S[x]) + epsT.S[1:(n_years - 1)] + epsR.S[x, 1:(n_years - 1)]
      } else{
        S[x, 1:(n_years - 1)] <- Mu.S[x]
      }
    } # x
    
    ###########
    # PRIORS  #
    ###########
    
    #-----------------------#
    # Intercepts / averages #
    #-----------------------#
    
    h.Mu.R  ~ dunif(0, 20) # Recruitment
    h.Mu.S ~ dunif(0, 1) # Survival
    log.mean.sigma ~ dunif(-10, 100) # Detection
    
    for (x in 1:n_areas) {
      ## Initial density
      Mu.D1[x] ~ dunif(0, 10)
      
      ## Recruitment
      epsA.R[x]  ~ dnorm(0, sd = h.sigma.R)
      log(Mu.R[x]) <- log(h.Mu.R) + epsA.R[x]
      
      ## Survival
      epsA.S[x]  ~ dnorm(0, sd = h.sigma.S)
      mu.S[x] <- logit(h.Mu.S) + epsA.S[x]
      
      ## Detection
      epsA.dd[x] ~ dnorm(0, sd = h.sigma.dd)
      log.mean.sigma.area[x] <- log.mean.sigma + epsA.dd[x]
    }
    
    
    #----------------#
    # Random effects #
    #----------------#
    
    ## Standard deviations
    
    # Recruitment
    h.sigma.R ~ dunif(0, 5)
    sigmaT.R ~ dunif(0, 5)
    sigmaR.R ~ dunif(0, 5)
    
    # Survival
    h.sigma.S ~ dunif(0, 5)
    Mu.S1 ~ dunif(0, 1)
    
    if (survVarT) {
      sigmaT.S ~ dunif(0, 5)
      sigmaR.S ~ dunif(0, 5)
      eps.S1.prop ~ dunif(0, 1)
    }
    
    # Detection
    h.sigma.dd ~ dunif(0, 5)
    sigmaT.dd ~ dunif(0, 20)
    sigmaR.dd ~ dunif(0, 20)
    
    # Initial density
    for (x in 1:n_areas) {
      sigma.D[x] ~ dunif(0, 20)
    }
    
    ## Random effect levels
    # Shared year variation
    for (t in 1:n_years) {
      epsT.R[t] ~ dnorm(0, sd = sigmaT.R) # Recruitment
      epsT.dd[t] ~ dnorm(0, sd = sigmaT.dd) # Detection
    }
    
    for (t in 1:(n_years - 1)) {
      if (survVarT) {
        epsT.S[t] ~ dnorm(0, sd = sigmaT.S) # Survival
      }
    }
    
    # Residual variation
    for (x in 1:n_areas) {
      for (t in 1:n_years) {
        epsR.R[x, t] ~ dnorm(0, sd = sigmaR.R)
        epsR.dd[x, t] ~ dnorm(0, sd = sigmaR.dd)
      }
    }
    
    for (x in 1:n_areas) {
      for (t in 1:(n_years - 1)) {
        if (survVarT) {
          epsR.S[x, t] ~ dnorm(0, sd = sigmaR.S)
        }
      }
    }
    
    # Site/transect variation
    for (x in 1:n_areas) {
      for (j in 1:n_sites[x]) {
        eps.D1[x, j] ~ dnorm(0, sd = sigma.D[x])
      }
    }
    
    #------------------#
    # Other parameters #
    #------------------#
    pi <- 3.141593
    
    ##################################
    # Module for age-at-harvest data #
    ##################################
    
    ## Parameters:
    # N = number of individuals in a given age class at a given time
    # h = age- and time-dependent probability of dying from hunting
    
    ## Data:
    # C = age-at-harvest matrix
    
    ## Priors
    h ~ dunif(0, 1)
    # for(t in 1:Tmax) {
    #   for(a in 1:age_max){
    #     h[a, t] ~ dunif(0, 1)
    #   }
    # }
    
    ## Likelihood
    for (x in 1:n_areas) {
      for(t in 1:Tmax){
        for(a in 1:age_max){
          N_area[x, a, t] <- meanDens[x, a, t] * size_area[x]
          C[x, a, t] ~ dbin(h, N_area[x, a, t]) #h[a, t]
        }
      }
    }

  })
  
  return(IDSM.code)
}