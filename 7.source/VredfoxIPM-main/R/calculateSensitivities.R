#' Calculate transient sensitivities and elasticities
#'
#' @param paramSamples a list of lists containing posterior samples for all vital rates and
#' population-level quantities.
#' @param Amax integer. Number of age classes. 
#' @param t_period integer vector. Optional argument for specifying the years
#' to use for calculating sensitivities. 
#' @return a list of lists containing posterior samples for transient 
#' sensitivities and elasticities for all vital rate parameters as well as
#' population structure (n) and population sizes per age class (N). 
#' 
#' @export
#'
#' @examples

calculateSensitivities <- function(paramSamples, Amax, t_period = NULL){
  
  
  #----------------------------------------#
  # CALCULATION OF TRANSIENT SENSITIVITIES #
  #----------------------------------------#
  
  ## Unpack objects from parameter list
  if(is.null(t_period)){
    for(i in 1:length(paramSamples$t_mean)){
      assign(names(paramSamples$t_mean)[i], 
             paramSamples$t_mean[[i]])
    } 
  }else{
    for(i in 1:length(paramSamples$t_mean)){
      
      paramName <- names(paramSamples$t)[i]
      t_offset <- ifelse(paramName %in% c("Psi", "rho", "S0", "m0"), 1, 0)
      
      focalParam <- paramSamples$t[[i]]
      
      if(length(dim(focalParam)) > 2){
        focalParam_mean <- apply(focalParam[, , t_period + t_offset], c(1, 2), mean)
      }else{
        focalParam_mean <- rowMeans(focalParam[, t_period + t_offset])
      }
      assign(paramName, focalParam_mean)
    } 
  }

  
  ## Set sample number
  nosamples <- length(lambda)
  
  ## Set up list of arrays for storing transient sensitivities
  sensList <- list(
    sens_S = matrix(NA, nrow = nosamples, ncol = Amax),
    sens_mH = matrix(NA, nrow = nosamples, ncol = Amax),
    sens_mO = matrix(NA, nrow = nosamples, ncol = Amax),
    
    sens_Psi = matrix(NA, nrow = nosamples, ncol = Amax),
    sens_rho = matrix(NA, nrow = nosamples, ncol = Amax),
    
    sens_S0 = rep(NA, nosamples),
    sens_m0 = rep(NA, nosamples),
    
    sens_Ss = matrix(NA, nrow = nosamples, ncol = Amax),
    sens_mHs = matrix(NA, nrow = nosamples, ncol = Amax),
    
    sens_immR = rep(NA, nosamples),
    
    sens_n = matrix(NA, nrow = nosamples, ncol = Amax),
    sens_N = matrix(NA, nrow = nosamples, ncol = Amax)
  )
  
  
  ## Calculate transient sensitivities for vital rates and population size/structure (evaluated at the temporal mean)
  
  for(i in 1:nosamples){
    for(a in 1:Amax){
      
      x <- ifelse(a < Amax, 1, 0)
      
      if(a == 1){
        
        sensList$sens_S[i, a] <- n[i, a]*Ss[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])*(1 + immR[i])
        sensList$sens_Ss[i, a] <- n[i, a]*S[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])*(1 + immR[i])
        
        sensList$sens_Psi[i, a] <- 0
        sensList$sens_rho[i, a] <- 0
        
        sensList$sens_n[i, a] <-   Ss[i, a]*S[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])*(1 + immR[i])
        
      }else{
        
        sensList$sens_S[i, a] <- n[i, a]*Ss[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])
        sensList$sens_Ss[i, a] <- n[i, a]*S[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])
        
        if(a < Amax){
          
          if(a == 2){
            sensList$sens_Psi[i, a] <- n[i, a-1]*Ss[i, a-1]*S[i, a-1]*0.5*rho[i, a]*S0[i]*(1 + immR[i])
            sensList$sens_rho[i, a] <- n[i, a-1]*Ss[i, a-1]*S[i, a-1]*0.5*Psi[i, a]*S0[i]*(1 + immR[i])
          }else{
            sensList$sens_Psi[i, a] <- n[i, a-1]*Ss[i, a-1]*S[i, a-1]*0.5*rho[i, a]*S0[i]
            sensList$sens_rho[i, a] <- n[i, a-1]*Ss[i, a-1]*S[i, a-1]*0.5*Psi[i, a]*S0[i]
          }
          
        }else{
          sensList$sens_Psi[i, a] <- (n[i, a-1]*Ss[i, a-1]*S[i, a-1] + n[i, a]*S[i, a]*Ss[i, a])*0.5*rho[i, a]*S0[i]
          sensList$sens_rho[i, a] <- (n[i, a-1]*Ss[i, a-1]*S[i, a-1] + n[i, a]*S[i, a]*Ss[i, a])*0.5*Psi[i, a]*S0[i]
        }
        
        sensList$sens_n[i, a] <-   Ss[i, a]*S[i, a]*(1 + Psi[i, a+x]*0.5*rho[i, a+x]*S0[i])
      }
      
      sensList$sens_mH[i, a] <- sensList$sens_mO[i, a] <- -exp(-(mH[i, a] + mO[i, a]))*sensList$sens_S[i, a]
      sensList$sens_mHs[i, a] <- -exp(-mHs[i, a])*sensList$sens_Ss[i, a]
      
      sensList$sens_N[i, a] <- (sensList$sens_n[i, a] - lambda[i]) / (sum(N[i, 1:Amax]))
    } 
    
    sensList$sens_S0[i] <- n[i, 1]*Ss[i, 1]*S[i, 1]*Psi[i, 2]*0.5*rho[i, 2]*(1 + immR[i]) +
                           sum(n[i, 2:(Amax-1)]*Ss[i, 2:(Amax-1)]*S[i, 2:(Amax-1)]*Psi[i, 3:Amax]*0.5*rho[i, 3:Amax]) +  
                           n[i, Amax]*Ss[i, Amax]*S[i, Amax]*Psi[i, Amax]*0.5*rho[i, Amax]
    
    sensList$sens_m0[i] <- -exp(-m0[i])*sensList$sens_S0[i]
    
    sensList$sens_immR[i] <- n[i, 1]*Ss[i, 1]*S[i, 1]*(1 + Psi[i, 2]*0.5*rho[i, 2]*S0[i])

  }
  
  
  ## Get posterior summaries for transient sensitivities
  postSum_sens <- data.frame()
  
  for(i in 1:length(sensList)){
    
    if(is.matrix(sensList[[i]])){
      
      # Extract quantiles for age-specific parameters
      quantiles <- apply(sensList[[i]], 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
      dimnames(quantiles)[[1]] <- c("lCI", "median", "uCI")
      
      data_temp <- cbind(data.frame(Parameter = names(sensList[i]), AgeClass = 1:Amax), t(quantiles))
      
      # Extract quantiles for sensitivities summed over age classes
      sum_temp <- cbind(data.frame(Parameter = names(sensList[i]), AgeClass = "summed"), 
                        t(quantile(rowSums(sensList[[i]]), probs = c(0.025, 0.5, 0.975))))
      colnames(sum_temp)[3:5] <- c("lCI", "median", "uCI")
      
      data_temp <- rbind(data_temp, sum_temp)
      
      # Merge into storage data frame
      postSum_sens <- rbind(postSum_sens, data_temp)
      
    }else{
      
      # Extract quantiles for age-specific parameters
      quantiles <- quantile(sensList[[i]], probs = c(0.025, 0.5, 0.975))
      names(quantiles) <- c("lCI", "median", "uCI")
      
      data_temp <- cbind(data.frame(Parameter = names(sensList[i]), AgeClass = NA), t(quantiles))
      
      # Merge into storage data frame
      postSum_sens <- rbind(postSum_sens, data_temp)
    }
  }
  
  
  #---------------------------------------#
  # CALCULATION OF TRANSIENT ELASTICITIES #
  #---------------------------------------#
  
  ## Calculate transient elasticities for vital rates and population size/structure (evaluated at the temporal mean)
  elasList <- list(
    elas_S = sensList$sens_S*(S/lambda),
    elas_mH = sensList$sens_mH*(mH/lambda),
    elas_mO = sensList$sens_mO*(mO/lambda),
    
    elas_Psi = sensList$sens_Psi*(Psi/lambda),
    elas_rho = sensList$sens_rho*(rho/lambda),
    
    elas_S0 = sensList$sens_S0*(S0/lambda),
    elas_m0 = sensList$sens_m0*(m0/lambda),
    
    elas_Ss = sensList$sens_Ss*(S/lambda),
    elas_mHs = sensList$sens_mHs*(mHs/lambda),
    
    elas_immR = sensList$sens_immR*(immR/lambda),
    
    elas_n = sensList$sens_n*(n/lambda),
    elas_N = sensList$sens_N*(N/lambda)
  )
  
  
  ## Get posterior summaries for transient elasticities
  postSum_elas <- data.frame()
  
  for(i in 1:length(elasList)){
    
    if(is.matrix(elasList[[i]])){
      
      # Extract quantiles for age-specific parameters
      quantiles <- apply(elasList[[i]], 2, stats::quantile, probs = c(0.025, 0.5, 0.975))
      dimnames(quantiles)[[1]] <- c("lCI", "median", "uCI")
      
      data_temp <- cbind(data.frame(Parameter = names(elasList[i]), AgeClass = 1:Amax), t(quantiles))
      
      # Extract quantiles for elasitivities summed over age classes
      sum_temp <- cbind(data.frame(Parameter = names(elasList[i]), AgeClass = "summed"), 
                        t(quantile(rowSums(elasList[[i]]), probs = c(0.025, 0.5, 0.975))))
      colnames(sum_temp)[3:5] <- c("lCI", "median", "uCI")
      
      data_temp <- rbind(data_temp, sum_temp)
      
      # Merge into storage data frame
      postSum_elas <- rbind(postSum_elas, data_temp)
      
    }else{
      
      # Extract quantiles for age-specific parameters
      quantiles <- quantile(elasList[[i]], probs = c(0.025, 0.5, 0.975))
      names(quantiles) <- c("lCI", "median", "uCI")
      
      data_temp <- cbind(data.frame(Parameter = names(elasList[i]), AgeClass = NA), t(quantiles))
      
      # Merge into storage data frame
      postSum_elas <- rbind(postSum_elas, data_temp)
    }
  }
  
  
  #-----------------------------------------------#
  # OPTIONAL: PLOT SENSITIVITIES AND ELASTICITIES #
  #-----------------------------------------------#
  
  
  #-----------------------------------------------------------------------------
  
  ## Save and return results
  sensResults <- list(sensitivity = list(samples = sensList, 
                                         summaries = postSum_sens),
                      elasticity = list(samples = elasList,
                                        summaries = postSum_elas))
  
  if(is.null(t_period)){
    saveRDS(sensResults, file = "RedFoxIPM_Sensitivities.rds")
  }
  
  return(sensResults)
  
}
