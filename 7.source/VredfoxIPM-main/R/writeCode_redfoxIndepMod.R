#' Write NIMBLE code for red fox independent models
#'
#' @param indLikelihood.genData logical. If TRUE, writes a model that includes 
#' an individual-level likelihood for genetic data which treats and individual's
#' true immigration status as the outcome of a Bernouilli trial with probability
#' corresponding to a "probability of belonging" derived from Geneclass 2 output.
#' If FALSE (default), writes a model that includes a group-level likelihood for
#' genetic data which treats the number of immigrants (identified using 
#' Geneclass 2 and a pre-defined p threshold) as the outcome of a Poisson trial
#' with lambda equal to the number of resident individuals times immigration rate.
#' For models that do not make use of genetic data, this choice is irrelevant.
#' 
#' @return an R call object specifying the model structure for the red fox IPM. 
#' @export
#'
#' @examples

writeCode_redfoxIndepMod <- function(indLikelihood.genData = FALSE){
  
  ## Check for incompatible toggles
  if(!imm.asRate & fitCov.immR){
    stop("Incompatible model settings. Rodent covariate effect on immigration can only be fit (fitCov.immR = TRUE) if immigration is estimated as a rate (imm.asRate = TRUE).")
  }
  
  if(fitCov.mO & rCov.idx){
    stop("Incompatible model settings. Rodent effects on natural mortality (fitCov.mO = TRUE) are only implemented with continuous covariates (rCov.idx = FALSE).")
  }
  
  if(fitCov.mO & !mO.varT){
    warning("Attempting to fit a model containing environmental covariates but no random year variation for natural mortality. This is not recommended as it may result in inflated effect size/precision due to due to pseudo-replication.")
  }
  
  ## Write model code (individua-level likelihood for genetic data)
  if(indLikelihood.genData){
    stop("Code for independent modelling including individual-genetic likelihood is not available. For running independent models, please select indLikelihood.genData = FALSE.")
  }
  
  #=============================================================================
  #=============================================================================
  #=============================================================================
  
  
  ## Write model code (group-level likelihood for genetic data)
  if(!indLikelihood.genData){
    redfox.code <- nimbleCode({
      
      ###############################
      #### PLACENTAL SCAR MODULE ####
      ###############################
      
      ### Parameters:
      # rho = expected number of placental scars (fetuses)
      # Psi = pregnancy rate
      
      ## Data:
      # P1 = individual placental scar counts
      # P1_age = individual ages associated with P1
      # P1_year = year associated with P1
      
      # P2 = individual presence/absence of placental scars
      # P2_age = individual ages associated with P2
      # P2_year = year associated with P2
      
      
      ### Likelihood (litter size)
      
      for(x in 1:X1){
        P1[x] ~ dpois(rho[P1_age[x], P1_year[x]])
      }
      
      ### Likelihood (pregnancy rate)
      
      for(x in 1:X2){
        P2[x] ~ dbern(Psi[P2_age[x], P2_year[x]])
      }
      
      #===============================================================================================
      
      
      ###########################################
      #### GENETIC IMMIGRATION STATUS MODULE ####
      ###########################################
      
      ### Parameters:
      # Mu.immR = average immigration rate
      # immR = annual immigration rates
      
      ## Data:
      # pImm = individual probability of being an immigrant
      
      
      ### Likelihood (immigration status of sampled individuals)
      if(imm.asRate & useData.gen){
        
        if(poolYrs.genData){
          
          genObs_Imm ~ dpois(genObs_Res*Mu.immR)
          
        }else{
          
          # Within study period 
          for(t in 1:Tmax_Gen){
            genObs_Imm[t] ~ dpois(genObs_Res[t]*immR[t])
          }
          
          # Outside study period
          for(t in 1:Tmax_Gen_pre){
            genObs_Imm_pre[t] ~ dpois(genObs_Res_pre[t]*immR_pre[t])
          }
        }
      }
      
      #===============================================================================================
      
      
      ################################
      #### PUP OBSERVATION MODULE ####
      ################################
      
      ## Likelihood (Number of pups)
      if(useData.pup){
        for(x in 1:X3){
          NoPups[x] ~ dpois(meanLS[NoPups_year[x]])
        }
        
        meanLS[1] <- 0
        for(t in 2:Tmax){
          meanLS[t] ~ dunif(0, 20)
        }
      }
      
      
      #===============================================================================================
      
      
      ################################
      #### PRIORS AND CONSTRAINTS ####
      ################################
      
      ## Survival and mortality
      
      for(t in 1:Tmax){ 
        
        # Summer harvest mortality hazard rate
        log(mHs[1:Amax, t]) <- log(Mu.mHs[1:Amax]) + epsilon.mHs[t]
        
        # Winter harvest mortality hazard rate
        if(fitCov.mH){
          log(mH[1:Amax, t]) <- log(Mu.mH[1:Amax]) + betaHE.mH*HarvestEffort[t] + epsilon.mH[t]
        }else{
          log(mH[1:Amax, t]) <- log(Mu.mH[1:Amax]) + epsilon.mH[t]
        }
        
        # Other (natural) mortality hazard rate
        if(fitCov.mO){
          log(mO[1:Amax, t]) <- log(Mu.mO[1:Amax]) + betaRd.mO*Reindeer[t] + betaR.mO*RodentAbundance[t+1] + betaRxRd.mO*Reindeer[t]*RodentAbundance[t+1] + epsilon.mO[t]
        }else{
          log(mO[1:Amax, t]) <- log(Mu.mO[1:Amax]) + epsilon.mO[t]
        }
        
        # Survival probability
        S[1:Amax, t] <- exp(-(mH[1:Amax, t] + mO[1:Amax,t]))
        
        # Proportion winter harvest mortality
        alpha[1:Amax, t] <- mH[1:Amax, t]/(mH[1:Amax, t] + mO[1:Amax, t])
        
        # Winter harvest rate
        h[1:Amax, t] <- (1-S[1:Amax, t])*alpha[1:Amax, t]
        
      }
      
      # Median harvest mortality hazard rates
      
      # Age-dependent
      for(a in 1:Amax){
        Mu.mH[a] ~ dunif(0, 5)
        Mu.mHs[a] ~ dunif(0, 5)
      }
      
      # Age-independent   
      #Mu.mH.all ~ dunif(0, 5) 
      #Mu.mH[1:Amax] <- Mu.mH.all
      
      # Median other (natural) cause mortality hazard rates
      #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE / HOENIG MODEL CALCULATION
      
      if(HoenigPrior){
        # Using prior distributions calculated with Hoenig model
        Mu.mO.ad ~ dlnorm(mnat.logmean, sdlog = mnat.logsd)
        Mu.mO[2:5] <- Mu.mO.ad
        Mu.mO[1] <- Mu.mO.ad*JuvAdRatio #* NOTE: Can be provided as constant or distribution
        
        #JuvAdRatio <- exp(JAratio.logmean)
        JuvAdRatio ~ dlnorm(ratioJA.logmean, sdlog = ratioJA.logsd)
        
      }else{
        # Using literature values on age-specific survival
        for(a in 1:Amax){
          Mu.mO[a] <- -log(Mu.Snat[a])
          Mu.Snat[a] ~ T(dnorm(Snat.mean[a], sd = Snat.sd[a]), 0, 1)   
        }
      }
      
      
      ## Covariate effects
      if(fitCov.mH){
        betaHE.mH ~ dunif(0, 5) # Effect of harvest effort on mH
      }
      
      if(fitCov.mO){
        betaRd.mO ~ dunif(-5, 0) # Effect of reindeer carcasses on mO
        betaR.mO ~ dunif(-5, 0) # Effect of rodent abundance on mO
        betaRxRd.mO ~ dunif(-5, 5) # Interactive effect of reindeer x rodent on mO
      }
      
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Pregnancy rate
      
      for(t in 1:(Tmax+1)){
        Psi[1, t] <- 0
        
        if(fitCov.Psi){
          if(rCov.idx){
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi[RodentIndex[t]] + epsilon.Psi[t] # Reindeer.rodent interaction not (yet) written in
          }else{
            logit(Psi[2:Amax,t]) <- logit(Mu.Psi[2:Amax]) + betaR.Psi*RodentAbundance[t] + epsilon.Psi[t]
          }
        }else{
          logit(Psi[2:Amax, t]) <- logit(Mu.Psi[2:Amax]) + epsilon.Psi[t]
        }
      }
      
      Mu.Psi[1] <- 0
      for(a in 2:Amax){	
        Mu.Psi[a] ~ dunif(0, 1)
      }
      
      if(fitCov.Psi){
        if(rCov.idx){
          betaR.Psi[1] <- 0 # --> Lowest level corresponds to intercept
          for(x in 2:nLevels.rCov){
            betaR.Psi[x] ~ dunif(-5, 5)
          }
        }else{
          betaR.Psi ~ dunif(-5, 5)
        }
      }
      
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Litter size
      
      
      for(t in 1:(Tmax+1)){
        rho[1, t] <- 0
        
        if(fitCov.rho){
          if(rCov.idx){
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho[RodentIndex[t]] + epsilon.rho[t]
          }else{
            log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + betaR.rho*RodentAbundance[t] + epsilon.rho[t]
          }
        }else{
          log(rho[2:Amax, t]) <- log(Mu.rho[2:Amax]) + epsilon.rho[t]
        }
      }
      
      Mu.rho[1] <- 0
      for(a in 2:Amax){
        Mu.rho[a] ~ dunif(0, maxPups) # Baseline number of pups 
      }
      
      if(fitCov.rho){
        if(rCov.idx){
          betaR.rho[1] <- 0 # --> Lowest level corresponds to intercept
          for(x in 2:nLevels.rCov){
            betaR.rho[x] ~ dunif(-5, 5)
          }
        }else{
          betaR.rho ~ dunif(-5, 5)
        }
      }
      
      #---------------------------------------------------------------------------------------------  
      
      
      ## Denning survival
      
      for(t in 1:(Tmax+1)){ 
        S0[t] <- Mu.S0
        #S0[t] <- exp(-m0[t])
        #log(m0[t]) <- log(-log(Mu.S0)) + epsilon.m0[t]
      }
      
      if(useInfPrior.S0){
        Mu.S0 ~ T(dnorm(S0.mean, sd = S0.sd), 0, 1)
      }else{
        Mu.S0 ~ dunif(0, 1)
      }
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Immigration (within the study period)
      
      if(imm.asRate){
        
        if(fitCov.immR){
          if(rCov.idx){
            for(t in 1:(Tmax+1)){
              log(immR[t]) <- log(Mu.immR) + betaR.immR[RodentIndex2[t]] + epsilon.immR[t]
            }
          }else{
            log(immR[1:(Tmax+1)]) <- log(Mu.immR) + betaR.immR*RodentAbundance2[1:(Tmax+1)] + epsilon.immR[1:(Tmax+1)]
          }
        }else{
          log(immR[1:(Tmax+1)]) <- log(Mu.immR) + epsilon.immR[1:(Tmax+1)]
        }
        
        
        for(t in 1:Tmax){ 
          Imm[t] ~ dpois(survN1[t]*immR[t])
        }
        
        Mu.immR ~ dunif(0, 10)
        
        
      }else{
        
        ## Lognormal prior for immigrant numbers
        for(t in 2:Tmax){
          Imm[t] <- round(ImmExp[t])
          ImmExp[t] ~ dlnorm(meanlog = log(Mu.Imm), sdlog = logsigma.Imm) 
        }
        
        Mu.Imm ~ dunif(1, uLim.Imm)
        logsigma.Imm ~ dunif(0, 10)
        
        ## Derivation of immigration rates
        immR[1] <- 0
        for(t in 2:Tmax){
          immR[t] <- Imm[t] / survN1[t]
        }
        
      }
      
      ## Prior for rodent effect
      if(fitCov.immR){
        if(rCov.idx){
          betaR.immR[1] <- 0 # --> Lowest level corresponds to intercept
          for(x in 2:nLevels.rCov){
            betaR.immR[x] ~ dunif(-5, 5)
          }
        }else{
          betaR.immR ~ dunif(-5, 5)
        }
      }
      
      
      ## Immigration outside the study period
      if(imm.asRate & useData.gen & !poolYrs.genData){
        
        for(t in 1:Tmax_Gen_pre){
          
          if(fitCov.immR){
            if(rCov.idx){
              immR_pre[t] ~ dlnorm(log(Mu.immR) + betaR.immR[RodentIndex2_pre[t]], sdlog = sigma.immR)
            }else{
              immR_pre[t] ~ dlnorm(log(Mu.immR) + betaR.immR*RodentAbundance2_pre[t], sdlog = sigma.immR)
            }
          }else{
            immR_pre[t] ~ dlnorm(log(Mu.immR), sdlog = sigma.immR)
          }
        }
      }
      
      #---------------------------------------------------------------------------------------------
      
      
      ## Random year variation
      for(t in 1:Tmax){  
        epsilon.mHs[t] ~ dnorm(0, sd = sigma.mHs)
        epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
        epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
      }
      
      for(t in 1:(Tmax+1)){
        epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
        epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 
        # epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
      }
      
      sigma.mHs ~ dunif(0, 5)
      sigma.mH ~ dunif(0, 5)
      sigma.Psi ~ dunif(0, 5)
      sigma.rho ~ dunif(0, 5)

      if(mO.varT){
        sigma.mO ~ dunif(0, 5)
      }else{
        sigma.mO <- 0
      }
      
      if(imm.asRate){
        for(t in 1:(Tmax+1)){
          epsilon.immR[t] ~ dnorm(0, sd = sigma.immR)
        }
        sigma.immR ~ dunif(0, 10)
      }
      
      #===============================================================================================
      
      
      
      ##############################
      #### COVARIATE IMPUTATION ####
      ##############################
      
      ## Missing covariate value(s) in number of successful hunters
      if(fitCov.mH){
        for(t in 1:Tmax){
          HarvestEffort[t] ~ dnorm(0, sd = 1)
        }
      }
      
      ## Missing covariate value(s) in rodent abundance
      if(rCov.idx){
        
        for(t in 1:Tmax+1){
          RodentIndex[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
          RodentIndex2[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
        }
        
        if(imm.asRate & fitCov.immR & !poolYrs.genData){
          for(t in 1:Tmax_Gen_pre){
            RodentIndex2_pre[t] ~ dcat(DU.prior.rCov[1:nLevels.rCov]) 
          }
        }
        
        DU.prior.rCov[1:nLevels.rCov] <- 1/nLevels.rCov
        
      }else{
        
        for(t in 1:Tmax+1){
          RodentAbundance[t] ~ dnorm(0, sd = 1)
          RodentAbundance2[t] ~ dnorm(0, sd = 1)
        }
        
        if(imm.asRate & fitCov.immR & !poolYrs.genData){
          for(t in 1:Tmax_Gen_pre){
            RodentAbundance2_pre[t] ~ dnorm(0, sd = 1) 
          }
        }
      }
      
      ## Missing covariate values in reindeer information
      if(fitCov.mO){
        for(t in 1:Tmax+1){
          Reindeer[t] ~ dnorm(0, sd = 1)
        }
      }
      
    })
  }
  
  
  return(redfox.code)
}

