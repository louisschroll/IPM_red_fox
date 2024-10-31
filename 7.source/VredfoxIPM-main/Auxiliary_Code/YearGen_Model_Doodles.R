calculateImmR <- nimbleFunction(
  run = function(ImmData = double(1), yearIdx = double(1),
                 Tmax = double(), skip_t1 = logical()){
    
    # Set up vector of time-dependent immigration rates
    immR <- rep(0, Tmax)
    
    if(skip_t1){
      Tstart <- 2
    }else{
      Tstart <- 1
    }
    
    for(t in Tstart:Tmax){
      
      # Extract relevant year indices
      idx.select <- which(yearIdx == t)
      
      # Calculate immigration rate for year t
      ImmData.select <- ImmData[idx.select]
      
      if(sum(ImmData.select) == length(ImmData.select)){
        immR[t] <- (sum(ImmData.select) - 1) / (length(ImmData.select) - sum(ImmData.select))
      }else{
        immR[t] <- sum(ImmData.select) / (length(ImmData.select) - sum(ImmData.select))
      }
    }
    
    return(immR)
    
    returnType(double(1))
  })



calculateLogSD <- nimbleFunction(
  run = function(immR = double(1), replace0 = double()){
    
    # Extract relevant time indices (immR = 0)
    idx.select <- which(immR == 0)
    
    # Replace 0 immigration rates with "replaceO" value
    immR[idx.select] <- replace0
    
    # Calculate log standard deviation
    immR_sdlog <- sd(log(immR))
    
    return(immR_sdlog)
    
    returnType(double())
  })


mod_test <- nimbleCode({
  
  ## Likelihoods
  for(x in 1:Xgen){
    ImmData[x] ~ dbern(pImm[x])
  }
  
  for(x in 1:Xgen_pre){
    ImmData_pre[x] ~ dbern(pImm_pre[x])
  }
  
  ## Calculation
  immR[1:Tmax_Gen] <- calculateImmR(ImmData = ImmData[1:Xgen], 
                                    yearIdx = pImm_yrs[1:Xgen],
                                    Tmax = Tmax_Gen, skip_t1 = FALSE)
  
  immR_pre[1:Tmax_Gen_pre] <- calculateImmR(ImmData = ImmData_pre[1:Xgen_pre], 
                                            yearIdx = pImm_yrs_pre[1:Xgen_pre],
                                            Tmax = Tmax_Gen_pre, skip_t1 = FALSE)
  
  ## Downstream calculation
  log(Mu.immR) <- log(mean(c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]))) 
  sigma.immR <- calculateLogSD(immR = c(immR[1:Tmax_Gen], immR_pre[1:Tmax_Gen_pre]),
                               replace0 = 0.01)

  
  ## Projection
  for(t in (Tmax_Gen+1):(Tmax+1)){
    immR[t] ~ dlnorm(meanlog = log(Mu.immR), sdlog = sigma.immR)
  }
  
})


test_run <- nimbleMCMC(code = mod_test, 
                       constants = input.data$nim.constants,
                       data = input.data$nim.data,
                       inits = model.setup$initVals[[1]], 
                       niter = 100,
                       nchains = 1,
                       monitors = c("immR", "immR_pre", "Mu.immR", "sigma.immR"),
                       samplesAsCodaMCMC = TRUE)

head(test_run)  
plot(test_run, ask = T)


