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