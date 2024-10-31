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