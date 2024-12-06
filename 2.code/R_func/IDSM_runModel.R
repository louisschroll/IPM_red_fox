


IDSM_runModel <- function(DS_data, 
                          harvest_data, 
                          dist_max, 
                          size_hunting_area,
                          testRun = FALSE,
                          mySeed = 42){
  input_data <- IDSM_prepareInputData(DS_data = DS_data,
                                      harvest_data = harvest_data,
                                      W = dist_max,
                                      size_hunting_area = size_hunting_area)
  
  IDSM_code <- IDSM_writeCode()
  
  model_setup <- IDSM_setupModel(nim.data = input_data$nim.data,
                                 nim.constants = input_data$nim.constants,
                                 niter = 10000, 
                                 nthin = 1, 
                                 nburn = 2000, 
                                 nchains = 3,
                                 testRun = testRun,
                                 initVals.seed = mySeed)
  
  IDSM_out <- nimbleMCMC(code = IDSM_code,
                         data = input_data$nim.data, 
                         constants = input_data$nim.constants,
                         inits = model_setup$initVals, 
                         monitors = model_setup$modelParams,
                         nchains = model_setup$mcmcParams$nchains, 
                         niter = model_setup$mcmcParams$niter, 
                         nburnin = model_setup$mcmcParams$nburn, 
                         thin = model_setup$mcmcParams$nthin, 
                         samplesAsCodaMCMC = TRUE, 
                         setSeed = mySeed)
  return(IDSM_out)
}


# # Define model
# myModel <- nimbleModel(code = IDSM_code,
#                        data = input_data$nim.data, 
#                        constants = input_data$nim.constants,
#                        inits = model_setup$initVals[[3]])
# 
# myModel$initializeInfo()
# myModel$calculate()
# 
# # Configure model
# confo <- configureMCMC(int.Nmixture.model, monitors = parameters.to.save)
# # confo$removeSampler()
# # confo$addSampler()
# 
# # Build and compile MCMC
# CmyModel <- compileNimble(myModel)
# 
# myMCMC <- buildMCMC(CmyModel, 
#                     monitors = model_setup$modelParams)
# 
# CmyMCMC <- compileNimble(myMCMC)
# 
# # Run
# IDSM_out <- runMCMC(CmyMCMC, 
#                    niter = model_setup$mcmcParams$niter, 
#                    nburnin = model_setup$mcmcParams$nburn, 
#                    thin = model_setup$mcmcParams$nthin, 
#                    samplesAsCodaMCMC = TRUE, 
#                    setSeed = per_chain_info$mySeed)
# 
