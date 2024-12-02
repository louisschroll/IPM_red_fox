


IDSM_runModel <- function(DS_data_list, 
                          harvest_data, 
                          W, 
                          size_hunting_area){
  input_data <- IDSM_prepareInputData(DS_data_list = DS_data_list,
                                      harvest_data = harvest_data,
                                      W = W,
                                      size_hunting_area = size_hunting_area)
  
  IDSM_code <- IDSM_writeCode()
  
  model_setup <- IDSM_setupModel(nim.data = input_data$nim.data,
                                 nim.constants = input_data$nim.constants,
                                 niter = 50000, 
                                 nthin = 5, 
                                 nburn = 20000, 
                                 nchains = 4,
                                 testRun = TRUE,
                                 initVals.seed = mySeed)
  
  IDSM_out <- nimbleMCMC(code = model_setup$modelCode,
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
