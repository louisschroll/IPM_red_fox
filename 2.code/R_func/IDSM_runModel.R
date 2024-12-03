


IDSM_runModel <- function(DS_data_list, 
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
                                 niter = 50000, 
                                 nthin = 5, 
                                 nburn = 20000, 
                                 nchains = 4,
                                 testRun = TRUE,
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
