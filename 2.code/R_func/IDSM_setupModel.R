#' Set up model, data, and initial values for running MCMC
#'
#' @param modelCode an R call object specifying the model structure for integrated
#' distance sampling model
#' @param nim.data list of input objects representing data
#' @param nim.constants list of input objects representing constants
#' @param R_perF logical. If TRUE, treats recruitment rate as juvenile per adult female.
#' If FALSE, treats recruitment rate as juvenile per adult (sum of both sexes).
#' @param niter integer. Number of MCMC iterations (default = 25000)
#' @param nthin integer. Thinning factor (default = 5)
#' @param nburn integer. Number of iterations to discard as burn-in (default = 5000)
#' @param nchains integer. Number of chains to run.
#' @param testRun logical. If TRUE, sets up for a test run with 10 iterations,
#' no thinning, and no burn-in (default = FALSE)
#' @param initVals.seed integer. Seed to use for inital value simulation.
#'
#' @return list of list containing all components necessary for running model
#' with `nimble::nimbleMCMC()`

IDSM_setupModel <- function(nim.data,
                            nim.constants,
                            niter = 20000,
                            nthin = 3,
                            nburn = 1000,
                            nchains = 3,
                            testRun = FALSE,
                            initVals.seed) {
  
  ## Set parameters to monitor
  params <- c("esw", "p", "mean.recruitment", "N_tot_gic", 
              "sigma", "meanDens", "mean.survival",
              "harvest_rate", "N_tot_gic_age", "N_tot_exp")
  
  ## Simulate initial values
  #set.seed(initVals.seed)
  initVals <- list()
  for (c in 1:nchains) {
    initVals[[c]] <- IDSM_simulateInits(
      nim.data = nim.data,
      nim.constants = nim.constants
    )
  }
  
  ## Adjust MCMC parameters if doing a test run
  if (testRun) {
    niter <- 50
    nthin <- 1
    nburn <- 0
  }
  
  ## Collate model setup variables in a list
  setup <- list(
    modelParams = params,
    initVals = initVals,
    mcmcParams = list(
      niter = niter,
      nthin = nthin,
      nburn = nburn,
      nchains = nchains
    )
  )
  
  return(setup)
}