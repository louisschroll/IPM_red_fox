library(ggplot2)
library(nimble)
library(tidyverse)

## Set toggle combos for models to run
prior_toggles <- list(
  HoenigPrior = c(TRUE, rep(FALSE, 4), rep(c(TRUE, FALSE), each = 4)),
  sPriorSource = c(NA, "Bristol", "NSweden", "metaAll", "metaSub", rep(c(NA, "NSweden"), each = 4))
)

cov_toggles <- list(
  fitCov.mH = c(rep(FALSE, 5), rep(c(TRUE, rep(FALSE, 3)), 2)),
  fitCov.Psi = c(rep(FALSE, 5), rep(c(FALSE, rep(TRUE, 3)), 2)),
  rCov.idx = c(rep(TRUE, 5), rep(c(rep(FALSE, 2), rep(TRUE, 2)), 2)),
  nLevels.rCov = c(rep(2, 5), rep(c(rep(2, 3), 3), 2))
)

## List model names
model_names <- c("Hoenig",
                 "Bristol",
                 "NSweden",
                 "metaAll", 
                 "metaSub",
                 "Hoenig-hunters",
                 "Hoenig-rodCont",
                 "Hoenig-rodIdx2",
                 "Hoenig-rodIdx3",
                 "NSweden-hunters",
                 "NSweden-rodCont",
                 "NSweden-rodIdx2",
                 "NSweden-rodIdx3")


for(i in 1:13){
  
  #**********#
  # 0) SETUP #
  #**********#
  
  ## Set seed
  mySeed <- 0
  
  ## Set general parameters
  Amax <- 5 # Number of age classes
  Tmax <- 15  # Number of years
  minYear <- 2004 # First year to consider
  
  maxAge_yrs <- 10 # Age of the oldest female recorded
  
  ## Source all functions in "R" folder
  sourceDir <- function(path, trace = TRUE, ...) {
    for (nm in list.files(path, pattern = "[.][RrSsQq]$")) {
      if(trace) cat(nm,":")
      source(file.path(path, nm), ...)
      if(trace) cat("\n")
    }
  }
  sourceDir('R')
  
  ## Set "switches" for running different model versions
  
  # Covariate toggles
  fitCov.mH <- cov_toggles$fitCov.mH[i]
  fitCov.Psi <- cov_toggles$fitCov.Psi[i]
  rCov.idx <- cov_toggles$rCov.idx[i]
  nLevels.rCov <- cov_toggles$nLevels.rCov[i]
  
  # Annual survival prior type toggles
  HoenigPrior <- prior_toggles$HoenigPrior[i]
  sPriorSource <- prior_toggles$sPriorSource[i]
  
  #*********************#
  # 1) DATA PREPARATION #
  #*********************#
  
  # 1a) Winter Age-at-Harvest data #
  #--------------------------------#
  
  ## Set data path/filename
  wAaH.datafile <- "Data/Cvar.tot_oct-mai_5age.txt"
  
  ## Prepare winter AaH data
  wAaH.data <- wrangleData_winterAaH(datafile = wAaH.datafile, 
                                     Amax = Amax)
  
  
  # 1b) Reproduction data #
  #-----------------------#
  
  ## Set data paths/filenames
  rep.datafiles <- c("Data/P1var_tot.txt", # Placental scar/embryo count
                     "Data/P2var_tot.txt") # Presence of placental scars/embryos/pregnancy signs
  
  ## Prepare reproduction data
  rep.data <- wrangleData_rep(datafiles = rep.datafiles, 
                              Amax = Amax, 
                              minYear = minYear)
  
  
  # 1c) Environmental data #
  #------------------------#
  
  ## Set data paths/filenames
  rodent.datafile <- "Data/stor_intensiv_04_20-year-var.txt"
  
  ## Prepare rodent abundance data
  rodent.data <- wrangleData_rodent(datafile = rodent.datafile,
                                    minYear = minYear,
                                    adjust = TRUE)
  
  ## Prepare harvest effort data
  hunter.data <- makeData_hunters()
  
  
  # 1d) Conceptual year information #
  #---------------------------------#
  
  YearInfo <- collate_yearInfo(minYear = minYear,
                               Tmax = Tmax)
  
  
  #**********************#
  # 2) PRIOR INFORMATION #
  #**********************#
  
  ## Make informative priors for natural mortality using Tom Porteus' Hoenig model approach
  mu.t.max <- 22.61062
  hoenig.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
  nsim <- 30
  
  
  ## Collate all prior information
  surv.priors <- collate_priorInfo(datafile = hoenig.datafile, 
                                   nsim = nsim, 
                                   mu.t.max = mu.t.max, 
                                   maxAge = maxAge_yrs)
  
  ## Define type of prior to use for annual survival
  survPriorType <- definePriorType_AnnSurv(HoenigPrior = HoenigPrior, 
                                           sPriorSource = sPriorSource)
  
  #****************#
  # 3) MODEL SETUP #
  #****************#
  
  # 3a) Write model code #
  #----------------------#
  
  redfox.code <- writeCode_redfoxIPM()
  
  
  # 3b) Assemble IPM input data #
  #-----------------------------#
  
  input.data <- assemble_inputData(Amax = Amax, 
                                   Tmax = Tmax, 
                                   minYear = minYear,
                                   maxPups = 14,
                                   uLim.N = 800,
                                   uLim.Imm = 800,
                                   nLevels.rCov = nLevels.rCov,
                                   wAaH.data = wAaH.data, 
                                   rep.data = rep.data, 
                                   rodent.data = rodent.data, 
                                   hunter.data = hunter.data, 
                                   surv.priors = surv.priors,
                                   survPriorType = survPriorType)
  
  
  # 3c) Set up for model run (incl. simulating initial values) #
  #------------------------------------------------------------#
  
  model.setup <- setupModel(modelCode = redfox.code, 
                            nim.data = input.data$nim.data, 
                            nim.constants = input.data$nim.constants, 
                            minN1 = c(600, 50, 50, 50, 50), 
                            maxN1 = c(800, 400, 400, 400, 400), 
                            minImm = 50, 
                            maxImm = 600,
                            fitCov.mH = fitCov.mH, 
                            fitCov.Psi = fitCov.Psi, 
                            rCov.idx = rCov.idx, 
                            HoenigPrior = HoenigPrior,
                            testRun = TRUE,
                            initVals.seed = mySeed)
  
  ####################
  # 4) MODEL FITTING #
  ####################
  
  IPM.out <- nimbleMCMC(code = model.setup$modelCode,
                        data = input.data$nim.data, 
                        constants = input.data$nim.constants,
                        inits = model.setup$initVals, 
                        monitors = model.setup$modelParams,
                        nchains = model.setup$mcmcParams$nchains, 
                        niter = model.setup$mcmcParams$niter, 
                        nburnin = model.setup$mcmcParams$nburn, 
                        thin = model.setup$mcmcParams$nthin, 
                        samplesAsCodaMCMC = TRUE, 
                        setSeed = 0)
  
  saveRDS(IPM.out, file = paste0(model_names[i], ".rds"))
  
}


#######################
# 5) MODEL COMPARISON #
#######################

post.filepaths <- list(
  paste0(model_names[1:5], ".rds"),
  paste0(model_names[2:5], ".rds"),
  paste0(model_names[c(3, 10:13)], ".rds"),
  paste0(model_names[c(1, 7:9)], ".rds")
)

model.names <- list(
  model_names[1:5],
  model_names[2:5],
  model_names[c(3, 10:13)],
  model_names[c(1, 7:9)]
)

plotFolder <- c(
  "Plots/Comp_Base",
  "Plots/Comp_BaseSurv",
  "Plots/Comp_SwedenCov",
  "Plots/Comp_HoenigCov"
)

for(n in 1:length(model.names)){
  compareModels(Amax = Amax, 
                Tmax = Tmax, 
                minYear = minYear,
                post.filepaths = post.filepaths[[n]], 
                model.names = model.names[[n]], 
                plotFolder = plotFolder[n])
}
