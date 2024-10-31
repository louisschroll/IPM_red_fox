library(ggplot2)
library(nimble)
library(sf)
library(reshape2)
library(remotes)
library(ckanr)
library(purrr)
library(dplyr)
library(metafor)
library(patchwork)

#**********#
# 0) SETUP #
#**********#

## Set seed
mySeed <- 10

## Set general parameters
Amax <- 5 # Number of age classes
Tmax <- 18  # Number of years
Tmax_sim <- 10
minYear <- 2004 # First year to consider
maxAge_yrs <- 10 # Age of the oldest female recorded
summer_removal <- c(6,7,8,9)    #removal of summer months: numerical months to be removed from winter age at harvest data
winter_removal <- c(1:6, 10:12) #removal of winter months: numerical months to be removed from summer age at harvest data
area_selection<- c("Inner", "BB",  "Tana")# choosing varanger sub area ("Inner" / "BB" / "Tana)     ((BB = Batsfjord and Berlevag areas))
# start and end of placental scars and embryo sample periods (julian day)
plac_start <- 180 #including
plac_end   <- 80  #until, not including
embr_start <- 100 #including
embr_end   <- 140 #until, not including

## set dataset names, versions, and directories, and access
carcass.dataset.name <- "v_redfox_carcass_examination_v3"
carcass.dataset.version <- 3

rodent.dataset.name <-"v_rodents_snaptrapping_abundance_regional_v5"
rodent.dataset.version <- 5

# Stijn
shapefile.dir <- "C:\\Users\\sho189\\OneDrive - UiT Office 365\\PhD\\RedfoxIPM\\Fox areas shapefile\\tana rest"
COAT_key <- Sys.getenv("API_COAT_Stijn") # Stijn's API key for the COAT dataportal is saved as an environmental variable on the computer 

# Chloe
#shapefile.dir <- "C:/Users/chloe.nater/OneDrive - NINA/Documents/Projects/RedFox_IPM/Data/shapefiles"
shapefile.dir <- "Data/shapefiles"
COAT_key <- Sys.getenv("COAT_API")

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
fitCov.mH <- FALSE # Fit covariates on mH (harvest effort)
fitCov.mO <- TRUE # Fit covariates on mO (rodent abundance x reindeer carcasses)
fitCov.Psi <- TRUE # Fit covariates on Psi (rodent abundance)
fitCov.rho <- TRUE # Fit covariates on rho (rodent abundance)
fitCov.immR <- TRUE # Fit covariates on immigration rate (rodent abundance) - only if immigration is estimated as a rate
rCov.idx <- FALSE # Use discrete vs. continuous rodent covariate
nLevels.rCov <- 2 # 2-level discrete rodent covariate
#nLevels.rCov <- 3 # 3-level discrete rodent covariate (data not currently prepared)
standSpec.rCov <- TRUE # standardize different rodent species before summing (offset catchability) v.s. simply sum all numbers
reinCov.VarTana <- TRUE # Calculate the reindeer carcass data count covariate using Varanger (+Tana) municipalities as geographical area. FALSE is for whole of Eastern Finnmark

# Random year effect toggles
mO.varT <- TRUE

# Age-at-harvest data toggles
add.sumr.unaged <- FALSE # Add summer harvested individuals as un-aged individuals to the total harvested individuals in winter
saAH.years <- c(2005:2012) # Years for which the summer age at harvest matrix should be constructed (e.g. years in which summer harvest was aged consistently)

# Annual survival prior type toggles
HoenigPrior <- FALSE # Use prior on natural mortality derived from Hoenig model
#sPriorSource <- "Bristol" # Base survival prior on data from Bristol (not hunted)
#sPriorSource <- "NSweden" # Base survival prior on data from North Sweden (lightly hunted)
sPriorSource <- "metaAll" # Base survival prior on meta-analysis including all populations
#sPriorSource <- "metaSub" # Base survival prior on meta-analysis including only not/lightly hunted populations

# Immigration parameters toggle
imm.asRate <- TRUE # Estimating immigration as a rate as opposed to numbers

# Genetic immigration data toggles (details in documentation of wrangleData_gen function
poolYrs.genData <- TRUE # Pool data across all years
useData.gen <- TRUE # Use genetic data for estimation of immigration rate
indLikelihood.genData <- FALSE # Apply an individual-level likelihood for genetic data
threshold <- 0.05
#threshold <- 0.2
#pImm.type <- "original"
pImm.type <- "rescaled"
#pImm.type <- "LL-based"

# Den survey prior and data toggles
useData.pup <- TRUE
useInfPrior.S0 <- FALSE

## Changes to denning survival prior
S0.mean.offset <- 0
S0.sd.factor <- 1

## Set up perturbation parameters for running standard scenarios
pert.mH <- FALSE
pert.mO <- FALSE
pert.S0 <- FALSE
pert.immR <- FALSE
pert.rodent <- FALSE
pert.reindeer <- FALSE

factor.mH <- 1
factor.mO <- 1
factor.S0 <- 1
factor.immR <- 1
factor.rodent <- 1
factor.reindeer <- 1

if(pert.mH & factor.mH == 0){
  pert.mHs <- TRUE
  factor.mHs <- 0
}else{
  pert.mHs <- FALSE
  factor.mHs <- 1
}

perturbVecs <- setupPerturbVecs_PVA(Tmax = Tmax, Tmax_sim = Tmax_sim,
                                    pert.mH = pert.mH, factor.mH = factor.mH,
                                    pert.mO = pert.mO, factor.mO = factor.mO,
                                    pert.S0 = pert.S0, factor.S0 = factor.S0,
                                    pert.mHs = pert.mHs, factor.mHs = factor.mHs,
                                    pert.immR = pert.immR, factor.immR = factor.immR,
                                    pert.rodent = pert.rodent, factor.rodent = factor.rodent,
                                    pert.reindeer = pert.reindeer, factor.reindeer = factor.reindeer)

## Set up perturbation parameters for running rodent-dependent harvest scenarios
factor.mH.rodent <- 1.5
threshold.rodent.mH <- 0
thresholdAbove <- FALSE


## Nimble function for determining perturbation factor based on covariate value
calculate_pertFac <- nimbleFunction(
  
  run = function(pertFactor = double(0),
                 covThreshold = double(0),
                 thresholdAbove = logical(0),
                 covValue = double(0)) {
    
    # Set conditional perturbation factor
    if((thresholdAbove & covValue > covThreshold) | (!thresholdAbove & covValue < covThreshold)){
      pertFac <- pertFactor
    }else{
      pertFac <- 1
    }
    
    # Return perturbation factor
    return(pertFac)
    returnType(double(0))
  }
)


#*********************#
# 1) DATA PREPARATION #
#*********************#

# 1a) Download and reformat carcass data
#-------------------------------#

## Download carcass data
carcass.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = carcass.dataset.name,
                                     COATdataset.version = carcass.dataset.version)

## Reformat carcass data
carcass.data <- reformatData_carcass(Amax = Amax,   
                                     summer_removal = summer_removal ,
                                     winter_removal = winter_removal ,
                                     area_selection = area_selection,
                                     plac_start = plac_start,
                                     plac_end = plac_end ,
                                     embr_start = embr_start ,
                                     embr_end = embr_end,
                                     carcass.dataset = carcass.data.raw,
                                     shapefile.dir = shapefile.dir,
                                     add.sumr.unaged = add.sumr.unaged, 
                                     saAH.years = saAH.years)


# 1b) Age-at-Harvest data #
#--------------------------------#

## Winter AaH data
wAaH.data <- wrangleData_AaH(AaH.datafile = carcass.data$WAaH.matrix, 
                             Amax = Amax)
## Summer AaH data
sAaH.data <- wrangleData_AaH(AaH.datafile = carcass.data$SAaH.matrix, 
                             Amax = Amax)


# 1c) Reproduction data #
#-----------------------#

## Set data paths/filenames
P1.datafile <- carcass.data$P1var # Placental scar/embryo count
P2.datafile <- carcass.data$P2var # Presence of placental scars/embryos/pregnancy signs

## Prepare reproduction data
rep.data <- wrangleData_rep(P1.datafile = P1.datafile, 
                            P2.datafile = P2.datafile,
                            Amax = Amax, 
                            minYear = minYear)


# 1d) Genetic data #
#------------------#

## Set data paths
genetics.datapath <- "Data/RedFox_genetics_immigrant_probabilities.txt"

## Prepare genetic data
gen.data <- wrangleData_gen(datapath = genetics.datapath,
                            minYear, 
                            onlyFemales = FALSE, 
                            poolYrs.genData = poolYrs.genData, 
                            threshold = threshold)


# 1e) Opportunistic pup observation data #
#----------------------------------------#

## Set data path
pups.datapath <- "Data/Rfox_early_litter_sizes.csv"

## Prepare pup observation data
pup.data <- wrangleData_pup(datapath = pups.datapath,
                            minYear = minYear)


# 1f) Harvest effort data #
#-------------------------#

## Prepare harvest effort data
hunter.data <- reformatData_hunters(area_selection = area_selection,
                                    carcass.dataset = carcass.data.raw,
                                    shapefile.dir = shapefile.dir)


# 1g) Environmental data #
#------------------------#

## Download rodent data
rodent.data.raw <- downloadData_COAT(COAT_key = COAT_key, 
                                     COATdataset.name = rodent.dataset.name,
                                     COATdataset.version = rodent.dataset.version)

## Reformat rodent data
rodent.data <- reformatData_rodent(rodent.dataset = rodent.data.raw,
                                          minYear = minYear)

## Reformat reindeer data
reindeer.data <- reformatData_reindeer(minYear = minYear,
                                       Tmax = Tmax,
                                       reinCov.VarTana = reinCov.VarTana)


# 1h) Conceptual year information #
#---------------------------------#

YearInfo <- collate_yearInfo(minYear = minYear,
                             Tmax = Tmax)


#**********************#
# 2) PRIOR INFORMATION #
#**********************#

## Parameters/paths for making informative priors for survival based on meta-analysis of literature data
meta.datafile <- "Data/RedFox_LiteratureData.csv"
simulateSD <- TRUE

## Parameters/paths for making informative priors for natural mortality using Tom Porteus' Hoenig model approach
mu.t.max <- 22.61062
hoenig.datafile <- "Data/HoenigMod_Posteriors_fromTomPorteus.txt"
nsim <- 30

## Collate all prior information
surv.priors <- collate_priorInfo(meta.datafile = meta.datafile,
                                 simulateSD = simulateSD,
                                 hoenig.datafile = hoenig.datafile, 
                                 nsim = nsim, 
                                 mu.t.max = mu.t.max, 
                                 maxAge = maxAge_yrs,
                                 S0.mean.offset = S0.mean.offset,
                                 S0.sd.factor = S0.sd.factor)

## Define type of prior to use for annual survival
survPriorType <- definePriorType_AnnSurv(HoenigPrior = HoenigPrior, 
                                         sPriorSource = sPriorSource)

#****************#
# 3) MODEL SETUP #
#****************#

# 3a) Write model code #
#----------------------#

redfox.code <- writeCode_redfoxIPM_PVA(indLikelihood.genData = indLikelihood.genData)


# 3b) Assemble IPM input data #
#-----------------------------#

input.data <- assemble_inputData_PVA(Amax = Amax, 
                                     Tmax = Tmax, 
                                     Tmax_sim = Tmax_sim,
                                     minYear = minYear,
                                     maxPups = 14,
                                     uLim.N = 800,
                                     uLim.Imm = 3000,
                                     nLevels.rCov = nLevels.rCov,
                                     standSpec.rCov = standSpec.rCov,
                                     poolYrs.genData = poolYrs.genData,
                                     pImm.type = pImm.type,
                                     wAaH.data = wAaH.data, 
                                     sAaH.data = sAaH.data,
                                     rep.data = rep.data, 
                                     gen.data = gen.data,
                                     pup.data = pup.data,
                                     rodent.data = rodent.data, 
                                     reindeer.data = reindeer.data,
                                     hunter.data = hunter.data, 
                                     surv.priors = surv.priors,
                                     survPriorType = survPriorType,
                                     perturbVecs = perturbVecs,
                                     factor.mH.rodent = factor.mH.rodent,
                                     threshold.rodent.mH = threshold.rodent.mH,
                                     thresholdAbove = thresholdAbove)


# 3c) Set up for model run (incl. simulating initial values) #
#------------------------------------------------------------#

model.setup <- setupModel_PVA(modelCode = redfox.code, 
                              nim.data = input.data$nim.data, 
                              nim.constants = input.data$nim.constants, 
                              minN1 = c(600, 50, 50, 50, 50), 
                              maxN1 = c(800, 400, 400, 400, 400), 
                              minImm = 50, 
                              maxImm = 600,
                              fitCov.mH = fitCov.mH, 
                              fitCov.mO = fitCov.mO,
                              fitCov.Psi = fitCov.Psi, 
                              fitCov.rho = fitCov.rho,
                              fitCov.immR = fitCov.immR,
                              rCov.idx = rCov.idx,
                              mO.varT = mO.varT,
                              HoenigPrior = HoenigPrior,
                              imm.asRate = imm.asRate,
                              testRun = FALSE,
                              initVals.seed = mySeed)


####################
# 4) MODEL FITTING #
####################

t1 <- Sys.time()
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
Sys.time() - t1

saveRDS(IPM.out, file = "RedFoxIPM_sim_baseline.rds") # No perturbation
#saveRDS(IPM.out, file = "RedFoxIPM_sim_noHarvest.rds") # pert.mH = TRUE, mH.factor = 0
#saveRDS(IPM.out, file = "RedFoxIPM_sim_higherHarvest_fac1.5.rds") # pert.mH = TRUE, mH.factor = 1.5 (initVals.seed = mySeed + 2 = 12)
#saveRDS(IPM.out, file = "RedFoxIPM_sim_lowRodentHarvestMatch_th0_fac1.50.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_sim_highRodentHarvestMatch_th0_fac1.50.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_sim_highRodentHarvestDelay_th0_fac1.50.rds")
#saveRDS(IPM.out, file = "RedFoxIPM_sim_lowRodentHarvestDelay_th0_fac1.50.rds")

#MCMCvis::MCMCtrace(IPM.out)


########################
# 5) MODEL COMPARISONS #
########################

## Baseline vs. no harvest
PVA1_comp <- compareModels(Amax = Amax, 
                           Tmax = Tmax, 
                           minYear = minYear, 
                           maxYear = 2030,
                           logN = TRUE,
                           post.filepaths = c("RedFoxIPM_sim_baseline.rds", 
                                              "RedFoxIPM_sim_noHarvest.rds",
                                              "RedFoxIPM_sim_higherHarvest_fac1.5.rds"), 
                           model.names = c("Baseline", 
                                           "No harvest",
                                           "50% higher harvest"), 
                           plotFolder = "Plots/ScenarioComp_PVA1_RodDyn",
                           returnSumData = TRUE)

# Extra plot for manuscript: 
maxYear <- 2030
pdf("Plots/ScenarioComp_PVA1_RodDyn/PosteriorSummaries_TimeSeries_Ntot.pdf", width = 8, height = 4)
print(
  PVA1_comp %>%
    dplyr::mutate(Action = dplyr::case_when(Model == "Baseline" ~ "None",
                                            Model == "50% higher harvest" ~ "+50% harvest",
                                            Model == "No harvest" ~ "Stop harvest")) %>%
    dplyr::filter(ParamName == "N.tot" & Year > minYear & Year <= maxYear) %>%
    ggplot(aes(x = Year, group = Model)) + 
    geom_line(aes(y = median, color = Action)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Action), alpha = 0.2) + 
    scale_fill_manual(values = c("#087792", "grey70", "#E3416F")) + 
    scale_color_manual(values =  c("#087792", "grey70", "#E3416F")) + 
    scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
    ylab("Log population size") +
    ggtitle("Female population size") +  
    theme_bw() + theme(panel.grid.minor = element_blank(), 
                       panel.grid.major.y = element_blank(), 
                       axis.text.x = element_text(angle = 45, vjust = 0.5))
)
dev.off()


## Different harvest scenario types
PVA2_comp <- compareModels(Amax = Amax, 
                           Tmax = Tmax, 
                           minYear = minYear, 
                           maxYear = 2030,
                           logN = TRUE,
                           post.filepaths = c("RedFoxIPM_sim_baseline.rds", 
                                              #"RedFoxIPM_sim_incHarvest_fac1.50.rds",
                                              "RedFoxIPM_sim_highRodentHarvestMatch_th0_fac1.50.rds",
                                              "RedFoxIPM_sim_lowRodentHarvestMatch_th0_fac1.50.rds",
                                              "RedFoxIPM_sim_highRodentHarvestDelay_th0_fac1.50.rds",
                                              "RedFoxIPM_sim_lowRodentHarvestDelay_th0_fac1.50.rds"#,
                           ), 
                           model.names = c("Baseline projection", 
                                           #"50% harvest increase",
                                           "Higher rodent +50% harvest, match",
                                           "Lower rodent +50% harvest, match",
                                           "Higher rodent +50% harvest, delay",
                                           "Lower rodent +50% harvest, delay"
                           ), 
                           plotFolder = "Plots/ScenarioComp_PVA2_RodDyn",
                           returnSumData = TRUE)

# Extra plot for manuscript: 
maxYear <- 2030
scenInfo <- data.frame(Action = c("None", 
                                  rep("High rodent +50% harvest", 2),
                                  rep("Low rodent +50% harvest", 2)), 
                       Timing = c("matched", 
                                  rep(c("delayed", "matched"), 2)),
                       Model = unique(PVA3_comp$Model))
scenInfo$Action <- factor(scenInfo$Action, levels = c("None", "High rodent +50% harvest", "Low rodent +50% harvest"))
scenInfo$Timing <- factor(scenInfo$Timing, levels = c("matched", "delayed"))

pdf("Plots/ScenarioComp_PVA2_RodDyn/PosteriorSummaries_TimeSeries_Ntot.pdf", width = 8, height = 4)
  print(
    PVA2_comp %>%
      dplyr::filter(ParamName == "N.tot" & Year > minYear & Year <= maxYear) %>%
      dplyr::left_join(scenInfo, by = "Model") %>%
      ggplot(aes(x = Year, group = Model)) + 
      geom_line(aes(y = median, color = Action, linetype = Timing)) + 
      geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Action), alpha = 0.1) + 
      scale_fill_manual(values = c("grey70", "#089392FF", "#E3756FFF")) + 
      scale_color_manual(values =  c("grey70", "#089392FF", "#E3756FFF")) + 
      scale_linetype_manual(values = c("solid", "dashed")) +
      scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
      ylab("Log population size") + 
      ggtitle("Female population size") +  
      theme_bw() + theme(panel.grid.minor = element_blank(), 
                         panel.grid.major.y = element_blank(), 
                         axis.text.x = element_text(angle = 45, vjust = 0.5))
  )
dev.off()
