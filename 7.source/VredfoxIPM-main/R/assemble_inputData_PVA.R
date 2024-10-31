#' Assemble demographic and environmental data for running the IPM-PVA
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param Tmax_sim integer. The number of years to consider for simulations 
#' beyond the data collection period. 
#' @param minYear integer. First year to consider in analyses.
#' @param maxPups integer. Upper prior bound for average litter size.
#' @param uLim.N integer. Upper prior bound for initial number of individuals per age class.
#' @param nLevels.rCov integer. Number of levels of categorical rodent abundance to use.
#' @param standSpec.rCov logical. If TRUE, standardises rodent numbers per species before summing 
#' to offset catchability, If FALSE simple sums alls rodent numbers. 
#' @param poolYrs.genData integer. Whether or not genetic immigration data is pooled across years.
#' @param pImm.type character. Which type of individual-level data to use for immigration. 
#' "original" = p values as output by Geneclass 2. "rescaled" = p values as output
#' by Geneclass 2 and standardized so that the minimum immigrant probability = 0.
#' "LL-based" = log likelihood other / log likelihood other + log likelihood Varanger. 
#' @param uLim.Imm integer. Upper prior bound for annual number of immigrants. 
#' @param wAaH.data a list containing a winter Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param sAaH.data a list containing a summer Age-at-Harvest matrix (C) and a vector of
#' yearly proportions of individuals aged/included in Age-at-Harvest data (pData).
#' @param rep.data a list containing formatted reproduction data in two data 
#' frames: 'P1' (counts) and 'P2' (presences/absences).
#' @param gen.data a list containing relevant data on genetically determined 
#' probabilities of individuals being immigrants.
#' @param pup.data a list containing data on numbers of pups observed on dens. 
#' @param rodent.data a list containing rodent abundance data as a continuous variable (cont),
#' and categorical variable with two (cat2) and three (cat3) levels.
#' @param reindeer.data a list containing reindeer carcass abundance and proportion
#' of foxes with reindeer in stomachs.
#' @param hunter.data a dataframe containing original and scaled counts of successful 
#' hunters per year.
#' @param surv.priors a list of lists containing parameters to define informative priors
#' for early survival, age-specific annual survival, and juvenile/adult natural
#' mortality hazard rate.
#' @param survPriorType a list containing information on prior for annual survival.
#' @param perturbVecs a list of perturbation vectors for vital rates, each 
#' with a length of Tmax+Tmax_sim or Tmax+Tmax_sim+1 and made up of positive 
#' numerics. The perturbation vectors have to be named pertFac.[X], where [X] =
#' mH, mO, S0, and immR. Output of function setupPerturbVecs. 
#' @param factor.mH.rodent numeric. Multiplicative perturbation factor for 
#' harvest mortality to be apply if rodent covariate is above a defined threshold.
#' @param threshold.rodent.mH numeric. Threshold for rodent covariate value
#' (z-standardized) above which to apply harvest perturbation.
#' @param thresholdAbove logical. If TRUE (default), applies harvest perturbation 
#' if rodent covariate value > threshold. If FALSE, applies harvest perturbation if rodent
#' covariate value < threshold.
#' @param save logical. If TRUE, saves assembled data as an .rds file in the 
#' working directory. Default = FALSE. 
#'
#' @return a list containing all data necessary for running the IPM. 
#' @export
#'
#' @examples

assemble_inputData_PVA <- function(Amax, Tmax, Tmax_sim, minYear,
                                   maxPups, uLim.N, uLim.Imm, 
                                   nLevels.rCov = NA, standSpec.rCov,
                                   poolYrs.genData, pImm.type,
                                   wAaH.data, sAaH.data, rep.data, gen.data, pup.data,
                                   rodent.data, reindeer.data, hunter.data, 
                                   surv.priors, survPriorType,
                                   perturbVecs, 
                                   factor.mH.rodent, threshold.rodent.mH,
                                   thresholdAbove = TRUE,
                                   save = FALSE){
  
  ## Select relevant years from observational data
  
  # Winter Age-at-Harvest data
  C_w <- wAaH.data$C[,which(colnames(wAaH.data$C) == minYear) + 1:Tmax - 1]
  pData_w <- wAaH.data$pData[which(colnames(wAaH.data$C) == minYear) + 1:Tmax - 1]
  
  # Summer Age-at-Harvest data
  C_s <- sAaH.data$C
  pData_s <- sAaH.data$pData
  XsH <- length(pData_s)
  sH_year <- as.numeric(colnames(C_s)) - minYear + 1
  
  # Reproduction data
  P1 <- subset(rep.data$P1, repryear %in% c(minYear + 0:Tmax))
  P2 <- subset(rep.data$P2, repryear %in% c(minYear + 0:Tmax))
  
  
  ## Select relevant categorical rodent covariate
  if(is.na(nLevels.rCov)){
    RodentIndex <- NA
    RodentIndex2 <- NA
  }else{
    
    if(nLevels.rCov == 2){
      RodentIndex <- rodent.data$cat2.wintvar
      RodentIndex2 <- rodent.data$cat2.fallstor
    }else{
      #RodentIndex <- rodent.data$cat3
      stop("3 level rodent covariate not currently supported")
    }
  }
  
  ## Select relevant continuous rodent covariate
  if(standSpec.rCov){
    RodentAbundance <- rodent.data$cont.wintvar.stsp
  }else{
    RodentAbundance <- rodent.data$cont.wintvar
  }
  
  if(standSpec.rCov){
    RodentAbundance2 <- rodent.data$cont.fallstor.stsp
  }else{
    RodentAbundance2 <- rodent.data$cont.fallstor
  }
  
  ## Select relevant reindeer covariates
  Reindeer <- reindeer.data$RDcarcass
  
  ## Add simulation years to covariates
  if(Tmax_sim > 0){
    RodentAbundance <- c(RodentAbundance, rep(NA, (Tmax+Tmax_sim+1-length(RodentAbundance))))
    RodentAbundance2 <- c(RodentAbundance2, rep(NA, (Tmax+Tmax_sim+1-length(RodentAbundance2))))
    RodentIndex <- c(RodentIndex, rep(NA, (Tmax+Tmax_sim+1-length(RodentIndex))))
    RodentIndex2 <- c(RodentIndex2, rep(NA, (Tmax+Tmax_sim+1-length(RodentIndex2))))
    Reindeer <- c(Reindeer, rep(NA, (Tmax+Tmax_sim+1-length(Reindeer))))
    HarvestEffort <- c(hunter.data$NHunters_std, rep(NA, (Tmax+Tmax_sim-length(hunter.data$NHunters_std))))
  }else{
    HarvestEffort <- hunter.data$NHunters_std
  }
  
  ## List all relevant data (split into data and constants as used by NIMBLE)
  # Data
  nim.data <- list(
    C_w = C_w,
    pData_w = pData_w,
    
    C_s = C_s,
    pData_s = pData_s,
    
    P1 = P1$P1,
    
    P2 = P2$P2,
    
    NoPups = pup.data$NoPups,
    
    HarvestEffort = HarvestEffort,
    RodentAbundance = RodentAbundance,
    RodentAbundance2 = RodentAbundance2,
    RodentIndex = RodentIndex,
    RodentIndex2 = RodentIndex2,
    Reindeer = Reindeer,
    
    pertFac.mH = perturbVecs$pertFac.mH,
    pertFac.mO = perturbVecs$pertFac.mO,
    pertFac.S0 = perturbVecs$pertFac.S0,
    pertFac.mHs = perturbVecs$pertFac.mHs,
    pertFac.immR = perturbVecs$pertFac.immR,
    pertFac.rodent = perturbVecs$pertFac.rodent,
    pertFac.reindeer = perturbVecs$pertFac.reindeer,
    
    factor.mH.rodent = factor.mH.rodent,
    threshold.rodent.mH = threshold.rodent.mH,
    thresholdAbove = thresholdAbove
  )
  
  # Constants
  nim.constants <- list(
    Amax = Amax,
    Tmax = Tmax,
    Tmax_sim = Tmax_sim,
    minYear = minYear,
    
    maxPups = maxPups,
    uLim.N = uLim.N,
    uLim.Imm = uLim.Imm,
    
    XsH = XsH,
    sH_year = sH_year,
    
    P1_age = P1$age_adj,
    P1_year = P1$RepYearIndex,
    X1 = length(P1$P1),
    
    P2_age = P2$age_adj,
    P2_year = P2$RepYearIndex,
    X2 = length(P2$P2),
    
    NoPups_year = pup.data$NoPups_year,
    X3 = pup.data$X3,
    
    nLevels.rCov = nLevels.rCov
  )
  
  ## Append relevant data from genetic immigration assignments
  if(poolYrs.genData){
    
    nim.data$pImm <- dplyr::case_when(pImm.type == "original" ~ gen.data$pImm,
                                      pImm.type == "rescaled" ~ gen.data$pImm_rescaled,
                                      pImm.type == "LL-based" ~ gen.data$pImm_LL)
    nim.constants$Xgen <- gen.data$Xgen
    
    nim.data$genObs_Imm <- gen.data$genObs_Imm
    nim.data$genObs_Res <- gen.data$genObs_Res
    
  }else{
    nim.data$pImm <- dplyr::case_when(pImm.type == "original" ~ gen.data$pImm_in,
                                      pImm.type == "rescaled" ~ gen.data$pImm_rescaled_in,
                                      pImm.type == "LL-based" ~ gen.data$pImm_LL_in)
    nim.data$pImm_pre <- dplyr::case_when(pImm.type == "original" ~ gen.data$pImm_pre,
                                          pImm.type == "rescaled" ~ gen.data$pImm_rescaled_pre,
                                          pImm.type == "LL-based" ~ gen.data$pImm_LL_pre)
    nim.constants$Xgen <- gen.data$Xgen_in
    nim.constants$Xgen_pre <- gen.data$Xgen_pre
    nim.constants$pImm_yrs <- gen.data$pImm_yrsB_in
    nim.constants$pImm_yrs_pre <- gen.data$pImm_yrsB_pre
    
    nim.constants$Tmax_Gen <- max(nim.constants$pImm_yrs)
    nim.constants$Tmax_Gen_pre <- max(nim.constants$pImm_yrs_pre)
    
    nim.data$genObs_Imm <- gen.data$genObs_Imm_in
    nim.data$genObs_Imm_pre <- gen.data$genObs_Imm_pre
    nim.data$genObs_Res <- gen.data$genObs_Res_in
    nim.data$genObs_Res_pre <- gen.data$genObs_Res_pre
    
    
  }
  
  ## Add relevant prior information
  nim.constants <- c(nim.constants, surv.priors$earlySurv)
  
  if(survPriorType$Parameter == "natMort"){
    nim.constants <- c(nim.constants, surv.priors$natMort)
  }else{
    sublistIdx <- which(names(surv.priors$annSurv) == survPriorType$Source)
    nim.constants <- c(nim.constants, surv.priors$annSurv[sublistIdx][[1]])
  }
  
  ## Combine data and constants in a list
  inputData <- list(nim.data = nim.data, 
                    nim.constants = nim.constants)
  
  ## Save (optional) and return data
  if(save){
    saveRDS("inputData_formatted.rds")
  }
  
  return(inputData)
  
}