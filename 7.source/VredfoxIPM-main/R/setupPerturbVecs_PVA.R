#' Set up perturbation vectors for Population Viability Analysis (PVA) scenarios
#'
#' This function sets up perturbation vectors for different vital rates and
#' environmental covariates for use in model simulations following the data
#' collection / study period. 
#'
#' For vital rates, the selected perturbation factor is implemented as a
#' proportional change in the predicted yearly values. For example, if a
#' perturbation factor of 1.2 is provided for harvest mortality, this will effect
#' a 20% increase in harvest mortality in the simulation period. Conversely, if
#' a perturbation factor of 0.8 is provided, this will result in a 20% decrease
#' in harvest mortality in the simulation period. 
#'
#' For environmental covariates, the perturbation factor does not directly work 
#' as a multiplier as this would only affect the spread of simulated covariate
#' values and not the average for covariates that are z-standardized like here. 
#' Consequently, the perturbation factor is interpreted as an offset to the mean
#' in the simulation of covariate values in the simulation period. For example, 
#' a perturbation factor of 1.2 will result in an increase of 0.2 (1.2-1) of the 
#' expected mean covariate value in the simulated period. Conversely, a 
#' perturbation factor of 0.8 will decrease the exoected mean by 0.2 (0.8-1).
#
#' @param Tmax integer. The number of years in the analysis.
#' @param Tmax_sim integer. The number of years to consider for simulations 
#' beyond the data collection period. 
#' @param pert.mH logical. Whether to apply a perturbation to winter harvest 
#' mortality hazard rate. Default = FALSE. 
#' @param factor.mH numeric. Relative change to winter harvest mortality hazard 
#' rate to apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.mO logical. Whether to apply a perturbation to natural mortality
#' hazard rate. Default = FALSE. 
#' @param factor.mO numeric. Relative change to natural mortality hazard rate to 
#' apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.S0 logical. Whether to apply a perturbation to denning survival. 
#' Default = FALSE. 
#' @param factor.S0 numeric. Relative change to denning survival to apply.
#' 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.mHs logical. Whether to apply a perturbation to summer harvest 
#' mortality hazard rate. Default = FALSE. 
#' @param factor.mHs numeric. Relative change to summer harvest mortality hazard 
#' rate to apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.immR logical. Whether to apply a perturbation to immigration rate. 
#' Default = FALSE. 
#' @param factor.immR numeric. Relative change to immigration rate to apply.
#' 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.rodent logical. Whether to apply a perturbation to rodent covariate. 
#' Default = FALSE. 
#' @param factor.rodent numeric. Change to standardized rodent covariate 
#' mean to apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#' @param pert.reindeer logical. Whether to apply a perturbation to reindeer 
#' covariate. Default = FALSE. 
#' @param factor.reindeer numeric. Change to standardized reindeer 
#' covariate mean to apply. 1 = no change (default). < 1 = decrease. > 1 = increase. 
#'
#' @return
#' @export
#'
#' @examples

setupPerturbVecs_PVA <- function(Tmax, Tmax_sim,
                                 pert.mH = FALSE, factor.mH = 1,
                                 pert.mO = FALSE, factor.mO = 1,
                                 pert.S0 = FALSE, factor.S0 = 1,
                                 pert.mHs = FALSE, factor.mHs = 1,
                                 pert.immR = FALSE, factor.immR = 1,
                                 pert.rodent = FALSE,  factor.rodent = 1,
                                 pert.reindeer = FALSE, factor.reindeer = 1){
  
  ## Check that there are no invalid perturbation factors
  if(any(c(factor.mH, factor.mO, factor.S0, factor.mHs, factor.immR, factor.rodent, factor.reindeer) < 0)){
    stop("Invalid perturbation factor provided. Perturbation factors have to be
         numerical values >= 0.")
  }
  
  if(any(!is.numeric(c(factor.mH, factor.mO, factor.S0, factor.mHs, factor.immR, factor.rodent, factor.reindeer)))){
    stop("Invalid perturbation factor provided. Perturbation factors have to be
         numerical values >= 0.")
  }
  
  ## Set up basics perturbation vectors during study/data period
  pertFac.mH <- pertFac.mO <- pertFac.mHs <- pertFac.immR <- rep(1, Tmax)
  pertFac.S0 <- rep(1, Tmax+1)
  pertFac.rodent <- pertFac.reindeer <- rep(1, Tmax+1)
  
  ## Add factors for perturbation period
  if(Tmax_sim > 0){
    pertFac.mH <- c(pertFac.mH, rep(ifelse(pert.mH, factor.mH, 1), Tmax_sim))
    pertFac.mO <- c(pertFac.mO, rep(ifelse(pert.mO, factor.mO, 1), Tmax_sim))
    pertFac.S0 <- c(pertFac.S0, rep(ifelse(pert.S0, factor.S0, 1), Tmax_sim))
    pertFac.mHs <- c(pertFac.mHs, rep(ifelse(pert.mHs, factor.mHs, 1), Tmax_sim))
    pertFac.immR <- c(pertFac.immR, rep(ifelse(pert.immR, factor.immR, 1), Tmax_sim))
    pertFac.rodent <- c(pertFac.rodent, rep(ifelse(pert.rodent, factor.rodent, 1), Tmax_sim))
    pertFac.reindeer <- c(pertFac.reindeer, rep(ifelse(pert.reindeer, factor.reindeer, 1), Tmax_sim))
  }
  
  ## List and return perturbation vectors
  pertVecs <- list(pertFac.mH = pertFac.mH, pertFac.mO = pertFac.mO, 
                   pertFac.S0 = pertFac.S0, pertFac.mHs = pertFac.mHs, 
                   pertFac.immR = pertFac.immR,
                   pertFac.rodent = pertFac.rodent, pertFac.reindeer = pertFac.reindeer)
  return(pertVecs)
}