#' Collate prior information for use in the model
#'
#' @param meta.datafile character string. Path file name for file containing literature data extracted from Devenish-Nelson et al. 2012
#' @param simulateSD logical. If TRUE (default), missing standard deviations for literature data are simulated using truncated normal distributions (with sd = 1).
#' If FALSE, missing standard deviations are instead calculated as within-age-class variation in raw literature data.
#' @param hoenig.datafile character string. Path/file name for file containing posterior
#' samples from Tom Porteus' run of the phylogenetic Hoenig model. Used by the
#' dependency function predict_mO_HoenigMod().
#' @param mu.t.max  numeric. Offset parameter for maximum age in the Hoenig 
#' model (= 22.61062). Provided by Tom Porteus. 
#' @param maxAge integer. Maximum recorded age of harvested foxes. 
#' @param nsim integer. Number of simulation replicates for each posterior sample.
#' @param S0.mean.offset numeric. Absolute change in average denning survival 
#' relative to estimate from arctic fox model. Default = 1 (no change). 
#' @param S0.sd.factor numeric.  Relative change in standard deviation of 
#' denning survival relative to estimate from arctic fox model. Default = 0 (no
#' change).
#'
#' @return a list of lists containing parameters to define informative priors
#' for early survival, age-specific annual survival, and juvenile/adult natural
#' mortality hazard rate.
#' @export
#'
#' @examples

collate_priorInfo <- function(meta.datafile, simulateSD = TRUE, hoenig.datafile, mu.t.max, maxAge, nsim, S0.mean.offset = 0, S0.sd.factor = 1){
  
  
  ## Simulate / load prior distribution parameters from Hoenig model
  if(file.exists("mO_prior_Parameters.rds")){
    mO.prior <- readRDS("mO_prior_Parameters.rds")
  }else{
    mO.prior <- predict_mO_HoenigMod(hoenig.datafile, mu.t.max, maxAge, nsim = nsim, plot = FALSE)
  }
  
  ## Simulate prior distribution parameters via meta-analysis
  metaS_all <- estimate_metaS(meta.datafile, lowHarvestOnly = FALSE, simulateSD = simulateSD)
  metaS_lowHarvest <- estimate_metaS(meta.datafile, lowHarvestOnly = TRUE, simulateSD = simulateSD)
  
  ## Assemble prior data
  prior.data <- list(
    
    # Early survival (denning period) - from arctic fox IPM
    earlySurv = list(
      S0.mean = 0.74 + S0.mean.offset, # From arctic fox IPM: 0.7428115
      S0.sd = 0.06 * S0.sd.factor #  IPM: 0.05983706
    ),

    
    # Annual natural survival
    annSurv = list(
      #* Literature, Bristol (non-hunted)
      Bristol = list(
        Snat.mean = c(0.48, 0.54, 0.53, 0.51, 0.51),
        Snat.sd = c(0.02, 0.03, 0.03, 0.03, 0.03)
      ),
      
      
      #* Literature, North Sweden (lightly hunted)
      NSweden = list(
        Snat.mean = c(0.33, 0.71, 0.50, 0.59, 0.59),
        Snat.sd = c(0.02, 0.04, 0.05, 0.04, 0.04)
      ),
      
      
      #*Literature meta-analysis (all)
      metaAll = list(
        Snat.mean = metaS_all$Estimate,
        Snat.sd = metaS_all$se
      ),
      
      
      #*Literature meta-analysis (non/lightly hunted)
      metaSub = list(
        Snat.mean = metaS_lowHarvest$Estimate,
        Snat.sd = metaS_lowHarvest$se
      )
    ),
    
    
    # Natural mortality (Base parameters from Hoenig model, ratio from arctic fox IPM)
    natMort = list(
      mnat.logmean = mO.prior$mnat.logmean,
      mnat.logsd = mO.prior$mnat.logsd,
      ratioJA.logmean = 0.4807439,
      ratioJA.logsd = 0.355224
    )
  )
  
  return(prior.data)
}