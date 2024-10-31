#' Run weighted regression on literature data for obtaining prior information on survival
#'
#' @param datafile character string. Path file name for file containing literature data extracted from Devenish-Nelson et al. 2012
#' @param lowHarvestOnly logical. If TRUE, only analyses studies that documented either none or only low levels of harvest.
#' Studies mentioning moderate to substantial harvest and studies for which no information on harvest level was available are excluded in this case.
#' @param simulateSD logical. If TRUE (default), standard deviations are simulated from a truncated normal distribution with mean = reported mean and sd = 1.
#' If FALSE, standard deviations are instead derived from variation of reported means for each age class. 
#' @param simSeed integer. The seed to use for simulating standard deviations. 
#'
#' @return
#' @export
#'
#' @examples

estimate_metaS <- function(datafile, lowHarvestOnly, simulateSD = TRUE, simSeed = 0){
  
  ## Set seed (for simulation)
  set.seed(simSeed) 
  
  ## Read in the data
  RedFox_LiteratureData <- readr::read_csv(datafile)
  
  ## Subset by include_MetaAnalysis (and lowHarvest)
  RedFox_MA <- RedFox_LiteratureData |> 
    dplyr::filter(include_MetaAnalysis==1)
  
  if(lowHarvestOnly){
    RedFox_MA <- RedFox_MA |>
      dplyr::filter(lowHarvest==1)
  }
  
  ## Impute/simulate SD

  if(simulateSD){
    
    RedFox_MA <- RedFox_MA |> 
      rowwise() |> 
      dplyr::mutate(Age0_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age0, lb=0, ub=1)),
                    Age1_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age1, lb=0, ub=1)),
                    Age2_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age2, lb=0, ub=1)),
                    Age3_sd=sd(TruncatedNormal::rtnorm(maxN_all,Age3, lb=0, ub=1))) 
    
    RedFox_MA$Age4_sd<-NULL
    
    for (i in 1:dim(RedFox_MA)[1]){
      if(is.na(RedFox_MA$Age4[i])){
        RedFox_MA$Age4_sd[i]=NA
      } else{
        RedFox_MA$Age4_sd[i]=sd(TruncatedNormal::rtnorm(RedFox_MA$maxN_all[i],RedFox_MA$Age4[i], lb=0, ub=1))
      }
    }
    
  }else{
    
    RedFox_MA <- RedFox_MA |>
      dplyr::mutate(Age0_sd=mean(Age0, na.rm=TRUE),
                    Age1_sd=mean(Age1, na.rm=TRUE),
                    Age2_sd=mean(Age2, na.rm=TRUE),
                    Age3_sd=mean(Age3, na.rm=TRUE),
                    Age4_sd=mean(Age4, na.rm=TRUE))
  }
  

  ## Weighted metaAnalysis
  Age0_MA<-rma(yi=Age0,vi=Age0_sd, weights=maxN_all, data=RedFox_MA)
  Age1_MA<-rma(yi=Age1,vi=Age1_sd, weights=maxN_all, data=RedFox_MA)
  Age2_MA<-rma(yi=Age2,vi=Age2_sd, weights=maxN_all, data=RedFox_MA)
  Age3_MA<-rma(yi=Age3,vi=Age3_sd, weights=maxN_all, data=RedFox_MA)
  Age4_MA<-rma(yi=Age4,vi=Age4_sd, weights=maxN_all, data=RedFox_MA)
  
  ## Plot estimates
  pdf.path <- ifelse(lowHarvestOnly, "Plots/Literature_MetaAnalysis_Estimates_lowHarvestOnly.pdf", "Plots/Literature_MetaAnalysis_Estimates.pdf")
  
  pdf(file = pdf.path, width = 5, height = 6)
    metafor::forest(Age0_MA, main = "Age 0")
    metafor::forest(Age1_MA, main = "Age 1")
    metafor::forest(Age2_MA, main = "Age 2")
    metafor::forest(Age3_MA, main = "Age 3")
    metafor::forest(Age4_MA, main = "Age 4")
  dev.off()

  ## Collect results in a dataframe
  out_dat <- tibble(Estimate=c(Age0_MA$beta,
                               Age1_MA$beta,
                               Age2_MA$beta,
                               Age3_MA$beta,
                               Age4_MA$beta),
                    se=c(Age0_MA$se,
                         Age1_MA$se,
                         Age2_MA$se,
                         Age3_MA$se,
                         Age4_MA$se),
                    ci.lb=c(Age0_MA$ci.lb,
                            Age1_MA$ci.lb,
                            Age2_MA$ci.lb,
                            Age3_MA$ci.lb,
                            Age4_MA$ci.lb),
                    ci.ub=c(Age0_MA$ci.ub,
                            Age1_MA$ci.ub,
                            Age2_MA$ci.ub,
                            Age3_MA$ci.ub,
                            Age4_MA$ci.ub),
                    SD_imputed=ifelse(simulateSD, "rtnorm sim", "within class var"),
                    Age_class=c("Age0","Age1","Age2", "Age3", "Age4"))
  
  
  ## Save and return results
  saveRDS(out_dat, "Data/weighted_MA.rds")
  return(out_dat)
}

