#' Run fixed design transient life table response experiment (LTRE) for all years
#'
#' @param paramSamples a list of lists containing posterior samples for all vital rates and
#' population-level quantities. The sublist "t" contains time-specific parameters
#' while the sublist "t_mean" contains time-average parameters. 
#' @param Amax integer. Number of age classes. 
#' @param Tmax integer. Number of years in the analysis.
#' @param HazardRates logical. If TRUE (default), runs LTRE with respect to mortality 
#' hazard rates. If FALSE, runs LTRE with respect to survival probabilities. 
#' @param PopStructure logical. If TRUE (default), runs LTRE with respect to population 
#' proportions (n). If FALSE, runs LTRE with respect to age-specific population numbers (N).
#'
#' @return  a list of lists containing results of the LTRE analysis. Object 'contList' 
#' collects posterior distributions for all parameters' LTRE contributions (sublist 'cont'),
#' as well as some auxiliary quantities (variances, co-variances, sublist 'other') as lists.
#' Object 'contData' is a dataframe consisting of posterior distributions for all parameters'
#' LTRE contributions. Object 'contData_summary' contains posterior summaries (medians and 
#' 95\% credible intervals) for the same information. 
#' @export
#'
#' @examples

runLTRE_fixedDesign_allYears <- function(paramSamples, Amax, Tmax, HazardRates = TRUE, PopStructure = TRUE){
  
  ## List all year pairs to run analyses for
  t_pairs <- cbind(2:(Tmax-2), 2:(Tmax-2) + 1)
  
  ## Set up results lists and dataframes
  contList <- list(cont = list(), other = list())
  contData <- data.frame()
  contData_summary <- data.frame()
  
  ## Run fixed design analyses for each pair of years
  for(t in 1:nrow(t_pairs)){
    LTRE_run <- runLTRE_fixedDesign(paramSamples = paramSamples, 
                                    t_pair = t_pairs[t,], 
                                    Amax = Amax,
                                    HazardRates = HazardRates, PopStructure = PopStructure)
    

    contList$cont[[t]] <- LTRE_run$contList$cont
    names(contList$cont)[t] <- paste0("t", t_pairs[t, 1], "t", t_pairs[t, 2])
    
    contList$other[[t]] <- LTRE_run$contList$other
    names(contList$other)[t] <- paste0("t", t_pairs[t, 1], "t", t_pairs[t, 2])
    
    contData_temp <- LTRE_run$contData
    contData_temp$t1 <- t_pairs[t, 1]
    contData_temp$t2 <- t_pairs[t, 2]
    contData <- rbind(contData, contData_temp)
    
    contData_summary_temp <- LTRE_run$contData_summary
    contData_summary_temp$t1 <- t_pairs[t, 1]
    contData_summary_temp$t2 <- t_pairs[t, 2]
    contData_summary <- rbind(contData_summary, contData_summary_temp)
    
  }
  
  ## Arrange all outputs in a new list
  results_allYears <- list(contList = contList, 
                           contData = contData,
                           contData_summary = contData_summary)
  
  
  ## Save and return
  saveRDS(results_allYears, file = "RedFoxIPM_LTREresults_fixedDesign.rds")
  return(results_allYears)
}
  
  