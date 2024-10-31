#' Run fixed design transient life table response experiment (LTRE) for a pair of years
#'
#' @param paramSamples a list of lists containing posterior samples for all vital rates and
#' population-level quantities. The sublist "t" contains time-specific parameters
#' while the sublist "t_mean" contains time-average parameters. 
#' @param t_pair integer vector of length 2 specifying the indices of the two years
#' to compare in the analysis. 
#' @param Amax integer. Number of age classes. 
#' @param HazardRates logical. If TRUE (default), runs LTRE with respect to mortality 
#' hazard rates. If FALSE, runs LTRE with respect to survival probabilities. 
#' @param PopStructure logical. If TRUE (default), runs LTRE with respect to population 
#' proportions (n). If FALSE, runs LTRE with respect to age-specific population numbers (N).
#' @param save logical. If TRUE, saves the results as an .rds file. If FALSE (default) 
#' results are only returned in the console. 
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

runLTRE_fixedDesign <- function(paramSamples, t_pair, Amax, HazardRates = TRUE, PopStructure = TRUE, save = FALSE){

  #-------#
  # SETUP #
  #-------#
  
  ## Extract number of samples
  nosamples <- length(paramSamples$t_mean$lambda)

  ## List parameters to ignore
  dropParams <- c("N_tot", "B", "B_tot", "L", "L_tot", "R", "R_tot")
  
  if(HazardRates){
    dropParams <- c(dropParams, "S", "S0", "Ss")
  }else{
    dropParams <- c(dropParams, "mH", "mO", "m0", "mHs")
  }
  
  if(PopStructure){
    dropParams <- c(dropParams, "N")
  }else{
    dropParams <- c(dropParams, "n")
  }
  
  dropIdx <- which(names(paramSamples$t) %in% dropParams)
  
  paramList <- paramSamples$t[-dropIdx]
  
  ## Set up list of arrays for storing calculated LTRE contributions
  contList <- list()
  
  for(x in 1:length(paramList)){
    
    if(names(paramList[x]) == "lambda"){
      next
    }
    
    if(length(dim(paramList[[x]])) == 3){
      
      tempList <- list()
      for(a in 1:Amax){
        tempList <- c(tempList, list(as.numeric(rep(NA, nosamples))))
      }
      names(tempList) <- paste0(names(paramList)[x], "_", 1:Amax)
      
    }else{
      
      tempList <- list(as.numeric(rep(NA, nosamples)))
      names(tempList) <- names(paramList)[x]
    }
    
    contList <- c(contList, tempList)
  }
  
  contCount <- length(contList)
  
  contList$delta_lambda <- rep(NA, nosamples)
  contList$delta_loglambda <- rep(NA, nosamples)
  
  
  #------------------------------------------------#
  # CALCULATE SENSITIVITIES FOR RELEVANT YEAR PAIR #
  #------------------------------------------------#
  
  sensitivities <- calculateSensitivities(paramSamples = paramSamples,
                                          Amax = Amax,
                                          t_period = t_pair)
  
  sensList <- sensitivities$sensitivity$samples[-dropIdx]
  
  
  #--------------------------------------------------------#
  # CALCULATE LTRE CONTRIBUTIONS PER SAMPLE - FIXED DESIGN #
  #--------------------------------------------------------#
  
  for(i in 1:nosamples){
    
    ## Make lists of vital rates/population structure and of sensitivities
    param_list <- list()
    sens_list <- list()
    
    for(x in 1:length(paramList)){
      
      # Skip lambda
      if(names(paramList)[x] == "lambda"){
        next
      }
      
      # Set time interval based on parameter
      if(names(paramList)[x] %in% c("Psi", "rho", "S0", "m0")){
        tInt <- t_pair + 1
      }else{
        tInt <- t_pair
      }
      
      # Expand age class if required and list relevant parameter estimates
      if(length(dim(paramList[[x]])) == 3){
        
        tempList <- list()
        tempListS <- list()
        
        for(a in 1:Amax){
          tempList <- c(tempList, list(diff(paramList[[x]][i, a, tInt])))
          tempListS <- c(tempListS, list(sensList[[x]][i, a]))
        }
        names(tempList) <- paste0(names(paramList)[x], "_", 1:Amax)
        names(tempListS) <- paste0(names(paramList)[x], "_", 1:Amax)
      }else{
        
        tempList <- list(diff(paramList[[x]][i, tInt]))
        tempListS <- list(sensList[[x]][i])
        names(tempList) <- names(paramList)[x]
        names(tempListS) <- names(paramList)[x]
      }
      
      param_list <- c(param_list, tempList)
      sens_list <- c(sens_list, tempListS)
    }
    
    ## Convert parameter and sensitivity lists to vectors
    paramvec <- do.call(c, param_list)
    sensvec <- do.call(c, sens_list)
    
    ## Calculate and save change in lambda 
    contList$delta_lambda[i] <- diff(paramList$lambda[i, t_pair])
    contList$delta_loglambda[i] <- diff(log(paramList$lambda[i, t_pair]))

    ## Calculate demographic contributions
    # NOTE: Here we multiply sensitivities and parameter differences
    cont <- paramvec * sensvec
    names(cont) <- names(sensvec)
    
    ## Insert contributions into storage list
    for(x in 1:length(cont)){
      contList[[x]][i] <- cont[x]
    }
    
  }
  
  ## Restructure results list
  contList <- list(cont = contList[1:contCount],
                   other = list(delta_lambda = contList$delta_lambda,
                                delta_loglambda = contList$delta_loglambda)
                   )
  
  #-------------------#
  # SUMMARIZE RESULTS #
  #-------------------#
  
  ## Calculate summed contributions for age-specific parameters
  sumParams <- names(paramList)[which(!(names(paramList) %in% c("S0", "m0", "immR", "lambda")))] 
  
  for(x in 1:length(sumParams)){
    subList <- contList$cont[which(grepl(paste0(sumParams[x], "_"), names(contList$cont)))]
    
    contList$cont$newSum <- rowSums(dplyr::bind_rows(subList))
    names(contList$cont)[which(names(contList$cont) == "newSum")] <- paste0(sumParams[x], "_sum")
  }
  
  ## Calculate overall survival/harvest contributions
  if(HazardRates){
    contList$cont$mHtot_sum <- contList$cont$mH_sum + contList$cont$mHs_sum
  }else{
    contList$cont$Stot_sum <- contList$cont$S_sum + contList$cont$Ss_sum
  }
  
  ## Arrange results as dataframe
  contData <- melt(dplyr::bind_rows(contList$cont, .id = "column_label")) 
  contData <- contList$cont %>%
    dplyr::bind_rows() %>%
    tidyr::pivot_longer(cols = everything())
  
  contData <- cbind(contData, stringr::str_split_fixed(contData$name, "_", 2))
  colnames(contData) <- c("Variable", "Contribution", "Parameter", "AgeClass")
  
  contData$AgeClass[which(contData$AgeClass == "")] <- NA
  
  
  ## Make posterior summaries
  contData_summary <- contData %>%
    dplyr::group_by(Variable) %>%
    dplyr::summarise(lCI = quantile(Contribution, 0.025),
                     median = median(Contribution),
                     uCI = quantile(Contribution, 0.975))
  
  ## Collect and return results
  
  results <- list(contList = contList,
                  contData = contData,
                  contData_summary = contData_summary)
  
  if(save){
    saveRDS(results, file = paste0("RedFoxIPM_LTREresults_fixedDesign_t", t_pair[1], "-t", t_pair[2], ".rds"))
  }
  
  return(results)
}
