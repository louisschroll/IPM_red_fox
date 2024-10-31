#' Defines and stores information on prior for annual survival
#'
#' @param HoenigPrior logical. If TRUE, prior is based on Hoenig model. If FALSE, prior is based on literature values for annual survival.
#' @param sPriorSource string. Has to be provided if HoenigPrior = FALSE and specifies which literature source is used to inform parameters for prior on annual survival. 
#' Currently supported options: Bristol (not hunted), NSweden (North Sweden, lightly hunted), metaAll (meta-analysis of all populations), metaSub (meta-analysis of not/lightly hunted populations).
#'
#' @return a list containing information on prior for annual survival.
#' @export
#'
#' @examples

definePriorType_AnnSurv <- function(HoenigPrior, sPriorSource = NULL){
  
  ## Check consistency of information provided
  if(!HoenigPrior & !exists("sPriorSource")){
    stop("When using literature priors for annual survival (HoenigPrior = FALSE), information on prior source (sPriorSource) has to be provided.")
  }
  
  if(!HoenigPrior & !(sPriorSource %in% c("Bristol", "NSweden", "metaAll", "metaSub"))){
    stop("Invalid prior source information provided. The following options are currently supported: Bristol, NSweden, metaAll, metaSub.")
  }
  
  ## List prior type information
  PriorType <- list(Parameter = ifelse(HoenigPrior, "natMort", "annSurv"),
                    Source = ifelse(HoenigPrior, NA, sPriorSource))
  
  return(PriorType)
}