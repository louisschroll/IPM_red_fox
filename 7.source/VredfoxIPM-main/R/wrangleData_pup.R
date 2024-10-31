#' Loading and reformatting opportunistic data on pups per den
#'
#' @param datapath character string with the path to the file containing the
#' data. 
#' @param minYear integer. The first year to consider in the analyses. 
#'
#' @return a list containing pup counts (NoPups), associated years (NoPups_year)
#' and the number of observations (X3). 
#' @export
#'
#' @examples
#' 
wrangleData_pup <- function(datapath, minYear){
  
  ## Load data file
  pupData <- read.csv(datapath, sep = ";")
  
  ## Reformat data
  NoPups <- pupData$nr.pups
  NoPups_year <- pupData$year - minYear + 1
  
  ## Return data
  return(list(NoPups = NoPups,
              NoPups_year = NoPups_year,
              X3 = length(NoPups)))
  
}