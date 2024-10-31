#' Adjusts age for reproduction data
#'
#' This function sets reproductive age by first adjusting sampling age and then 
#' adding 1 to match indexing in the model. 
#' For placental scars, reproductive age and sampling age are the same. 
#' For embryos and/or pregnancy signs, hypothetical reproductive age is one year 
#' higher than sampling age (as reproductive output of the individual would 
#' have resulted in recruitment to the population at the subsequent census).

#' @param data dataframe containing records of placental scars/embryos/pregnancy
#' signs.
#' @param ageCol integer. Column index where sampling age is stored in 'data'.
#' @param Amax integer. Number of age classes to consider in analyses.
#'
#' @return dataframe 'data' with an additional column 'age_adj'.
#' @export
#'
#' @examples

adjustAge_repData <- function(data, ageCol, Amax){
  
  data$age_adj <- NA
  for(x in 1:nrow(data)){
    
    if(data$type[x] == "pl.scars"){
      data$age_adj[x] <- ifelse(data[x, ageCol] < Amax-1, data[x, ageCol], Amax-1) + 1
    }
    
    if(data$type[x] %in% c("embryos", "pregnant", "combined")){
      data$age_adj[x] <- ifelse(data[x, ageCol] + 1 < Amax-1, data[x, ageCol] + 1, Amax-1) + 1
    }
  }
  
  return(data)
  
}