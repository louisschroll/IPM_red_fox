#' Prepare reproduction data
#'
#' This function relies on a helper function adjustAge_repData. 
#' 
#' @param P1.datafile # Placental scar/embryo count data from: P1.datafile <- carcassData$P1var
#' @param P2.datafile # Presence of placental scars/embryos/pregnancy signs from: P2.datafile <- carcassData$P2var
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#'
#' @return a list containing two formatted dataframes 'P1' (counts) and 'P2'
#' (presences/absences).
#' @export
#'
#' @examples

wrangleData_rep <- function(P1.datafile, P2.datafile, Amax, minYear){

  ## Adjust age indices
  P1.data <- adjustAge_repData(data = P1.datafile, ageCol = 2, Amax = Amax)
  P2.data <- adjustAge_repData(data = P2.datafile, ageCol = 2, Amax = Amax)
  
  ## Remove data for age 0 (index 1) individuals
  age0.rep_P1 <- which(P1.data$age_adj == 1)
  if(length(age0.rep_P1) > 0){
    P1.data <- P1.data[!P1.data$age_adj ==1, ]  #I had to change this part (used to be: P1.data <- P1.data[-age0.rep_P1, ] ). This is because if there are no age 1 with placenta, age0.rep_P1 = integer(0), P1.data[-integer(0)] will result in 0 rows for P1.data
    warning(paste0(length(age0.rep_P1), " instance(s) of age 0 individuals with placental scars removed from data (P1)"))
  }

  age0.rep_P2 <- subset(P2.data, age_adj == 1 & P2 == 1)
  P2.data <- P2.data[!(P2.data$age_adj == 1), ] #I had to change this part (used to be: P2.data <- P2.data[-which(P2.data$age_adj == 1), ] ). This is because if there are no age 1 with embryo, age0.rep_P2 = integer(0), P2.data[-integer(0)] will result in 0 rows for P2.data
  if(nrow(age0.rep_P2) > 0){
    warning(paste0(nrow(age0.rep_P2), " instance(s) of age 0 individuals with placental scars removed from data (P2)")) #Actually I already removed 0 year olds from the placental scar record in the data.reformatting script, so I suppose this step is double
  }

  ## Add year index
  P1.data$RepYearIndex <- P1.data$repryear - minYear + 1
  P2.data$RepYearIndex <- P2.data$repryear - minYear + 1

  ## List and return data
  return(list(P1 = P1.data, P2 = P2.data))
  
}
