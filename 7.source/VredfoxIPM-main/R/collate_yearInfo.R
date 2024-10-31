#' Collate information on year indexing
#'
#' @param minYear integer. First year to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#'
#' @return a dataframe containing year indices as used in the model with 
#' corresponding reproduction years and winter/summer harvest seasons. 
#' This dataframe is not required to run the analyses but the information is
#' useful for understanding and interpreting model and results. 
#' @export
#'
#' @examples

collate_yearInfo <- function(minYear, Tmax){

  ## Set up dataframe with year index information
  YearInfo <- data.frame(
    index = 1:Tmax, # year index as used in model
    ReproductionYear = (1:Tmax) + minYear - 1, # reproduction year
    WinterHarvestSeason = paste0("Oct ", (1:Tmax) + minYear - 1, " - May ", (1:Tmax) + minYear), # winter harvest season
    SummerHarvestSeason = paste0("Jul ", (1:Tmax) + minYear - 1, " - Sep ", (1:Tmax) + minYear - 1)) 

  return(YearInfo)
}

