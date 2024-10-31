#' Makes yearly reindeer carcass abundance variable
#'
#' This function relies on two datafiles "KadaverOstFinnmark-22.rds" and
#' "KadaverOstF-GamvikLebesby-22.rds". These files originate from Rovbase and
#' have subsequently been converted to .rds in an R session where the locale was
#' set to "Norwegian Bokmål_Norway".
#'  
#' @param minYear integer. First year to consider in analyses.
#' @param Tmax integer. Number of years to consider in analysis
#' @param reinCov.VarTana logical. If TRUE uses varanger (+Tana) municipalities to count reindeer carcasses, if false uses all of eastern finnmark
#' 
#'
#' @return a list containing reindeer carcass abundance data as scaled continuous variables,
#' for each winter
#' @export
#'
#' @examples

reformatData_reindeer <- function(minYear, Tmax, reinCov.VarTana ){

#---------Load data-----------
#Carcasses found in Varanger municipalities (+Tana) + Sor varanger
carc_VTS <- readRDS("Data/KadaverOstFinnmark-22.rds")
#Carcasses found in Gamvik and Lebesby
carc_GL <- readRDS("Data/KadaverOstF-GamvikLebesby-22.rds")

#Combine to carcasses found in whole area:
#East finnmark in this case, but If we get the file directly from Rovbase, it will have all municipalities.
carc <- rbind(carc_VTS, carc_GL)

#-----Select and format data ----------
#select only reindeer
carc  <- carc[carc$Skadetype == "Rein", ]

#Format dates
carc$dato <- as.Date(carc$Funnetdato, format="%d.%m.%Y")
carc$year <- as.numeric(format(carc$dato, format="%Y"))
carc$month <- as.numeric(format(carc$dato, format="%m"))

#select years above 2000
carc <- carc[carc$year > 2000, ]
#select november-may months
carc <- carc[carc$month %in% c(1:6, 11:12),]

#assemble per winter, winter year of december, as with fox hunting year.
carc$winter <- carc$year
carc$winter[carc$month %in% c(1:6)] <- carc$year[carc$month %in% c(1:6)] - 1

#aggregate sum of carcasses for all of eastern Finnmark
carc$nb <- 1
cwinter <- aggregate(carc$nb, by=list(carc$winter), sum)
names(cwinter) <- c("winter", "EFinnmark")

# If we get the file directly from Rovbase, it will have all municipalities from the start. So this is to sort
#Varanger (+Tana), Sor-Varanger
carc_VTS <- carc[carc$Kommune %in% c("Deatnu - Tana (N)", "BerlevC%g (N)", "Soer-Varanger (N)",
                                  "BC%tsfjord (N)", "VardC8 (N)", "VadsC8 (N)", "UnjC!rga - Nesseby (N)"), ]
#Varanger (+Tana)
carc_VT <- carc[carc$Kommune %in% c("Deatnu - Tana (N)", "BerlevC%g (N)", 
                                  "BC%tsfjord (N)", "VardC8 (N)", "VadsC8 (N)", "UnjC!rga - Nesseby (N)"), ]

carc_VTS$nb <- 1
carc_VT$nb <- 1
cwinter$SVaranger <- aggregate(carc_VTS$nb, by=list(carc_VTS$winter), sum)$x
cwinter$Varanger <- aggregate(carc_VT$nb, by=list(carc_VT$winter), sum)$x

## Define year range
years <- minYear:(minYear + Tmax - 1)

## Assemble yearly data + Toggle to select Varanger (+Tana) (smallest scale) or whole of eastern finnmark (largest scale)
RDcarcass <- rep(NA, Tmax+1)

if(reinCov.VarTana){
  for(t in 1:Tmax){
  year <- t + minYear - 1
    if(year %in% cwinter$winter){
    RDcarcass[t] <- cwinter$Varanger[which(cwinter$winter == year)]
    }
  }
} else{
  for(t in 1:Tmax){
    year <- t + minYear - 1
    if(year %in% cwinter$winter){
      RDcarcass[t] <- cwinter$EFinnmark[which(cwinter$winter == year)]
    }
  }
}

## Arrange relevant data in a list
reindeer.data <- list(RDcarcass = as.vector(scale(RDcarcass)))

## Return data
return(reindeer.data)

}