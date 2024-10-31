#' Reformat the carcass data to get hunting effort --> nr of succesfull hunters per year
#'
#' @param area_selection vector of study-area sub-area names to consider in the analyses: c("Inner", "BB",  "Tana")
#' @param carcass.dataset dataframe containing the carcass dataset downloaded from the COAT dataportal
#' @param shapefile.dir string. Directory of shapefile data with study areas sub-areas
#'
#' @return a dataframe containing numbers of succesfull hunters per year (regular and standardised)
#' @export
#'
#' @examples


reformatData_hunters <- function(area_selection,
                                      carcass.dataset, shapefile.dir){
  

# === 1 thing that I use later down ====
'%notin%' <- Negate('%in%')

#========= LOAD DATA ==============
## load data
allf <- carcass.dataset

shapefile <- st_read(paste(shapefile.dir, sep = "/"))
shapefile <- st_transform(shapefile, crs = 4326)

#check if loading carcass data worked
if(!exists("allf") || !length(allf)==42){
  stop("Carcass data not loaded properly or not 42 columns")
}

#check if loading shapefile data worked
if(!exists("shapefile")){
  stop("Shapefile not loaded properly")
}

#========= PREPARE DATA ==============
#we need hunting year as an actual number sometimes so we can add and substract
allf$start_hunting_year<-substr(allf$t_hunting_year, start = 1, stop = 4)
allf$start_hunting_year<- as.numeric(allf$start_hunting_year)

#Select Varanger only
fvar1 <- allf[is.na(allf$sn_region)== F & allf$sn_region == "varanger" , ]

#Remove the ones that were not shot
fvar1 <- fvar1[fvar1$v_hunter_id != "road_kill" & fvar1$v_hunter_id != "found_dead", ]

#Remove 2 foxes shot in nessenby_karlebotn (outside study area)
fvar1 <- fvar1[fvar1$v_place_name != "nessenby_karlebotn" | is.na(fvar1$v_place_name) == TRUE,]

#Assigning study area - sub area names with shapefile
fvar1   <- fvar1[is.na(fvar1$e_dd)==FALSE & is.na(fvar1$n_dd)==FALSE,] #only foxes with valid location
sf_fvar1<- st_as_sf(fvar1, coords = c("e_dd", "n_dd"), crs="+proj=longlat +datum=WGS84") # change dataset in to sf object with geometry (instead of lat/long)
sf_fvar1<- st_join(sf_fvar1, shapefile["name"], join = st_within) #overlap of dataset positions with shapefile of study area - sub areas, taking the name of the subarea
colnames(sf_fvar1)[colnames(sf_fvar1) == 'name'] <- 'sub_area'
coords  <- data.frame(st_coordinates(sf_fvar1)) #get back coordinates column from geometry
sf_fvar1$e_dd <- coords$X # getting back lat long columns
sf_fvar1$n_dd <- coords$Y
fvar1 <- st_drop_geometry(sf_fvar1) #remove the geometry again and get back to fvar1, #fvar1 now has shape file area names

#check if shapefile subarea assignment worked and subarea names match those given in script
if(sum(is.na(fvar1$sub_area)) > 0 || !identical(unique(fvar1$sub_area), c("Inner", "Tana",  "BB"))){  
  stop("Sub-area assignment with shapefile produced NA's, or sub-area names do not match those in script")
}

#Selecting for study area - sub areas
fvar1 <- subset(fvar1,(fvar1$sub_area %in% area_selection))

#======= HUNTING EFFORT ========
#nr of successfull hunters per year
hunter.data<- data.frame()

y <- unique(fvar1$start_hunting_year[!is.na(fvar1$start_hunting_year)])
for(i in y){
  newdata <- data.frame(
    start_hunting_year =i,
    NHunters = length(unique(fvar1$v_hunter_id[fvar1$start_hunting_year==i ]))
  )
  hunter.data<-rbind(hunter.data, newdata)
}

#set 2004 to NA because start project
hunter.data$NHunters[hunter.data$start_hunting_year==2004]<-NA

## Standardize and return
hunter.data$NHunters_std <- (hunter.data$NHunters - mean(hunter.data$NHunters, na.rm = T))/sd(hunter.data$NHunters, na.rm = T)

return(hunter.data)

}
