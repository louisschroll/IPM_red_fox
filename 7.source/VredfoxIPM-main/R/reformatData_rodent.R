#' Makes yearly rodent abundance variables
#'
#' @param rodent.dataset dataframe containing the carcass dataset downloaded from the COAT dataportal
#'
#' @param minYear integer. First year to consider in analyses.
#'
#' @return a list containing rodent abundance data as continuous variables (cont),
#' and categorical variables with two (cat2). Data are provided for winter 
#' (fall + spring) in Varanger (.wintvar) and for fall for the larger area (.fallstor).
#' Continuous data are provided as total sums of individuals across all species
#' and as sums weighed by species (voles vs. lemmings, .stsp). 
#' Note that the time indices are shifted forward to represent that reproduction
#' is a function of past rodent abundance.
#' @export
#'
#' @examples


reformatData_rodent <- function(rodent.dataset, minYear) {

#========= LOAD DATA ==============
allrod <- rodent.dataset
#================ Select rodents from storskala areas ==============================

storskala <- reshape2::dcast(allrod, sn_locality + sn_site + t_year + t_season ~ v_species, value.var = "v_abundance", fun.aggregate = sum) 

#remove row where season is NA
storskala <- storskala[!is.na(storskala$t_season) ,]

# remove aves sorex and ranidae and neo_fod (keep "cricetidae","lem_lem", "mic_oec", "myo_ruf", "myo_rut")
storskala<- subset(storskala, select = -c(aves,neo_fod, ranidae, sor_ara, sor_cae, sor_min, sor_sp))

#rename some columns
names(storskala)[1:4] <- c("region", "plot", "year", "season")
names(storskala)[names(storskala) == 'cricetidae'] <- 'rodsp'
names(storskala)[names(storskala) == 'mic_oec'] <- 'Moec'
names(storskala)[names(storskala) == 'myo_ruf'] <- 'Mruf'
names(storskala)[names(storskala) == 'myo_rut'] <- 'Mrut'
names(storskala)[names(storskala) == 'lem_lem'] <- 'Llem'

#================================================================
storskala$session <- as.factor(paste(storskala$year, storskala$season, sep=""))

# aggregate by region, plot, session, season, year
stor <- aggregate(storskala$Llem, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T) 
sum(storskala$Llem, na.rm=T)
sum(stor$x) # OK except for one with lacking plot number
names(stor) <- c("reg", "plot", "session", "season", "year", "Llem") 
stor$Moec <- aggregate(storskala$Moec, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$Mruf <- aggregate(storskala$Mruf, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$Mrut <- aggregate(storskala$Mrut, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x 
stor$rodsp <- aggregate(storskala$rodsp, by=list(storskala$region, storskala$plot, storskala$session, storskala$season, storskala$year), sum, na.rm=T)$x
stor$vole <- stor$Mruf + stor$Moec + stor$Mrut + stor$rodsp
stor$tot <- stor$vole + stor$Llem

# Aggregate annual means into three regions; Varanger; Nordkyn; Ifjord-gaissene
stor$foxreg <- NA
stor$foxreg[stor$reg %in% c("komagdalen",  "stjernevann", "vestre_jakobselv")] <- "varanger"
stor$foxreg[stor$reg %in% c("nordkynn",  "bekkarfjord")] <- "nordkynn"
stor$foxreg[stor$reg == "ifjordfjellet"] <- "ifjordfjellet"

#--- continuous rodent variables WINTER (VARANGER only) ---
agvar <-  aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot) ~ year + season + foxreg, stor, mean)  #the mean nr of rodents per plot, for each year and season
agvar <- agvar[agvar$foxreg=="varanger",] # we only use varanger

agvar$start_hunting_year <- ifelse(agvar$season == "spring", agvar$year-1, agvar$year) # -1 year from spring, so we can add fall to form a hunting year representative rodent variable instead of adding spring and fall from the same year

agvar_fall <- agvar[agvar$season=="fall",]
agvar_fall$st.tot <- (agvar_fall$tot-mean(agvar_fall$tot))/sd(agvar_fall$tot) #vole and lemming actual numbers together
agvar_fall$st.lem <- (agvar_fall$Llem-mean(agvar_fall$Llem))/sd(agvar_fall$Llem)  #standardise lemming
agvar_fall$st.vole <- (agvar_fall$vole-mean(agvar_fall$vole))/sd(agvar_fall$vole) #standardise vole
agvar_fall$st.lemvole <- (agvar_fall$st.lem+ agvar_fall$st.vole) # vole and lemming together after standardised

agvar_spring <- agvar[agvar$season=="spring",]
agvar_spring$st.tot <- (agvar_spring$tot-mean(agvar_spring$tot))/sd(agvar_spring$tot)#vole and lemming actual numbers together
agvar_spring$st.lem <- (agvar_spring$Llem-mean(agvar_spring$Llem))/sd(agvar_spring$Llem) #standardise lemming
agvar_spring$st.vole <- (agvar_spring$vole-mean(agvar_spring$vole))/sd(agvar_spring$vole)#standardise vole
agvar_spring$st.lemvole <- (agvar_spring$st.lem+ agvar_spring$st.vole) # vole and lemming together after standardised

agvar_winter <- rbind(agvar_fall, agvar_spring)
agvar_winter <- aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot, st.tot, st.lem, st.vole, st.lemvole) ~ start_hunting_year , agvar_winter, mean)

#remove data for start_hunting year 2003 (only spring 2004, no rodents from 2003 fall)
agvar_winter[ agvar_winter$start_hunting_year==2003, -1] <-NA

#if spring and fall have dissimilar "start_hunting_year" in the final row, this "start_hunting_year" should also be NA (because spring rodents not collected yet for that year)
if (agvar_spring[nrow(agvar_spring), 11] != agvar_fall[nrow(agvar_fall), 11]) {
 agvar_winter[nrow(agvar_winter), -1] <-NA    
}
  

#--- categorical rodent variable WINTER (VARANGER only)---
#if we split the continuous abundance into 2 -> high and low, at threshold -0.34, 9 years are high and 9 are low, 
# and the categories are the same whether we split based on st.tot or st.lemvole (as long as you choose -0.34 as your threshold)
agvar_winter$cat2 <- ifelse(agvar_winter$st.tot>-0.34, 1,0)


#--- continuous rodent variables FALL (Varanger + IFJORD + NORDKYNN), because immigration ---
agstor <-  aggregate(cbind(Llem, Moec, Mruf, Mrut, rodsp, vole, tot) ~ year + season , stor, mean)  #the mean nr of rodents per plot, for each year and season

agstor$start_hunting_year <- ifelse(agstor$season == "spring", agstor$year-1, agstor$year) # -1 year from spring, so we can add fall to form a hunting year representative rodent variable instead of adding spring and fall from the same year
agstor_fall <- agstor[agstor$season=="fall",]

agstor_fall$st.tot <- (agstor_fall$tot-mean(agstor_fall$tot))/sd(agstor_fall$tot) #vole and lemming actual numbers together
agstor_fall$st.lem <- (agstor_fall$Llem-mean(agstor_fall$Llem))/sd(agstor_fall$Llem)  #standardise lemming
agstor_fall$st.vole <- (agstor_fall$vole-mean(agstor_fall$vole))/sd(agstor_fall$vole) #standardise vole
agstor_fall$st.lemvole <- (agstor_fall$st.lem+ agstor_fall$st.vole) # vole and lemming together after standardised

#--- categorical rodent variable FALL (Varanger + IFJORD + NORDKYNN), because immigration---
#if we split the continuous abundance into 2 -> high and low, at threshold 0.25, 8 years are high and 11 are low, 
# and the categories are the same whether we split based on st.tot or st.lemvole (as long as you choose 0.25 as your threshold)
agstor_fall$cat2 <- ifelse(agstor_fall$st.tot>0.25, 1,0)


## List and return
return(list(cont.wintvar          =   agvar_winter$st.tot,      #winter varanger continuous, only standardised for seasons
            cont.wintvar.stsp     =   agvar_winter$st.lemvole,  #winter varanger continuous, standardised for seasons and species
            cat2.wintvar          =   agvar_winter$cat2 + 1,    #winter varanger 2 categories
            
            cont.fallstor         =   agstor_fall$st.tot,     #fall storskala continuous
            cont.fallstor.stsp    =   agstor_fall$st.lemvole, #fall storskala continuous, standardised for species
            cat2.fallstor         =   agstor_fall$cat2 + 1,   #fall storskala 2 factors
            
            YearInfo.wint         =   paste0("fall ", agvar_winter$start_hunting_year, " - spring ", agvar_winter$start_hunting_year + 1),
            YearInfo.fall         =   paste0("fall ", agstor_fall$start_hunting_year)))

}


