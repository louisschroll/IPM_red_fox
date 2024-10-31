rm(list=ls())

## Read in the carcass examination files. This could probably be modified to Stijn's 
## version to retrieve this from the data base

in.dir<-"C:/Users/deh000/Box Sync/KOAT/Arctic fox/rodrev/coatdb" 
dataset<-"V_redfox_carcass_examination"

setwd(paste(in.dir, "carcass_examination", sep = "/"))
filenames <- dir()

mylist<-c()
for (i in 1:length(filenames)){
  mylist[[i]]<-read.table(paste(in.dir,"carcass_examination", filenames[i], sep = "/"), header=T, sep = ";")
}

allf <- do.call(rbind, mylist)  # combine all files

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Read the genetic data and merge them with the carcass examination file

setwd("C:/Users/deh000/Box Sync/KOAT/Arctic fox/rodrev/Genetikk/analyses2023")
gen <- read.delim("REDFOX_data1.txt", stringsAsFactors = F)

summary(gen$ind3 %in% allf$v_individual_id) # 42 are false
gen$ind3[gen$ind3 %in% allf$v_individual_id == F]  
# 1379 needs to be removed! Fox without hunting year
gen <- gen[gen$ind3 != "1379", ]

library(plyr)
gen$ind4 <- gen$ind3
gen$ind4 <- revalue(gen$ind4, c("2009_50juv"="50", "2010_111" = "111", "2010_112"="112", "2010_119"="119", "2010_120"="120"))

summary(gen$ind4 %in% allf$v_individual_id) # 37 are false
gen$ind4[gen$ind4 %in% allf$v_individual_id == F] # These are OK, Nordkynn and sør Varanger

names(gen)[26] <- "v_individual_id"
head(gen)
gena <- merge(gen, allf, by="v_individual_id", all.x = T) # 699 lines

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Tranform genetic data to genepop file, one column per locus and three digits format
# and assemble with necessary red fox data

geno <- gena[ , 3:26]
# replace 2 digits genotypes by 3 digits genotypes
for (i in 1:dim(geno)[1]){
  for(j in 1:dim(geno)[2]){
    if(nchar(geno[i, j]) < 3) {geno[i, j] <- paste(0, geno[i, j], sep="") } 
  }}
# merge to one column per locus
geno1 <- matrix(nrow=dim(gena)[1], ncol = 12)
for (i in 1:dim(geno1)[2]) {
  geno1[ ,i] <- paste(geno[ , i*2-1],geno[ , i*2], sep="")
  }

geno1 <- data.frame(geno1)

GEN <- cbind(id = gena[ , 1], geno1, gena[, c(27, 35, 36, 46, 47, 49)])
locnames <- c("AHT133", "C08.618", "CPH2", "FH2001", "FH2054", "FH2328", "CPH7", "CPH11", "CPH18", "CXX468", "FH2848", "REN54P1")
names(GEN)[2:13] <- locnames
head(GEN)
GEN$reg <- GEN$sn_region
GEN$reg[is.na(GEN$reg)] <- "other"
GEN$reg[GEN$reg == "finnmarksvidda"] <- "other"
table(GEN$reg, useNA="ifany") # 194 other, 505 varanger

# function to write genepop file

Write.genepop <- function(var, other=NA, ye, vrs, nbpop){
  filename <- paste(vrs, "_", ye,".dat", sep="")
  cat(paste("Genepop file", vrs, "year", ye, "\n"), file=filename)
  for(i in 1:length(locnames)){
    cat(paste(locnames[i], "\n"), file=filename, append = T)
  }
  cat(paste("POP", "\n"), file=filename, append = T)
  for (i in 1:dim(var)[1]){
    cat(paste(var[i, 1], ",", var[i, 2], var[i, 3], var[i, 4], var[i, 5], var[i, 6], var[i, 7], 
              var[i, 8], var[i, 9], var[i, 10], var[i, 11], var[i, 12], var[i, 12], "\n"), file=filename, append = T)
  }
  if(nbpop == 2){
  cat(paste("POP", "\n"), file=filename, append = T)
  for (i in 1:dim(other)[1]){
    cat(paste(other[i, 1], ",", other[i, 2], other[i, 3], other[i, 4], other[i, 5], other[i, 6], other[i, 7], 
              other[i, 8], other[i, 9], other[i, 10], other[i, 11], other[i, 12], other[i, 12], "\n"), file=filename, append = T)
}}}


###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Make input files for Geneclass to identify immigrants 

huntyears <- sort(unique(GEN$t_hunting_year[is.na(GEN$t_hunting_year)==F]))
nbyears <- length(huntyears)
source <- GEN[GEN$reg == "other", ]
GEN$yearhunt <- as.numeric(substring(GEN$t_hunting_year, 1, 4))

## Two versions were retained to use in the IPM

## Version 1
# Moving base line based on the foxes shot every year. Foxes from one year as assigned to two possible 
# reference populations: Varanger (all foxes shot in the current year and all foxes shot in the two previous 
# years) or the potential source populations (all foxes samples available from Iesjavri, Neiden and Norkynn, 
# independent of the year in which they were shot).

yearhunts <- c(2004:2014)

# year 1
varanger1 <- GEN[GEN$reg == "varanger" & is.na(GEN$yearhunt)==F & GEN$yearhunt == yearhunts[1], ] # all the ones shot that year
Write.genepop(rbind(varanger1), source, ye = huntyears[1], vrs = "Vimm2", nbpop=2)
rm(varanger1)

# year 2
varanger1 <- GEN[GEN$reg == "varanger" & is.na(GEN$yearhunt)==F & GEN$yearhunt %in% yearhunts[1:2], ] # all the ones shot in the previous year
Write.genepop(rbind(varanger1), source, ye = huntyears[2], vrs = "Vimm2", nbpop=2)
rm(varanger1)

for(i in 3:nbyears){
  varanger1 <- GEN[GEN$reg == "varanger" & is.na(GEN$yearhunt)==F & GEN$yearhunt %in% yearhunts[(i-2):(i)], ] # all the ones shot in these two year year
  Write.genepop(rbind(varanger1), source, ye = huntyears[i], vrs = "Vimm2", nbpop=2)
  rm(varanger1)
}

# Run the files with the option detect first generation immigrants with 1000 simulations in Geneclass 2.0.


## Version 2
# Moving baseline assuming all foxes immigrated in the first year of their life. Here all foxes are 
# classified into cohorts based on their birth year. Then they were assigned to two possible reference 
# populations: Varanger (all foxes from Varanger present in the population in the focal year and the 
# previous year - those shot then a and those present then based on the fact that they were born earlier 
# and shot later) and the potential source populations (all foxes samples available from Iesjavri, Neiden 
# and Norkynn, independent of the year in which they were shot).

# We assume that the foxes without age were born in their first year
GEN$v_year_birth[is.na(GEN$v_year_birth) & GEN$reg == "varanger"] <- GEN$yearhunt[is.na(GEN$v_year_birth) & GEN$reg == "varanger"]

# The first focal year will be 2000 to have a sufficient sample size
birthyears <- c(2000:2014)
nbbirth <- length(birthyears)

for(i in 1:nbbirth){
  varanger1 <- GEN[GEN$reg == "varanger" & is.na(GEN$yearhunt)==F & is.na(GEN$v_year_birth) == F &
                     GEN$v_year_birth <= birthyears[i] & GEN$yearhunt >= birthyears[i], ] 
  # all the ones born before the focal year, including the focal year
  # and all the ones shot the focal season or later

  Write.genepop(varanger1, source, ye = birthyears[i], vrs = "Vimm62", nbpop=2)
  rm(varanger1)
}

# Run the files with the option detect first generation immigrants with 1000 simulations in Geneclass 2.0.

###~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Assemble the output file

# read the output files for version 1
Readresults1 <- function(vrs, yrs) {
  resfiles <- vector("character", length(yrs))
  for (i in 1:length(resfiles)) {
    resfiles[i] <- paste("result-", vrs, "_", yrs[i], ".csv", sep="")}
  
  reslist<-c()
  for (i in 1:length(resfiles)){
    temp <- read.table(resfiles[i], header=F,skip = 15, sep = ";")
    temp$id <- substring(temp$V1, 2)
    reslist[[i]] <- temp[temp$id %in% GEN$id[GEN$t_hunting_year == yrs[i] & GEN$reg == "varanger"], ]
  }
  
  resfile<-do.call(rbind, reslist)  # combine all files
  resfile$id <- substring(resfile$V1, 2)
  return(data.frame(id = resfile$id, prob = resfile$V4, llvar = resfile$V5, llother = resfile$V6, logratio = resfile$V3))
}

# read the output files for version 2
Readresults2 <- function(vrs, yrs) {
  resfiles <- vector("character", length(yrs))
  for (i in 1:length(resfiles)) {
    resfiles[i] <- paste("result-", vrs, "_", yrs[i], ".csv", sep="")}
  
  reslist<-c()
  for (i in 1:length(resfiles)){
    temp <- read.table(resfiles[i], header=F,skip = 15, sep = ";")
    temp$id <- substring(temp$V1, 2)
    reslist[[i]] <- temp[temp$id %in% GEN$id[GEN$v_year_birth == yrs[i] & GEN$reg == "varanger"], ]
  }
  
  resfile<-do.call(rbind, reslist)  # combine all files
  resfile$id <- substring(resfile$V1, 2)
  return(data.frame(id = resfile$id, prob = resfile$V4, llvar = resfile$V5, llother = resfile$V6, logratio = resfile$V3))
}


resvers1 <- Readresults1(vrs = "Vimm2", yrs = huntyears)
names(resvers1) <- c("id", "pvarV1", "llvar", "llother", "logratio")
head(resvers1)
dim(resvers1) # 505 - all OK!

resvers2 <- Readresults2(vrs = "Vimm62", yrs = birthyears)
head(resvers2)
dim(resvers2) # 487 - the ones born earlier are lacking
names(resvers2) <- c("id", "pvarV2", "llvar2", "llother2", "logratio2")

# Assemble a file with the resulting estimates

immigration <- merge(GEN[ , c("id", "t_hunting_year", "v_age", "v_year_birth", "v_sex")], resvers1[ ,c("id","pvarV1")], by = "id", all.y = T)
immigration <- merge(immigration, resvers2[ ,c("id","pvarV2")], by = "id", all.x = T)

names(immigration)[1] <- "v_individual_id"
write.table(immigration, file = "RedFox_genetics_immigrant_probabilities.txt", sep="\t", quote=F, row.names = F)

