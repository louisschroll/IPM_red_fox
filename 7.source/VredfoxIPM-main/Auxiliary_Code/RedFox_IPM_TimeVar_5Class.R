#################################################
#### INTEGRATED POPULATION MODEL FOR RED FOX ####
#################################################

library(coda)
library(nimble)

mySeed <- 0


## Load data
IPM.data <- readRDS("RedFoxData_IPM.rds")
IPM.data <- IPM.data$total

#*******************#
# PRIOR INFORMATION #
#*******************#

## Prior information

# Early survival (denning period)
S0.mean <- 0.74 # AF IPM: 0.7428115
S0.sd <- 0.06 # AF IPM: 0.05983706

# Annual natural survival - Bristol (non-hunted)
#Snat.mean <- c(0.48, 0.54, 0.53, 0.51, 0.51)
#Snat.sd <- c(0.02, 0.03, 0.03, 0.03, 0.03)

# Annual natural survival - North Sweden (lightly hunted)
Snat.mean <- c(0.33, 0.71, 0.50, 0.59, 0.59)
Snat.sd <- c(0.02, 0.04, 0.05, 0.04, 0.04)

# Annual natural survival - Literature meta-analysis (all)
#Snat.mean <- c(0.38, 0.53, 0.57, 0.50, 0.52)
#Snat.sd <- c(0.10, 0.16, 0.18, 0.20, 0.21)

# Annual natural survival - Literature meta-analysis (non/lightly hunted)
#Snat.mean <- c(0.38, 0.53, 0.54, 0.46, 0.55)
#Snat.sd <- c(0.09, 0.10, 0.18, 0.25, 0.25)
# TODO: Check with Matt if we can re-do this meta-analysis properly

#***************#
# DATA OVERVIEW #
#***************#

## Data for Age-at-Harvest module 

IPM.data$C # Age-at-Harvest matrix (5 age classes, 15 seasons)
IPM.data$A # Oldest age class (5)
IPM.data$Tmax # Number of years in harvest carcass data
IPM.data$pData # Annual proportion of harvests with (complete) carcass information


## Data for Placental Scar module 

IPM.data$P1 # Vector of observed numbers of placental scars / fetuses (across individuals and years)
IPM.data$P1_age # Ages associated with elements of P1
IPM.data$P1_year # Years associated with elements of P1
IPM.data$X1 # Length of P1

IPM.data$P2 # Vector of observed presence/absence of placental scars / pregnancies (across individuals and years)
IPM.data$P2_age # Ages associated with elements of P2
IPM.data$P2_year # Years associated with elements of P2
IPM.data$X2 # Length of P2


## Arranging into data and constants for NIMBLE
RF.data <- list(C = IPM.data$C, P1 = IPM.data$P1, P2 = IPM.data$P2)

RF.constants <- list(
  A = IPM.data$A, Tmax = IPM.data$Tmax, pData = IPM.data$pData,
  P1_age = IPM.data$P1_age, P1_year = IPM.data$P1_year, X1 = IPM.data$X1,
  P2_age = IPM.data$P2_age, P2_year = IPM.data$P2_year, X2 = IPM.data$X2,
  maxPups = max(IPM.data$P1),
  S0.mean = S0.mean, S0.sd = S0.sd,
  
  # Survival prior
  Snat.mean = Snat.mean, Snat.sd = Snat.sd #,
  
  # Hoening model prior - fixed juv/ad ratio
  #mnat.logmean = mnat.logmean, mnat.logsd = mnat.logsd,
  #ratioJA.logmean = ratioJA.logmean
  
  # Hoening model prior - uncertainty in juv/ad ratio
  #mnat.logmean = mnat.logmean, mnat.logsd = mnat.logsd,
  #ratioJA.logmean = ratioJA.logmean, ratioJA.logsd = ratioJA.logsd
  
  #RodentAbundance = RodentAbundance,
  #RodentIndex2 = RodentIndex2,
  #RodentIndex3 = RodentIndex3
)


#********************#
# PARAMETER OVERVIEW #
#********************#

## Population Model
##-----------------

# N[a,t] = Number of individuals in age class a in year t
# R[a,t] = Number of recruits produced by age class a in year t
# L[a,t] = Number of offspring conceive by age class a in year t
# Imm[t] = Number of juvenile individuals immigrating into the popualtion in year t

# S[a,t] = Survival of age class a from year t to t+1

# rho[a,t] = number of placental scars (fetuses) of age class a in year t t
# Psi[a,t] = pregnancy rate of age class a in year t

# S0[t] = Denning survival of offspring (conception to emergence from the den) in year t

# Mu.Imm = average number of immigrants per year
# sigma.Imm = SD of among-year variation in number of immigrants


## Age-at-Harvest module
##----------------------

# h[a,t] = harvest rate of age class a (interval t to t+1)

# alpha[a,t] = proportion of deaths due to harvest in age class a (interval t to t+1)

# mH[a,t] = harvest mortality hazard rate of age class a (interval t to t+1)
# mO[a,t] = other (natural) cause mortality hazard rate of age class a (interval t to t+1)

# betaHE.mH = effect of harvest effort on mH
# {betaX.mO = effect of covariate X on mO (e.g. rodent/reindeer abundance)}

# {sigma.mH = SD of random year variation in mH}
# {sigma.mO = SD of random year variation in mO}


## Placental Scar module
##----------------------

# Mu.rho = baseline number of placental scars (fetuses)
# a.eff1 = linear age effect on the number of placental scars (fetuses)

# par.a = pregancy rate for oldest age class females
# par.b = slope for age effect on logit(pregnany rate)
# par.c = age when pregancy rate = par.a/2

# betaR.Psi = Effect of previous year rodent abundance on Psi

# {sigma.rho = SD of random year variation in rho}
# {sigma.Psi = SD of random year variation in rho}



#***************#
# NIMBLE - CODE #
#***************#

redfox.code <- nimbleCode({
  
  
  ##########################  
  #### POPULATION MODEL ####
  ##########################
  
  ### Likelihood (age classes: 1, 2, 3+)
  
  ## Survival
    
  for(t in 1:(Tmax-1)){ 
    
    # Age class 0 (index = 1): sum of local reproduction & immigrants
    N[1,t+1] <- sum(R[2:A,t+1]) + Imm[t+1]     
    
    # Age classes 1 to 3 (indeces = 2, 3, 4): age classes 0, 1, and 2 survivors    
    for(a in 1:(A-2)){
    	N[a+1,t+1] ~ dbin(S[a,t], N[a,t])
    }			
           
    # Age class 5+ (index = A = 5): age class 4 and 5+ survivors
    N[A,t+1] ~ dbin(S[A,t], N[A-1,t] + N[A,t])
  }
  
  ## Reproduction
  
  # Age class 0 (young of the year --> do not reproduce in year of birth)
  B[1,1:Tmax] <- 0
  L[1,1:Tmax] <- 0
  R[1,1:Tmax] <- 0
  
  # Age classes 1 to 3+    	    
  for(t in 1:Tmax){        				
    for(a in 2:A){
      
      # Breeding Population Size: Number of females that reproduce
      B[a,t] ~ dbin(Psi[a,t], N[a,t])
      
      # Litter Size (in utero): Number of pups produced by females of age class a
      L[a,t] ~ dpois(B[a,t]*rho[a,t]*0.5)
      
      # Number Recruits: Number of pups surviving to emerge from the den
      R[a,t] ~ dbin(S0[t], L[a,t])
    } 
  }
  
  #===============================================================================================
  
  
  
  ############################
  #### DERIVED QUANTITIES ####
  ############################
  
  for(t in 1:Tmax){
    N.tot[t] <- sum(N[1:A,t])
    R.tot[t] <- sum(R[1:A,t])		
    B.tot[t] <- sum(B[1:A,t])
  }
  
  #===============================================================================================

  
  
  ######################################
  #### WINTER AGE-AT-HARVEST MODULE ####
  ######################################
  
  ### Parameters:
  # N = number of individuals in a given age class at a given time
  # h = time-dependent probability of dying from hunting ([1] = adults, [2] = juveniles)
  
  ### Data:
  # C = age-at-harvest matrix
  # pData = annual proportion of harvests with (complete) carcass data
  
  ### Likelihood
  
  for(t in 1:Tmax){
    
    for(a in 1:A){
      C[a,t] ~ dbin(h[a,t]*pData[t], N[a,t])
    }
  }
  
  #===============================================================================================

  
  
  ###############################
  #### PLACENTAL SCAR MODULE ####
  ###############################
  
  ### Parameters:
  # rho = expected number of placental scars (fetuses)
  # Psi = pregnancy rate
  
  ## Data:
  # P1 = individual placental scar counts
  # P1_age = individual ages associated with P1
  # P1_year = year associated with P1

  # P2 = individual presence/absence of placental scars
  # P2_age = individual ages associated with P2
  # P2_year = year associated with P2
  
  
  ### Likelihood (litter size)
  
  for(x in 1:X1){
    P1[x] ~ dpois(rho[P1_age[x], P1_year[x]])
  }
  
  ### Likelihood (pregnancy rate)
  
  for(x in 1:X2){
    P2[x] ~ dbern(Psi[P2_age[x],P2_year[x]])
  }
  
  #===============================================================================================


  
  ################################
  #### PRIORS AND CONSTRAINTS ####
  ################################
  
  ## Survival and mortality
  
  for(t in 1:Tmax){ 
  	
  	# Harvest mortality hazard rate
  	log(mH[1:A,t]) <- log(Mu.mH[1:A]) + epsilon.mH[t]
  	#log(mH[1:A,t]) <- log(Mu.mH[a]) #+ betaHE.mH*HarvestEffort[t] + epsilon.mH[t]
  	
  	# Other (natural) mortality hazard rate
  	log(mO[1:A,t]) <- log(Mu.mO[1:A]) #+ betaX.mO*CovariateX[t] #+ epsilon.mO[t]
  	
  	# Survival probability
  	S[1:A,t] <- exp(-(mH[1:A,t] + mO[1:A,t]))
  	
  	# Proportion harvest mortality
  	alpha[1:A,t] <- mH[1:A,t]/(mH[1:A,t] + mO[1:A,t])
  	
    # Harvest rate
    h[1:A,t] <- (1-S[1:A,t])*alpha[1:A,t]
    
  }
  
  # Median harvest mortality hazard rates
  
  # Age-dependent
  for(a in 1:A){
  	Mu.mH[a] ~ dunif(0, 5) 
  }
  
  # Age-independent   
  #Mu.mH.all ~ dunif(0, 5) 
  #Mu.mH[1:A] <- Mu.mH.all
  
  # Median other (natural) cause mortality hazard rates
  #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE / HOENING MODEL CALCULATION

  # Using literature values on age-specific survival
  for(a in 1:A){
  	Mu.mO[a] <- -log(Mu.Snat[a])
  	Mu.Snat[a] ~ T(dnorm(Snat.mean[a], sd = Snat.sd[a]), 0, 1)   
  }
  
  # Using prior distributions calculated with Hoening model
  #Mu.mO.ad ~ dlnorm(mnat.logmean, logsd = mnat.logsd)
  #Mu.mO[2:5] <- Mu.mO.ad
  #Mu.mO[1] <- Mu.mO.ad*JuvAdRatio_mO #* NOTE: Can be provided as constant or distribution
  
  # JuvAdRatio <- exp(JAratio.logmean)
  # JuvAdRatio ~ dlnorm(JAratio.logmean, logsd = JAratio.logsd)
  
  # Covariate effects
  #betaHE.mH ~ dunif(0, 5) # Effect of harvest effort on mH
  #betaX.mO ~ dunif(-5, 5) # Effect of covariate X on mO (e.g. rodent/reindeer abundance)
  
 #---------------------------------------------------------------------------------------------

 
  ## Pregnancy rate

  for(t in 1:(Tmax+1)){
  	Psi[1,t] <- 0
  	
  	logit(Psi[2:A,t]) <- logit(Mu.Psi[2:A]) + epsilon.Psi[t]
    
    #logit(Psi[2:A,t]) <- logit(Mu.Psi[2:A]) + betaR.Psi*RodentAbundance[t] + epsilon.Psi[t]
  	#logit(Psi[2:A,t]) <- logit(Mu.Psi[2:A]) + betaRI2.Psi*RodentIndex2[t] + epsilon.Psi[t]
  	#logit(Psi[2:A,t]) <- logit(Mu.Psi[2:A]) + betaRI3.Psi[RodentIndex[t]+1] + epsilon.Psi[t]
  	
  }
  
  for(a in 2:A){	
  	Mu.Psi[a] ~ dunif(0, 1)
  }


  #betaR.Psi ~ dunif(-5, 5)
  #betaRI2.Psi ~ dunif(-5, 5)
  
  #betaRI3.Psi[1] <- 0
  
  #for(r in 2:3){
  #	betaRI3.Psi[r] ~ dunif(-5, 5)
  #}
  

    
 #---------------------------------------------------------------------------------------------

 
  ## Litter size

  for(t in 1:(Tmax+1)){
    rho[1,t] <- 0
  	log(rho[2:A,t]) <- log(Mu.rho[2:A]) + epsilon.rho[t]	
  }
  
  for(a in 2:A){
  	Mu.rho[a] ~ dunif(0, maxPups) # Baseline number of pups 
    #TODO:  ADJUST UPPER LIMIT 
  }
  
  
  #---------------------------------------------------------------------------------------------  

  
  ## Denning survival
  #* INFORMATIVE PRIOR REQUIRED: LITERATURE VALUE
  
  for(t in 1:Tmax){ 
    S0[t] <- mean.S0
    #S0[t] <- exp(-m0[t])
    #log(m0[t]) <- log(-log(mean.S0)) + epsilon.m0[t]
  }
  
  mean.S0 ~ T(dnorm(S0.mean, sd = S0.sd), 0, 1)  
  #---------------------------------------------------------------------------------------------
  
  
  ## Immigration
  
  Imm[1] <- 0 # (Immigration in the first year cannot be disentangled from reproduction)
  #ImmT[1] <- 0 
    
  for(t in 2:Tmax){
    Imm[t] ~ dcat(DU.prior[1:800]) 
    #Imm[t] ~ dpois(ImmT[t])
    #ImmT[t] ~ T(dnorm(Mu.Imm, sd = sigma.Imm), 0, 700)
  }
  
  #Mu.Imm ~ dunif(0, 400) #TODO: UPPER LIMIT NEEDS TO BE ADJUSTED
  #sigma.Imm ~ dunif(0, 500) #TODO: UPPER LIMIT NEEDS TO BE ADJUSTED
  
  # NOTE: Try constraining this again once the rest of the model is structured!
  
  #---------------------------------------------------------------------------------------------

  
  ## Initial population size (discrete uniform prior) 
  #TODO: UPPER LIMIT NEEDS TO BE ADJUSTED
  
  N[1:A,1] <- initN[1:A]
  
  for(a in 1:A){
    initN[a] ~ dcat(DU.prior[1:800]) 
  }
  
  DU.prior[1:800] <- 1/800

#---------------------------------------------------------------------------------------------


  ## Random year variation
  
  for(t in 1:Tmax){  
    epsilon.mH[t] ~ dnorm(0, sd = sigma.mH)
  # epsilon.mO[t] ~ dnorm(0, sd = sigma.mO)
    epsilon.Psi[t] ~ dnorm(0, sd = sigma.Psi)
    epsilon.rho[t] ~ dnorm(0, sd = sigma.rho) 
  # epsilon.m0[t] ~ dnorm(0, sd = sigma.m0)
  }
  
  sigma.mH ~ dunif(0, 5)
  #sigma.mO ~ dunif(0, 5) 
  sigma.Psi ~ dunif(0, 5)
  sigma.rho ~ dunif(0, 5)
  #sigma.m0 ~ dunif(0, 5) 
  

})


#****************#
# INITIAL VALUES #
#****************#

## Source function for setting initial values
source("RedFox_InitalValuesSim.R")

# Sampling initial values
#initVals <- RF.IPM.inits(IPM.data = RF.data, IPM.constants = RF.constants, 
#                         minN1 = c(600, 50, 50, 50, 50), 
#                         maxN1 = c(800, 400, 400, 400, 400), 
#                         minImm = 50, maxImm = 600)

initVals <- list(
  RF.IPM.inits(IPM.data = RF.data, IPM.constants = RF.constants, 
               minN1 = c(600, 50, 50, 50, 50), 
               maxN1 = c(800, 400, 400, 400, 400), 
               minImm = 50, maxImm = 600), 
  RF.IPM.inits(IPM.data = RF.data, IPM.constants = RF.constants, 
               minN1 = c(600, 50, 50, 50, 50), 
               maxN1 = c(800, 400, 400, 400, 400), 
               minImm = 50, maxImm = 600), 
  RF.IPM.inits(IPM.data = RF.data, IPM.constants = RF.constants, 
               minN1 = c(600, 50, 50, 50, 50), 
               maxN1 = c(800, 400, 400, 400, 400), 
               minImm = 50, maxImm = 600),
  RF.IPM.inits(IPM.data = RF.data, IPM.constants = RF.constants, 
               minN1 = c(600, 50, 50, 50, 50), 
               maxN1 = c(800, 400, 400, 400, 400), 
               minImm = 50, maxImm = 600))


#***********************#
# MCMC SETUP & TEST RUN #
#***********************#

## Setting parameters to monitor
IPM.params <- c(
"Mu.mH", "Mu.mO",
"Mu.Psi", "Mu.rho", "mean.S0",
"sigma.mH", "sigma.Psi", "sigma.rho",
#"betaR.Psi", "betaRI2.Psi", "betaRI3.Psi",
#"Mu.Imm", "sigma.Imm", 
"initN",
"N.tot", "B.tot", "R.tot", "N", "B", "L", "Psi", "rho", "mH", "Imm")

## MCMC settings
#niter <- 10
#nburnin <- 0
#nchains <- 4
#nthin <- 1

niter <- 30000
nburnin <- 5000
nchains <- 4
nthin <- 4

## Test run
t1 <- Sys.time()
RF.IPM.test <- nimbleMCMC(redfox.code, RF.constants, RF.data, 
                          initVals, monitors = IPM.params, 
                          niter = niter, nburnin = nburnin, 
                          nchains = nchains, thin = nthin, 
                          setSeed = mySeed, samplesAsCodaMCMC = TRUE)
Sys.time() - t1

saveRDS(RF.IPM.test, file = "RedFox_IPM_NSweden.rds")


#*******************#
# PRELIMINARY PLOTS #
#*******************#

library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape)

out.mat <- as.matrix(RF.IPM.test)

## Age-specific breeding probability
Psi.data <- data.frame(Estimate = c(out.mat[, paste0("Mu.Psi[", 2:5, "]")]), 
                       AgeClass = rep(c(1:3, "4+"), each = nrow(out.mat)))

Psi.plot <- ggplot(Psi.data, aes(x = AgeClass, y = Estimate, group = AgeClass)) + 
  geom_boxplot(aes(fill = AgeClass)) + 
  ylab("Breeding probability") + 
  scale_fill_viridis_d(option = "plasma", direction = -1) + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(), 
                     axis.ticks.x = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

## Age-specific fetus number
rho.data <- data.frame(Estimate = c(out.mat[, paste0("Mu.rho[", 2:5, "]")]), 
                       AgeClass = rep(c(1:3, "4+"), each = nrow(out.mat)))

rho.plot <- ggplot(rho.data, aes(x = AgeClass, y = Estimate, group = AgeClass)) + 
  geom_boxplot(aes(fill = AgeClass)) + 
  scale_fill_viridis_d(option = "plasma", direction = -1) + 
  ylab("Fetus number") + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())

pdf("Plots/Basic/Psi&rho.pdf", width = 12, height = 5)
grid.arrange(Psi.plot + theme(legend.position = "none"), rho.plot, nrow = 1, widths = c(0.8, 1))
dev.off()

## Age-specific harvest mortality
mH.data <- data.frame(Estimate = c(out.mat[, paste0("Mu.mH[", 2:5, "]")]), 
                       AgeClass = rep(c(1:3, "4+"), each = nrow(out.mat)))

mH.plot <- ggplot(mH.data, aes(x = AgeClass, y = Estimate, group = AgeClass)) + 
  geom_boxplot(aes(fill = AgeClass)) + 
  ylab("Harvest mortality") + 
  scale_fill_viridis_d(option = "plasma", direction = -1) + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())

## Age-specific natural mortality
mO.data <- data.frame(Estimate = c(out.mat[, paste0("Mu.mO[", 2:5, "]")]), 
                      AgeClass = rep(c(1:3, "4+"), each = nrow(out.mat)))

mO.plot <- ggplot(mO.data, aes(x = AgeClass, y = Estimate, group = AgeClass)) + 
  geom_boxplot(aes(fill = AgeClass)) + 
  ylab("Natural mortality") + 
  scale_fill_viridis_d(option = "plasma", direction = -1) + 
  theme_bw() + theme(axis.title.x = element_blank(), 
                     axis.text.x = element_blank(),
                     axis.ticks.x = element_blank(), 
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())


pdf("Plots/Basic/mH&mO.pdf", width = 12, height = 5)
grid.arrange(mH.plot + theme(legend.position = "none"), mO.plot, nrow = 1, widths = c(0.8, 1))
dev.off()


## Total population size over time
data <- melt(out.mat)
colnames(data) <- c("index", "parameter", "value")
data.sum <- ddply(data, .(parameter), summarise, median = median(value, na.rm = T), lCI_90 = quantile(value, probs = 0.05, na.rm = T), uCI_90 = quantile(value, probs = 0.95, na.rm = T), lCI_50 = quantile(value, probs = 0.25, na.rm = T), uCI_50 = quantile(value, probs = 0.75, na.rm = T))

popN <- paste("N.tot[", c(1:15), "]", sep = "")
data.Ntot <- subset(data.sum, parameter%in%popN)

data.Ntot$indexT <- c(10:15, 1:9)
data.Ntot$Year <- data.Ntot$indexT+2003
data.Ntot <- data.Ntot[order(data.Ntot$Year),]

pdf("Plots/Basic/Ntot_time.pdf", width = 6.5, height = 4)
ggplot(data.Ntot, aes(x = Year, y = median)) + geom_line(color = "#5C566B") + geom_ribbon(aes(ymin = lCI_90, ymax = uCI_90), alpha = 0.5, fill = "#5C566B") + ylab("Total population size") + scale_x_continuous(breaks = c(1997:2019)) + theme_bw() + theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
dev.off()
