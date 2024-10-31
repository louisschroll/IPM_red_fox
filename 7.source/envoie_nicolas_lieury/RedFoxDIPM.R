### INTEGRATED POPULATION MODEL FOR VITAL RATE ESTIMATE OF RED FOX POPULATIONS

# Load R2WinBUGS package
library(R2WinBUGS)

# Specify directory where WinBUGS is stored
bugs.dir <- "c:/Users/etu-devillard/Documents/7. LOGICIELS/WinBUGS14"

# Specify working directory
# setwd("c:/Users/etu-devillard/Documents/1. THESE/data/renard/FDC35/Vital_rate/WinBugs")
setwd("~/IPM_red_fox/7.source")
#######################################################
## integrated distance sampling analysis in winbugs
#######################################################

fox <-read.table ("DOM.txt",h=T)       ## import line transect complete data
foxNA <- na.omit(fox)                  ## remove NA value for detection function estimation
# Data
nind  <-as.vector(tapply(fox$distance, fox$year, length))       ## number of total individual seen per transect per year
n <- as.vector(tapply(foxNA$distance,foxNA$year,length))       ## number of fox with distance value available for detection function estimation

vec1<-c(foxNA$distance[foxNA$year=="2002"],
        rep(NA, max(n)-length(foxNA$distance[foxNA$year=="2002"])))
vec2<-c(foxNA$distance[foxNA$year=="2003"],
        rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2003"])))
vec3<-c(foxNA$distance[foxNA$year=="2004"],
        rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2004"])))
vec4<-c(foxNA$distance[foxNA$year=="2005"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2005"])))
vec5<-c(foxNA$distance[foxNA$year=="2006"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2006"])))
vec6<-c(foxNA$distance[foxNA$year=="2007"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2007"])))
vec7<-c(foxNA$distance[foxNA$year=="2008"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2008"])))
vec8<-c(foxNA$distance[foxNA$year=="2009"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2009"])))
vec9<-c(foxNA$distance[foxNA$year=="2010"],rep(NA,max(n)-length(foxNA$distance[foxNA$year=="2010"])))

mat<-cbind(vec1,vec2,vec3,vec4,vec5,vec6,vec7,vec8,vec9)        ## matrix construction of distance distribution per year
C <- as.matrix(mat)
x <- C/1000                                                       ## translate in km
t.census <- ncol(x)                                                  ## number of year

cev <- as.vector( tapply(fox$length[fox$year=="2002"],fox$transect[fox$year=="2002"],mean))
cev1<-as.vector( tapply(fox$length[fox$year=="2003"],fox$transect[fox$year=="2003"],mean))
cev2<-as.vector( tapply(fox$length[fox$year=="2004"],fox$transect[fox$year=="2004"],mean))
cev3<-as.vector( tapply(fox$length[fox$year=="2005"],fox$transect[fox$year=="2005"],mean))
cev4<-as.vector( tapply(fox$length[fox$year=="2006"],fox$transect[fox$year=="2006"],mean))
cev5<-as.vector( tapply(fox$length[fox$year=="2007"],fox$transect[fox$year=="2007"],mean))
cev6<-as.vector( tapply(fox$length[fox$year=="2008"],fox$transect[fox$year=="2008"],mean))
cev7<-as.vector( tapply(fox$length[fox$year=="2009"],fox$transect[fox$year=="2009"],mean))
cev8<-as.vector( tapply(fox$length[fox$year=="2010"],fox$transect[fox$year=="2010"],mean))

L<-c(sum(cev),sum(cev1),sum(cev2),sum(cev3),sum(cev4),sum(cev5),sum(cev6),sum(cev7),sum(cev8))     ## vector construction of Total transect length per year

# Likelihood model from Gimenez et al 2009

sink("distmodelyear.txt")
cat("
model {
  # Prior
  theta ~ dgamma(15, 0.5)                                                        ## time invariant theta

  # likelihood
  for (t in 1:t.census) { 
    for (i in 1:n[t]) {
      zeros[i, t] <- 0
      zeros[i, t] ~ dpois(phi[i, t]) # likelihood is exp(-phi[i,t])
      # -log(likelihood)
      phi[i, t] <- -(log(2 * theta / 3.14) / 2 - theta * pow(x[i, t], 2) / 2)     ## estimation of detection function from available value
    } #i

  # derived parameters
  D[t] <- nind[t] * sqrt(2 * theta / 3.14) / (2 * L[t])                       ## density estimation per year with all number of encountered foxes
  N[t] <- D[t] * 366                                                              ## effective estimation per year given the surface of the area
 } #t
}
",fill = TRUE)
sink()

# Data (R 'list' format)

data = list(x = x, t.census = t.census, n = n, nind = nind, L = L)

# MCMC details
nb.iterations = 10000
nb.burnin = 1000
nt = 2

# Initial values (prior theta ~ Gamma[0.1,0.1])
init1 = list(theta = 0.1)
init2 = list(theta = 0.05)
init3 = list(theta = 0.01)

inits = list(init1,init2,init3)
nb.chains = length(inits)

# Parameters to be monitored
parameters <- c("theta", "D", "N")

library(jagsUI)
# MCMC simulations
res.sim <- jags(
  data,
  inits,
  parameters,
  "distmodelyear.txt",
  n.thin = nt,
  #bugs.directory = bugs.dir,
  n.chains = nb.chains,
  n.iter = nb.iterations,
  #debug = TRUE,
  n.burnin = nb.burnin
)

# Save the output
save(res.sim, file="DSdom.Rdata")

# Summarize results
print(res.sim,6)


# WRITE THE WinBUGS MODEL

sink("vulpes.bug")
cat("
##########################################################################
#  Integrated population model for the red fox population in Domagne France
#
#  Age structured, female-based model (2 age classes-1-year)
#  Pre-breeding census
#
#  Combination of:
#	- Distance sampling census (2002-2010): --> state-space model
#	- productivity for harvest data (2002-2006) --> Poisson regression
#	- age-at-harvest data from trapping and hunting(2002-2006) --> yearling and adult survival:
#		model by Udevitz & Gonan (2012)
#
##########################################################################
model
{
############################################################
# 1. Priors for the parameters
############################################################

    # Survival probabilities from the age-at-harvest data
for (t in 1:(nyear-1)) {
    Sy[t] ~ dunif(0, 0.8)
	}

for (t in 1:(nyear-1)) {
    Sa[t] ~ dunif(0.2, 1)
	}

    # Juvenile Survival probabilities
for (t in 1:(nyear-1)) {
    Sj[t] ~ dnorm(0.377,0.2)I(0,1)    # Sj prior from Devenish Nelson average of 8 population
	}

    # Productivity
for (t in 1:nyear) {
    muY[t] ~ dnorm(1.4, 0.2)I(0,2)   # constrained, because larger fecundity than exp(2) = 7 is not possible
    muA[t] ~ dnorm(1.6, 0.2)I(0,2)
  }

    # Proportion of breeding female
for (t in 1:nyear) {
    pY[t] ~ dunif(0,1)
    pA[t] ~ dunif(0.5,1)
  }

    # Immigration rate
for (t in 1:nyear) {
    im[t] ~ dgamma(0.1,0.5)
	}

    # Initial population sizes  # based on population count and literature age structure
    Ny[1] ~ dnorm(110,0.01)I(0,)
    Na[1] ~ dnorm(65,0.01)I(0,)

###################################################################
# 2. Derived and fixed parameters
###################################################################
    # Annual population growth rate
for (t in 2:nyear){
    lambda[t] <- Ntot[t] / Ntot[t-1]        # Ntot : total population size estimated from Distance Sampling
        }

    # Annual harvest rate
for (t in 1:nyear) {
    H[t] <- har[t]/ (Ntot[t]*2)                  # har : total removal quantity each year
	}

######################################################################
# 3. The Integrated Population Model
######################################################################

	####################################################
	# 3.1 Likelihood for the age-at-harvest data
	####################################################

	# binomial likelihood based on formula 6 and 13 in Udevitz & Gogan

  for (t in 2:nyear) {
		zy[t] <- Xit[1,t-1] * har[t]
		wy[t] <- har[t-1] * lmbda[t]
    totaly[t] <- zy[t] / wy[t]      # totaly[t] <- (Xit[1,t-1] * har[t]) / (har[t-1] * lmbda[t])

		za[t] <- sum(Xit[2:(nage-1),t-1]) * har[t]
		wa[t] <- har[t-1] * lmbda[t]
    totala[t] <- za[t] / wa[t]
	}

  for (t in 2:nyear) {
		X[2,t] ~ dbin(Sy[t-1], totaly[t])
		X[3,t] ~ dbin(Sa[t-1], totala[t])
	}

	####################################################
	# 3.2 Likelihood for productivity data
	####################################################

    ########################################
		# 3.2.1 Productivity data
		########################################
  for (t in 1:nyear) {

 			newborns[1,t] ~ dpois(rhoY[t])
			log(rhoY[t]) <- log(breeds[1,t]) + muY[t]

      newborns[2,t] ~ dpois(rhoA[t])
			log(rhoA[t]) <- log(breeds[2,t]) + muA[t]

      FecY[t] <- exp(muY[t])
      FecA[t] <- exp(muA[t])
  } #t



    ########################################
		# 3.2.1 Proportion of breeding female
		########################################
for (t in 1:nyear) {

 			breeds[1,t] ~ dbin(pY[t],females[1,t])
      breeds[2,t] ~ dbin(pA[t],females[1,t]) #prq mm indice (1) ?

} #t


	###################################################
	# 3.3 Likelihood for population survey data
	###################################################

		#############################################################
		# 3.3.1 System process:  2 age class matrix population model
		#############################################################

	for (t in 2:nyear){

		psi1[t] <- (pY[t-1]* FecY[t-1]* 0.5 * Ny[t-1])+ (pA[t-1]* FecA[t-1]* 0.5 * Na[t-1])
		JUV[t] ~ dpois(psi1[t])
       R[t] ~ dbin(Sj[t-1], JUV[t])		# number of local recruits

		psi2[t] <- (pY[t-1]* Ny[t-1] + pA[t-1] * Na[t-1]) * im[t-1]
		  IM[t] ~ dpois(psi2[t])		  # number of immigrants

		Ny[t] <- R[t] + IM[t]      # number of yearlings (age 1 class)

		Na1[t] ~ dbin(Sy[t-1], Ny[t-1])
    Na2[t] ~ dbin(Sa[t-1], Na[t-1])
      Na[t] <- Na1[t] + Na2[t]   # number of adults (age 2 class)

		} # t

		###############################################################
		# 3.3.2 Observation process : Spotlight Count Distance Sampling
		###############################################################

    for (t in 1:nyear){
    
     Ntot[t] <- Ny[t]+Na[t]
    census[t] ~ dpois(Ntot[t])             ## female effective estimation per year given the surface of the area
                         ## Total population size observe in population survey
       } #t

}  # End Model
",fill=TRUE)
sink()


# IMPORT DATA

# 1. Census data (2002-2010)

N <- round(res.sim$mean$N)
t.census <- length(N)

lam<-rep(NA,t.census)
for (t in 1:t.census) {
  lam[t + 1] <- N[t + 1] / N[t]              # Eq (2) from Udevitz 2012
}		# t

census <- c(N[1:5])                     # work only on 2002-2006 for the moment
census <- census / 2                    # work only on females
lmbda <- c(lam[1:5])

# 2. Age-at-death data (2002-2006)

tab <- read.table("DATA.txt", h = T, dec = ",")
tab1 <- subset(tab, tab$Mode != "D" &
                 tab$Mode != "R" &
                 tab$Age != "0")                   # remove biaised Digging out Data and juvenile class
DOM <- as.matrix(table(tab1$Age[tab1$GIC == "D"], tab1$Year[tab1$GIC == "D"]))        # work on Domagne population

Xit <- DOM[, 2:6]                                                              # keep only 2002-2006 for a better quality data
har <- colSums(Xit)
nyear <- ncol(Xit)
nage <- nrow(Xit)

X <- rbind(Xit[1:2, ], colSums(Xit[3:nage, ]))

# 3. Data on productivity (2002-2006)
tab2 <- subset(tab1, tab1$Info == "1" & tab1$GIC == "D")                        # keep Domagne female with information on repro
dom<- subset(tab2, tab2$Year!="1"& tab2$Year!="7")                              # keep only 2002-2006 for a better quality data

newborns<-cbind(tapply(dom$Count[dom$Year=="2"],dom$Age2.[dom$Year=="2"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="3"],dom$Age2.[dom$Year=="3"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="4"],dom$Age2.[dom$Year=="4"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="5"],dom$Age2.[dom$Year=="5"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="6"],dom$Age2.[dom$Year=="6"],sum,na.rm=T))
			## Nb of newborns estimated from embryos count and placental scars

breeds <- table(dom$Age2.[dom$Status == "R"], dom$Year[dom$Status == "R"])
## Matrix of number of breeding female per year and age class
			## assuming each breeding female have one brood

females <- table(dom$Age2.[dom$Info == "1"], dom$Year[dom$Info == "1"])
## Matrix of total number of female sampled per year and age class


# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification (time ~ 10 mn)
thin <- 20		# thinning
chain <- 3		# number of parallel chains
iter <- 100000	# number of iterations
burn <- 50000		# number of burn-in

# 2. Define parameters to be sampled
parameters <- c("Sj","Sy","Sa","im","pY","pA","FecY","FecA","Ntot","lambda","Ny","Na","R","IM","H")

# 3. Create random initial values to start the chains

inits <- function() {
  list(
    pY = c(runif(nyear, 0, 1)),
    pA = c(runif(nyear, 0.5, 1)),
    R = c(NA, rep(1, nyear - 1)),
    JUV = c(NA, rep(1, nyear - 1)),
    IM = c(NA, rep(1, nyear - 1)),
    muY = c(rnorm(nyear, 1.4, 0.4)),
    muA = c(rnorm(nyear, 1.6, 0.4)),
    Sj = c(runif((nyear - 1), 0, 0.6)),
    Sy = c(runif((nyear - 1), 0.2, 0.8)),
    Sa = c(runif((nyear - 1), 0.2, 1))
  )}

# 4. Combine all data into a list
vulpes <- list(
  census = census,
  lmbda = lmbda,
  nage = nage,
  nyear = nyear,
  newborns = newborns,
  breeds = breeds,
  females = females,
  Xit = Xit,
  X = X,
  har = har
)

# 5. Run the model from R
out <- bugs(
  vulpes,
  inits,
  parameters,
  "vulpes.bug",
  n.chains = chain,
  n.iter = iter,
  n.burn = burn,
  n.thin = thin,
  debug = T,
  bugs.directory = bugs.dir
)

print(out,6)

# 6. Save the output
save(out, file="vulpes.Rdata")


## FIGURES

#  Population size
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
year <- 2002:2006
for (i in 1:nyear){
   lower[i] <- quantile(out$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(out$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(out$mean$Ntot, census, lower), na.rm = T)
m2 <- max(c(out$mean$Ntot, census, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, nyear), ylab = "Population size", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:nyear, labels = year)
polygon(x = c(1:nyear, nyear:1), y = c(lower, upper[nyear:1]), col = "grey90", border = "grey90")
points(census, type = "l", col = "grey30", lwd = 2)
points(out$mean$Ntot, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 50, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "blue"), bty = "n", cex = 1)


lower <- upper <- numeric()
T <- nyear
for (t in 1:T){
   lower[t] <- quantile(out$sims.list$Ny[,t], 0.025)
   upper[t] <- quantile(out$sims.list$Ny[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = out$mean$Ny, x = (1:T)+0.5, xlim= c(1, 6), type = "b", pch = 1, ylim = c(0, 200), ylab = "Population size", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$Na[,t], 0.025)
   upper[t] <- quantile(out$sims.list$Na[,t], 0.975)}
points(y=out$mean$Na, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 30, legend = c("Yearlings", "Adults"), lty = c(1, 1),lwd = c(2, 2), pch = c(1, 16), bty = "n", cex = 1)

#  Survival rate
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- nyear-1
for (t in 1:T){
   lower[t] <- quantile(out$sims.list$Sy[,t], 0.025)
   upper[t] <- quantile(out$sims.list$Sy[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = out$mean$Sy, x = (1:T)+0.5, xlim= c(1, 5), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2006)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$Sa[,t], 0.025)
   upper[t] <- quantile(out$sims.list$Sa[,t], 0.975)}
points(y=out$mean$Sa, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$Sj[,t], 0.025)
   upper[t] <- quantile(out$sims.list$Sj[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = out$mean$Sj, x = (1:T)+0.5, xlim= c(1, 5), type = "b", pch = 17, ylim = c(0, 1), ylab = "Annual juvenile survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2006)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

#  Productivity rate : warning !! transform mu to have fec with fec=exp(mu)

par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- nyear
for (t in 1:T){
   lower[t] <- quantile(out$sims.list$FecY[,t], 0.025)
   upper[t] <- quantile(out$sims.list$FecY[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = out$mean$FecY, x = (1:T)+0.5, xlim= c(1, 6), type = "b", pch = 1, ylim = c(0, 7), ylab = "Annual Productivity (offspring per female)", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$FecA[,t], 0.025)
   upper[t] <- quantile(out$sims.list$FecA[,t], 0.975)}
points(y=out$mean$FecA, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$pY[,t], 0.025)
   upper[t] <- quantile(out$sims.list$pY[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = out$mean$pY, x = (1:T)+0.5, xlim= c(1, 6), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual Breeding Probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(out$sims.list$pA[,t], 0.025)
   upper[t] <- quantile(out$sims.list$pA[,t], 0.975)}
points(y=out$mean$pA, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 0.2, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")


#  immigration rate :

 par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)

lower <- upper <- numeric()
T <- nyear
for (t in 1:T){
   lower[t] <- quantile(out$sims.list$im[,t], 0.025)
   upper[t] <- quantile(out$sims.list$im[,t], 0.975)}
plot(y = out$mean$im, x = (1:T)+0.5, xlim = c(1, 6), type = "b", pch = 16, ylim = c(0, 1.1), ylab = "Immigration rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

#  growth rate :


lower <- upper <- numeric()
T <- nyear
for (t in 1:T-1){
   lower[t] <- quantile(out$sims.list$lambda[,t], 0.025)
   upper[t] <- quantile(out$sims.list$lambda[,t], 0.975)}
plot(y = lam[2:6], x = (1:T)+0.5, xlim = c(1,6), type = "b", pch = 16, ylim = c(0.5, 1.5), ylab = "Growth rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
points(y = out$mean$lambda, x = (1:4)+0.5, type = "b", pch = 1, cex = 1.5, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

legend(x = 1, y = 0.7, legend = c("Bugs", "Distance"), pch = c( 16, 1), bty = "n")

#  removal rate :

lower <- upper <- numeric()
T <- nyear 
for (t in 1:T){
   lower[t] <- quantile(out$sims.list$H[,t], 0.025)
   upper[t] <- quantile(out$sims.list$H[,t], 0.975)}
plot(y = out$mean$H, x = (1:T)+0.5, xlim = c(1,6), type = "b", pch = 16, ylim = c(0, 1), ylab = "Removal rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T), labels = 2002:2006)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)



# DESCRIPTIVE STATISTICS

# Survival  

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyear-1, ncol = 4)

for (i in 1:(nyear-1)){
   lambda.h[i] <- mean(out$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(out$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(out$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,1] <- mean(out$sims.list$Sj[,i])
   lower.h[i,1] <- quantile(out$sims.list$Sj[,i], 0.025)
   upper.h[i,1] <- quantile(out$sims.list$Sj[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,2] <- mean(out$sims.list$Sy[,i])
   lower.h[i,2] <- quantile(out$sims.list$Sy[,i], 0.025)
   upper.h[i,2] <- quantile(out$sims.list$Sy[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,3] <- mean(out$sims.list$Sa[,i])
   lower.h[i,3] <- quantile(out$sims.list$Sa[,i], 0.025)
   upper.h[i,3] <- quantile(out$sims.list$Sa[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,4] <- mean(out$sims.list$im[,i])
   lower.h[i,4] <- quantile(out$sims.list$im[,i], 0.025)
   upper.h[i,4] <- quantile(out$sims.list$im[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 7500)
for (i in 1:7497){
   correl.h[i,1] <- cor(out$sims.list$lambda[i,], out$sims.list$Sj[i,], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(out$sims.list$lambda[i,], out$sims.list$Sy[i,], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(out$sims.list$lambda[i,], out$sims.list$Sa[i,], use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(out$sims.list$lambda[i,], out$sims.list$im[i,1:4], use = "pairwise.complete.obs")
   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,2], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,3], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,4], c(0.05, 0.5, 0.95), na.rm = TRUE)


# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,2], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,3], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,4], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/7500
sum(correl.h[!is.na(correl.h[,2]),2]>0)/7500
sum(correl.h[!is.na(correl.h[,3]),3]>0)/7500
sum(correl.h[!is.na(correl.h[,4]),4]>0)/7500


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.4, y = 1.6, "r = 0.85 (-0.58,0.56, 0.96)", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 1.5, "P(r>0) = 0.81", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Yearling survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.4, y = 1.6, "r = 0.62 (-0.87,0.07, 0.90)", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 1.5, "P(r>0) = 0.53", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.4, y = 1.6, "r = -0.86 (-0.96,-0.50, 0.55)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.4, y = 1.5, "P(r>0) = 0.21", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Immigration rate", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.4, y = 1.6, "r = 0.83 (-0.86,0.12, 0.92)", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 1.5, "P(r>0) = 0.56", pos = 4, font = 3, cex = 0.8)

# Productivity  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyear-1, ncol = 4)

for (i in 1:(nyear-1)){
   lambda.h[i] <- mean(out$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(out$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(out$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,1] <- mean(out$sims.list$pY[,i])
   lower.h[i,1] <- quantile(out$sims.list$pY[,i], 0.025)
   upper.h[i,1] <- quantile(out$sims.list$pY[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,2] <- mean(out$sims.list$pA[,i])
   lower.h[i,2] <- quantile(out$sims.list$pA[,i], 0.025)
   upper.h[i,2] <- quantile(out$sims.list$pA[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,3] <- mean(out$sims.list$FecY[,i])
   lower.h[i,3] <- quantile(out$sims.list$FecY[,i], 0.025)
   upper.h[i,3] <- quantile(out$sims.list$FecY[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,4] <- mean(out$sims.list$FecA[,i])
   lower.h[i,4] <- quantile(out$sims.list$FecA[,i], 0.025)
   upper.h[i,4] <- quantile(out$sims.list$FecA[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 7497)
for (i in 1:7497){
   correl.h[i,1] <- cor(out$sims.list$lambda[i,], out$sims.list$pY[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(out$sims.list$lambda[i,], out$sims.list$pA[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(out$sims.list$lambda[i,], out$sims.list$FecY[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(out$sims.list$lambda[i,], out$sims.list$FecA[i,1:4], use = "pairwise.complete.obs")
   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,2], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,3], c(0.05, 0.5, 0.95), na.rm = TRUE)
quantile(correl.h[,4], c(0.05, 0.5, 0.95), na.rm = TRUE)


# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,2], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,3], na.rm = TRUE)
m$x[which(m$y==max(m$y))]

m <- density(correl.h[,4], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/7497
sum(correl.h[!is.na(correl.h[,2]),2]>0)/7497
sum(correl.h[!is.na(correl.h[,3]),3]>0)/7497
sum(correl.h[!is.na(correl.h[,4]),4]>0)/7497


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Yearling breeding probability", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = -0.87 (-0.96,-0.54, 0.48)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.18", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult breeding probability", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.59 (-0.76,0.25, 0.90)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.64", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Yearling productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.11 (-0.88,-0.04, 0.82)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.47", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.71 (-0.82,0.14, 0.90)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.58", pos = 4, font = 3, cex = 0.8)

# Removal intensity

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyear-1, ncol = 4)

for (i in 1:(nyear-1)){
   lambda.h[i] <- mean(out$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(out$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(out$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,1] <- mean(out$sims.list$H[,i])
   lower.h[i,1] <- quantile(out$sims.list$H[,i], 0.025)
   upper.h[i,1] <- quantile(out$sims.list$H[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 7500)
for (i in 1:7497){
   correl.h[i,1] <- cor(out$sims.list$lambda[i,], out$sims.list$H[i,1:4], use = "pairwise.complete.obs")

   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)



# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/7500


# Plot Fig. 11-8

linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.4), ylab = "Population growth rate", xlab = "Removal", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")

text(x = 0.1, y = 0.7, "r = 0.93 (-0.01,0.73,0.98)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 0.65, "P(r>0) = 0.94", pos = 4, font = 3, cex = 0.8)
