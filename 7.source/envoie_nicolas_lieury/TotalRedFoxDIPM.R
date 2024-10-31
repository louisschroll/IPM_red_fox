### INTEGRATED POPULATION MODEL FOR VITAL RATE ESTIMATE OF RED FOX POPULATIONS

# Load R2WinBUGS package
library(R2WinBUGS)

# Specify directory where WinBUGS is stored
bugs.dir <- "c:/Users/etu-devillard/Documents/7. LOGICIELS/WinBUGS14"

# Specify working directory
setwd("c:/Users/etu-devillard/Documents/1. THESE/data/renard/FDC35/Vital_rate/WinBugs")

# WRITE THE WinBUGS MODEL

sink("vulpestot.bug")
cat("
##########################################################################
#  Integrated population model for the red fox population in Domagne France
#
#  Age structured, female-based model (2 age classes-1-year )
#  Pre-breeding census
#
#  Combination of:
#	- Distance sampling census (2002-2010): --> state-space model
#	- productivity for harvest data (2002-2006) --> Poisson regression
#	- age-at-harvest data from trapping and hunting (2002-2006) --> yearling and adult survival:
#		model by Udevitz & Gonan (2012)
#
##########################################################################
model
{
############################################################
# 1. Priors for the parameters
############################################################

    # Survival probabilities from the age-at-harvest data
for (t in 1:(t.census-1)) {
    Sy[t] ~ dunif(0,0.8)
	}

for (t in 1:(t.census-1)) {
    Sa[t] ~ dunif(0.2,1)
	}

    # Juvenile Survival probabilities
for (t in 1:(t.census-1)) {
    Sj[t] ~ dnorm(0.377,0.2)I(0,1)    # Sj prior from Devenish Nelson average of 8 population
	}

    # Productivity
for (t in 1:t.census) {
    muY[t] ~ dnorm(1.4,0.2)I(0,2)   # constrained, because larger fecundity than exp(2) = 7 is not possible
    muA[t] ~ dnorm(1.6,0.2)I(0,2)
  }

    # Proportion of breeding female
for (t in 1:t.census) {
    pY[t] ~ dunif(0,1)
    pA[t] ~ dunif(0.5,1)
  }

    # Immigration rate
for (t in 1:t.census) {
    im[t] ~ dgamma(0.1,0.5)
	}

    # Initial population sizes
    Ny[1] ~ dnorm(110,0.01)I(0,)
    Na[1] ~ dnorm(65,0.01)I(0,)

###################################################################
# 2. Derived and fixed parameters
###################################################################
    # Annual population growth rate
for (t in 2:t.census){
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
    totaly[t] <- zy[t] / wy[t]

		za[t] <- sum(Xit[2:(nage-1),t-1]) * har[t]
		wa[t] <- har[t-1] * lmbda[t]
    totala[t] <- za[t] / wa[t]
		}

for (t in 2:nyear) {

		X[2,t] ~ dbin(Sy[t-1],totaly[t])

		X[3,t] ~ dbin(Sa[t-1],totala[t])

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

} #t

    ########################################
		# 3.2.1 Proportion of breeding female
		########################################
for (t in 1:nyear) {

 			breeds[1,t] ~ dbin(pY[t],females[1,t])
      breeds[2,t] ~ dbin(pA[t],females[1,t])

} #t


	###################################################
	# 3.3 Likelihood for population survey data
	###################################################

		#############################################################
		# 3.3.1 System process:  2 age class matrix population model
		#############################################################

	for (t in 2:t.census){

		psi1[t] <- (pY[t-1] * exp(muY[t-1])* 0.5 * Ny[t-1])+ (pA[t-1]* exp(muA[t-1])* 0.5 * Na[t-1])
		JUV[t] ~ dpois(psi1[t])
    R[t] ~ dbin(Sj[t-1], JUV[t])		# number of local recruits

		psi2[t] <- (pY[t-1] * Ny[t-1] + pA[t-1] * Na[t-1]) * im[t-1]
		IM[t] ~ dpois(psi2[t])		  # number of immigrants

		Ny[t] <- R[t] + IM[t]      # number of yearlings (age 1 class)

		Na1[t] ~ dbin(Sy[t-1], Ny[t-1])
    Na2[t] ~ dbin(Sa[t-1], Na[t-1])
    Na[t] <- Na1[t] + Na2[t]   # number of adults (age 2 class)

		} # t

		###############################################################
		# 3.3.2 Observation process : Spotlight Count Distance Sampling
		###############################################################

    for (t in 1:t.census){

     Ntot[t] <- Ny[t]+Na[t]
    census[t] ~ dpois(Ntot[t])             ## female effective estimation per year given the surface of the area
                         ## Total population size observe in population survey
       } #t

}  # End Model
",fill=TRUE)
sink()

###################
#### DOMAGNE
###################

# IMPORT DATA

# 1. Census data (2002-2010)

domD <- c(0.851, 0.973, 1.032, 1.043, 1.075, 1.35, 1.048, 1.355, 1.268)
domVD<-c( 0.007,0.01,0.017,0.0165,0.011,0.031,0.012,0.019,0.015	)

N<-round(domD*366 )
t.census<-length(N)

lam<-rep(NA,t.census)
for (t in 1:(t.census-1)) {
  lam[t+1]<- N[t+1]/N[t]              # Eq (2) from Udevitz 2012
  }		# t

census<-c(N/2)                      # work only on females
lmbda<-c(lam[1:5])

# 2. Age-at-death data (2002-2006)

tab<-read.table("DATA.txt",h=T,dec=",")
tab1<-subset(tab, tab$Mode!="D"& tab$Mode!="R"& tab$Age!="0")                   # remove biaised Digging out Data and juvenile class
DOM<-as.matrix(table(tab1$Age[tab1$GIC=="D"], tab1$Year[tab1$GIC=="D"]))        # work on Domagne population

	Xit <- DOM[,2:6]                                                              # keep only 2002-2006 for a better quality data
	har<- colSums(Xit)
	nyear<- ncol(Xit)
  nage<-nrow(Xit)

  X <- rbind(Xit[1:2,],colSums(Xit[3:nage,]))

# 3. Data on productivity (2002-2006)
tab2<-subset(tab1, tab1$Info=="1"& tab1$GIC=="D")                               # keep Domagn? female with information on repro
dom<- subset(tab2, tab2$Year!="1"& tab2$Year!="7")                              # keep only 2002-2006 for a better quality data

newborns<-cbind(tapply(dom$Count[dom$Year=="2"],dom$Age2.[dom$Year=="2"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="3"],dom$Age2.[dom$Year=="3"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="4"],dom$Age2.[dom$Year=="4"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="5"],dom$Age2.[dom$Year=="5"],sum,na.rm=T),
            +tapply(dom$Count[dom$Year=="6"],dom$Age2.[dom$Year=="6"],sum,na.rm=T))
			## Nb of newborns estimated from embryos count and placental scars

breeds<-table(dom$Age2.[dom$Status=="R"],dom$Year[dom$Status=="R"])
			## Matrix of number of breeding female per year and age class
			## assuming each breeding female have one brood

females<-table(dom$Age2.[dom$Info=="1"],dom$Year[dom$Info=="1"])
			## Matrix of total number of female sampled per year and age class


# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification   (time ~ 30mn)
thin <- 20		# thinning
chain <- 3		# number of parallel chains
iter <- 100000	# number of iterations
burn <- 50000		# number of burn-in

# 2. Define parameters to be sampled
parameters <- c("Sj","Sy","Sa","im","pY","pA","muY","muA","Ntot","lambda","Ny","Na","R","IM","H")

# 3. Create random initial values to start the chains

initstot <- function(){list(pY=c(runif(t.census,0,1)), pA=c(runif(t.census,0.5,1)), R=c(NA,rep(1,t.census-1)), JUV=c(NA,rep(1,t.census-1)), IM=c(NA,rep(1,t.census-1)),  muY=c(rnorm(t.census,1.4,0.4)), muA=c(rnorm(t.census,1.6,0.4)), Sj=c(runif((t.census-1),0,0.6)), Sy=c(runif((t.census-1),0.2,0.8)), Sa=c(runif((t.census-1),0.2,1)))}

# 4. Combine all data into a list
vulpestot <- list(census=census,t.census=t.census,lmbda=lmbda, nage=nage,nyear=nyear,newborns=newborns,breeds=breeds,females=females,Xit=Xit,X=X,har=har)

# 5. Run the model from R
outot <- bugs(vulpestot, initstot, parameters, "vulpestot.bug", n.chains=chain, n.iter=iter, n.burn=burn, n.thin=thin, debug=T, bugs.directory = bugs.dir)

print(outot,6)

# 6. Save the output
save(outot, file="vulpestot.Rdata")



## FIGURES

#  Population size
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
year <- 2002:2010
for (i in 1:t.census){
   lower[i] <- quantile(outot$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(outot$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(outot$mean$Ntot, census, lower), na.rm = T)
m2 <- max(c(outot$mean$Ntot, census, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, t.census), ylab = "Population size", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:t.census, labels = year)
polygon(x = c(1:t.census, t.census:1), y = c(lower, upper[t.census:1]), col = "grey90", border = "grey90")
points(census, type = "l", col = "grey30", lwd = 2)
points(outot$mean$Ntot, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 50, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "blue"), bty = "n", cex = 1)


lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$Ny[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$Ny[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = outot$mean$Ny, x = (1:T)+0.5, xlim= c(1, 10), type = "b", pch = 1, ylim = c(0, 200), ylab = "Population size", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$Na[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$Na[,t], 0.975)}
points(y=outot$mean$Na, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 30, legend = c("Yearlings", "Adults"), lty = c(1, 1),lwd = c(2, 2), pch = c(1, 16), bty = "n", cex = 1)

#  Survival rate
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census-1
for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$Sy[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$Sy[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = outot$mean$Sy, x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$Sa[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$Sa[,t], 0.975)}
points(y=outot$mean$Sa, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$Sj[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$Sj[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = outot$mean$Sj, x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 17, ylim = c(0, 1), ylab = "Annual juvenile survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

#  Productivity rate : warning !! transform mu to have fec with fec=exp(mu)

FecY<-exp(outot$mean$muY)
FecA<-exp(outot$mean$muA)

par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- exp(quantile(outot$sims.list$muY[,t], 0.025))
   upper[t] <- exp(quantile(outot$sims.list$muY[,t], 0.975))}
par(mgp=c(3.8,1,0))
plot(y = exp(outot$mean$muY), x = (1:T)+0.5, xlim= c(1, 10), type = "b", pch = 1, ylim = c(0, 7), ylab = "Annual Productivity (offspring per female)", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- exp(quantile(outot$sims.list$muA[,t], 0.025))
   upper[t] <- exp(quantile(outot$sims.list$muA[,t], 0.975))}
points(y=exp(outot$mean$muA), x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$pY[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$pY[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = outot$mean$pY, x = (1:T)+0.5, xlim= c(1, 10), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual Breeding Probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$pA[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$pA[,t], 0.975)}
points(y=outot$mean$pA, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 0.2, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")


#  immigration rate :

 par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)

lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$im[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$im[,t], 0.975)}
plot(y = outot$mean$im, x = (1:T)+0.5, xlim = c(1, 10), type = "b", pch = 16, ylim = c(0, 1.1), ylab = "Immigration rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

#  growth rate :


lower <- upper <- numeric()
T <- t.census -1
for (t in 1:T){
   lower[t] <- quantile(outot$sims.list$lambda[,t], 0.025)
   upper[t] <- quantile(outot$sims.list$lambda[,t], 0.975)}
plot(y = outot$mean$lambda, x = (1:T)+0.5, xlim = c(1,9), type = "b", pch = 16, ylim = c(0.5, 1.5), ylab = "Growth rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
points(y = lam[2:9], x = (1:T)+0.5, type = "b", pch = 1, cex = 1.5, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2002:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

legend(x = 1, y = 0.7, legend = c("Bugs", "Distance"), pch = c( 16, 1), bty = "n")



# DESCRIPTIVE STATISTICS

# Survival  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(outot$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(outot$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(outot$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(outot$sims.list$Sj[,i])
   lower.h[i,1] <- quantile(outot$sims.list$Sj[,i], 0.025)
   upper.h[i,1] <- quantile(outot$sims.list$Sj[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(outot$sims.list$Sy[,i])
   lower.h[i,2] <- quantile(outot$sims.list$Sy[,i], 0.025)
   upper.h[i,2] <- quantile(outot$sims.list$Sy[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(outot$sims.list$Sa[,i])
   lower.h[i,3] <- quantile(outot$sims.list$Sa[,i], 0.025)
   upper.h[i,3] <- quantile(outot$sims.list$Sa[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(outot$sims.list$im[,i])
   lower.h[i,4] <- quantile(outot$sims.list$im[,i], 0.025)
   upper.h[i,4] <- quantile(outot$sims.list$im[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 59976)
for (i in 1:59976){
   correl.h[i,1] <- cor(outot$sims.list$lambda[i,], outot$sims.list$Sj[i,], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(outot$sims.list$lambda[i,], outot$sims.list$Sy[i,], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(outot$sims.list$lambda[i,], outot$sims.list$Sa[i,], use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(outot$sims.list$lambda[i,], outot$sims.list$im[i,1:8], use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/59976
sum(correl.h[!is.na(correl.h[,2]),2]>0)/59976
sum(correl.h[!is.na(correl.h[,3]),3]>0)/59976
sum(correl.h[!is.na(correl.h[,4]),4]>0)/59976


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.48 (-0.35,0.35, 0.80)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.80", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Yearling survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = -0.13 (-0.70,-0.07, 0.66)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.44", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = -0.17 (-0.71,-0.06, 0.62)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.44", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Immigration rate", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.20 (-0.55, 0.11,0.66)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.60", pos = 4, font = 3, cex = 0.8)

# Productivity  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(outot$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(outot$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(outot$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(outot$sims.list$pY[,i])
   lower.h[i,1] <- quantile(outot$sims.list$pY[,i], 0.025)
   upper.h[i,1] <- quantile(outot$sims.list$pY[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(outot$sims.list$pA[,i])
   lower.h[i,2] <- quantile(outot$sims.list$pA[,i], 0.025)
   upper.h[i,2] <- quantile(outot$sims.list$pA[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(exp(outot$sims.list$muY[,i]))
   lower.h[i,3] <- quantile(exp(outot$sims.list$muY[,i]), 0.025)
   upper.h[i,3] <- quantile(exp(outot$sims.list$muY[,i]), 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(exp(outot$sims.list$muA[,i]))
   lower.h[i,4] <- quantile(exp(outot$sims.list$muA[,i]), 0.025)
   upper.h[i,4] <- quantile(exp(outot$sims.list$muA[,i]), 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 59976)
for (i in 1:59976){
   correl.h[i,1] <- cor(outot$sims.list$lambda[i,], outot$sims.list$pY[i,1:8], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(outot$sims.list$lambda[i,], outot$sims.list$pA[i,1:8], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(outot$sims.list$lambda[i,], exp(outot$sims.list$muY[i,1:8]), use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(outot$sims.list$lambda[i,], exp(outot$sims.list$muA[i,1:8]), use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/59976
sum(correl.h[!is.na(correl.h[,2]),2]>0)/59976
sum(correl.h[!is.na(correl.h[,3]),3]>0)/59976
sum(correl.h[!is.na(correl.h[,4]),4]>0)/59976


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Yearling breeding probability", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.34 (-0.58,0.18, 0.72)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.65", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult breeding probability", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.34 (-0.29,0.26, 0.71)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.77", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Yearling productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.22 (-0.66,0.07, 0.66)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.56", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.33 (-0.46,0.27, 0.75)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.74", pos = 4, font = 3, cex = 0.8)





##########################################
## VENDELAIS POPULATION
##########################################

 
# WRITE THE WinBUGS MODEL

sink("vulpestot.bug")
cat("
##########################################################################
#  Integrated population model for the red fox population in Vendelais France
#
#  Age structured, female-based model (2 age classes-1-year )
#  Pre-breeding census
#
#  Combination of:
#	- Distance sampling census (2003-2010): --> state-space model
#	- productivity for harvest data(2003-2006) --> Poisson regression
#	- age-at-harvest data from trapping and hunting(2003-2006) --> yearling and adult survival:
#		model by Udevitz & Gonan (2012)
#
##########################################################################
model
{
############################################################
# 1. Priors for the parameters
############################################################

    # Survival probabilities from the age-at-harvest data
for (t in 1:(t.census-1)) {
    Sy[t] ~ dunif(0,0.8)
	}

for (t in 1:(t.census-1)) {
    Sa[t] ~ dunif(0.2,1)
	}

    # Juvenile Survival probabilities
for (t in 1:(t.census-1)) {
    Sj[t] ~ dnorm(0.377,0.2)I(0,1)    # Sj prior from Devenish Nelson average of 8 population
	}

    # Productivity
for (t in 1:t.census) {
    muY[t] ~ dnorm(1.4,0.2)I(0,2)   # constrained, because larger fecundity than exp(2) = 7 is not possible
    muA[t] ~ dnorm(1.6,0.2)I(0,2)
  }

    # Proportion of breeding female
for (t in 1:t.census) {
    pY[t] ~ dunif(0,1)
    pA[t] ~ dunif(0.5,1)
  }

    # Immigration rate
for (t in 1:t.census) {
    im[t] ~ dgamma(0.1,0.5)
	}

    # Initial population sizes
    Ny[1] ~ dnorm(65,0.01)I(0,)
    Na[1] ~ dnorm(45,0.01)I(0,)

###################################################################
# 2. Derived and fixed parameters
###################################################################
    # Annual population growth rate
for (t in 2:t.census){
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
    totaly[t] <- zy[t] / wy[t]

		za[t] <- sum(Xit[2:(nage-1),t-1]) * har[t]
		wa[t] <- har[t-1] * lmbda[t]
    totala[t] <- za[t] / wa[t]
		}

for (t in 2:nyear) {

		X[2,t] ~ dbin(Sy[t-1],totaly[t])

		X[3,t] ~ dbin(Sa[t-1],totala[t])

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

} #t



    ########################################
		# 3.2.1 Proportion of breeding female
		########################################
for (t in 1:nyear) {

 			breeds[1,t] ~ dbin(pY[t],females[1,t])
      breeds[2,t] ~ dbin(pA[t],females[1,t])

} #t


	###################################################
	# 3.3 Likelihood for population survey data
	###################################################

		#############################################################
		# 3.3.1 System process:  2 age class matrix population model
		#############################################################

	for (t in 2:t.census){

		psi1[t] <- (pY[t-1]* exp(muY[t-1])* 0.5 * Ny[t-1])+ (pA[t-1]* exp(muA[t-1])* 0.5 * Na[t-1])
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

    for (t in 1:t.census){

     Ntot[t] <- Ny[t]+Na[t]
    census[t] ~ dpois(Ntot[t])             ## female effective estimation per year given the surface of the area
                         ## Total population size observe in population survey
       } #t

}  # End Model
",fill=TRUE)
sink()


# IMPORT DATA

# 1. Census data (2003-2010)
venD<- c(0.93,0.958,0.865,1.032,0.818,0.9676,0.781,0.846)
venVD<-c(0.012,0.016,0.016,0.015,0.016,0.017,0.010,0.015)

Nven<-round(venD*238)
t.census<-length(Nven)

lam<-rep(NA,t.census)
for (t in 1:(t.census-1)) {
  lam[t+1]<- Nven[t+1]/Nven[t]              # Eq (2) from Udevitz 2012
  }		# t

census<-c(Nven/2)                      # work only on females
lmbda<-c(lam[1:5])

# 2. Age-at-death data (2003-2007)

tab<-read.table("DATA.txt",h=T,dec=",")
tab1<-subset(tab, tab$Mode!="D"& tab$Mode!="R"& tab$Age!="0")                   # remove biaised Digging out Data and juvenile class
VEN<-as.matrix(table(tab1$Age[tab1$GIC=="V"], tab1$Year[tab1$GIC=="V"]))        # work on Vendelais population

	Xit <- VEN[,3:7]                                                              # keep only 2003-2007 for a better quality data
	har<- colSums(Xit)
	nyear<- ncol(Xit)
  nage<-nrow(Xit)

  X <- rbind(Xit[1:2,],colSums(Xit[3:nage,]))

# 3. Data on productivity (2002-2006)
tab2<-subset(tab1, tab1$Info=="1"& tab1$GIC=="V")                               # keep Vendelais female with information on repro
ven<- subset(tab2, tab2$Year!="1"& tab2$Year!="2")                              # keep only 2003-2007 for a better quality data

a<-tapply(ven$Count[ven$Year=="3"],ven$Age2.[ven$Year=="3"],sum,na.rm=T)
b<-tapply(ven$Count[ven$Year=="4"],ven$Age2.[ven$Year=="4"],sum,na.rm=T)
d<-tapply(ven$Count[ven$Year=="5"],ven$Age2.[ven$Year=="5"],sum,na.rm=T)
e<-tapply(ven$Count[ven$Year=="6"],ven$Age2.[ven$Year=="6"],sum,na.rm=T)
f<-tapply(ven$Count[ven$Year=="7"],ven$Age2.[ven$Year=="7"],sum,na.rm=T)
newborns<-cbind(a,b,d,e,f)
			## Nb of newborns estimated from embryos count and placental scars

breeds<-table(ven$Age2.[ven$Status=="R"],ven$Year[ven$Status=="R"])
			## Matrix of number of breeding female per year and age class
			## assuming each breeding female have one brood

females<-table(ven$Age2.[ven$Info=="1"],ven$Year[ven$Info=="1"])
			## Matrix of total number of female sampled per year and age class


# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification   (time ~ 30mn)
thin <- 20		# thinning
chain <- 3		# number of parallel chains
iter <- 100000	# number of iterations
burn <- 50000		# number of burn-in

# 2. Define parameters to be sampled
parameters <- c("Sj","Sy","Sa","im","pY","pA","muY","muA","Ntot","lambda","Ny","Na","R","IM","H")

# 3. Create random initial values to start the chains

initstot <- function(){list(pY=c(runif(t.census,0,1)), pA=c(runif(t.census,0.5,1)), R=c(NA,rep(1,t.census-1)), JUV=c(NA,rep(1,t.census-1)), IM=c(NA,rep(1,t.census-1)),  muY=c(rnorm(t.census,1.4,0.4)), muA=c(rnorm(t.census,1.6,0.4)), Sj=c(runif((t.census-1),0,0.6)), Sy=c(runif((t.census-1),0.2,0.8)), Sa=c(runif((t.census-1),0.2,1)))}

# 4. Combine all data into a list
vulpestot <- list(census=census,t.census=t.census,lmbda=lmbda, nage=nage,nyear=nyear,newborns=newborns,breeds=breeds,females=females,Xit=Xit,X=X,har=har)

# 5. Run the model from R
venout <- bugs(vulpestot, initstot, parameters, "vulpestot.bug", n.chains=chain, n.iter=iter, n.burn=burn, n.thin=thin, debug=T, bugs.directory = bugs.dir)

print(venout,6)

# 6. Save the output
save(venout, file="venout.Rdata")



## FIGURES

#  Population size
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
year <- 2003:2010
for (i in 1:t.census){
   lower[i] <- quantile(venout$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(venout$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(venout$mean$Ntot, census, lower), na.rm = T)
m2 <- max(c(venout$mean$Ntot, census, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, t.census), ylab = "Population size", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:t.census, labels = year)
polygon(x = c(1:t.census, t.census:1), y = c(lower, upper[t.census:1]), col = "grey90", border = "grey90")
points(census, type = "l", col = "grey30", lwd = 2)
points(venout$mean$Ntot, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 50, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "blue"), bty = "n", cex = 1)


lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$Ny[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$Ny[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = venout$mean$Ny, x = (1:T)+0.5, xlim= c(1, 10), type = "b", pch = 1, ylim = c(0, 150), ylab = "Population size", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$Na[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$Na[,t], 0.975)}
points(y=venout$mean$Na, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 6, y = 150, legend = c("Yearlings", "Adults"), lty = c(1, 1),lwd = c(2, 2), pch = c(1, 16), bty = "n", cex = 1)

#  Survival rate
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census-1
for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$Sy[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$Sy[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = venout$mean$Sy, x = (1:T)+0.5, xlim= c(1, 8), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$Sa[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$Sa[,t], 0.975)}
points(y=venout$mean$Sa, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$Sj[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$Sj[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = venout$mean$Sj, x = (1:T)+0.5, xlim= c(1, 8), type = "b", pch = 17, ylim = c(0, 1), ylab = "Annual juvenile survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

#  Productivity rate : warning !! transform mu to have fec with fec=exp(mu)

FecY<-exp(venout$mean$muY)
FecA<-exp(venout$mean$muA)

par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- exp(quantile(venout$sims.list$muY[,t], 0.025))
   upper[t] <- exp(quantile(venout$sims.list$muY[,t], 0.975))}
par(mgp=c(3.8,1,0))
plot(y = exp(venout$mean$muY), x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 1, ylim = c(0, 7), ylab = "Annual Productivity (offspring per female)", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- exp(quantile(venout$sims.list$muA[,t], 0.025))
   upper[t] <- exp(quantile(venout$sims.list$muA[,t], 0.975))}
points(y=exp(venout$mean$muA), x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$pY[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$pY[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = venout$mean$pY, x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual Breeding Probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$pA[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$pA[,t], 0.975)}
points(y=venout$mean$pA, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 0.2, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")


#  immigration rate :

 par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)

lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$im[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$im[,t], 0.975)}
plot(y = venout$mean$im, x = (1:T)+0.5, xlim = c(1, 9), type = "b", pch = 16, ylim = c(0, 1.1), ylab = "Immigration rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

#  growth rate :


lower <- upper <- numeric()
T <- t.census -1
for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$lambda[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$lambda[,t], 0.975)}
plot(y = venout$mean$lambda, x = (1:T)+0.5, xlim = c(1,9), type = "b", pch = 16, ylim = c(0.5, 1.5), ylab = "Growth rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
points(y = lam[2:8], x = (1:T)+0.5, type = "b", pch = 1, cex = 1.5, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

legend(x = 1, y = 0.7, legend = c("Bugs", "Distance"), pch = c( 16, 1), bty = "n")



# DESCRIPTIVE STATISTICS

# Survival  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(venout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(venout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(venout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(venout$sims.list$Sj[,i])
   lower.h[i,1] <- quantile(venout$sims.list$Sj[,i], 0.025)
   upper.h[i,1] <- quantile(venout$sims.list$Sj[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(venout$sims.list$Sy[,i])
   lower.h[i,2] <- quantile(venout$sims.list$Sy[,i], 0.025)
   upper.h[i,2] <- quantile(venout$sims.list$Sy[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(venout$sims.list$Sa[,i])
   lower.h[i,3] <- quantile(venout$sims.list$Sa[,i], 0.025)
   upper.h[i,3] <- quantile(venout$sims.list$Sa[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(venout$sims.list$im[,i])
   lower.h[i,4] <- quantile(venout$sims.list$im[,i], 0.025)
   upper.h[i,4] <- quantile(venout$sims.list$im[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 52500)
for (i in 1:52500){
   correl.h[i,1] <- cor(venout$sims.list$lambda[i,1:4], venout$sims.list$Sj[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(venout$sims.list$lambda[i,1:4], venout$sims.list$Sy[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(venout$sims.list$lambda[i,1:4], venout$sims.list$Sa[i,1:4], use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(venout$sims.list$lambda[i,1:7], venout$sims.list$im[i,1:7], use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/52500
sum(correl.h[!is.na(correl.h[,2]),2]>0)/52500
sum(correl.h[!is.na(correl.h[,3]),3]>0)/52500
sum(correl.h[!is.na(correl.h[,4]),4]>0)/52500


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.38 (-0.36,0.32, 0.77)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.79", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Yearling survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.007 (-0.60,0.04, 0.63)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.54", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = -0.03 (-0.58,0.03, 0.64)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.52", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Immigration rate", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.25 (-0.53, 0.17,0.75)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.63", pos = 4, font = 3, cex = 0.8)

# Productivity  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(venout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(venout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(venout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(venout$sims.list$pY[,i])
   lower.h[i,1] <- quantile(venout$sims.list$pY[,i], 0.025)
   upper.h[i,1] <- quantile(venout$sims.list$pY[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(venout$sims.list$pA[,i])
   lower.h[i,2] <- quantile(venout$sims.list$pA[,i], 0.025)
   upper.h[i,2] <- quantile(venout$sims.list$pA[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(exp(venout$sims.list$muY[,i]))
   lower.h[i,3] <- quantile(exp(venout$sims.list$muY[,i]), 0.025)
   upper.h[i,3] <- quantile(exp(venout$sims.list$muY[,i]), 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(exp(venout$sims.list$muA[,i]))
   lower.h[i,4] <- quantile(exp(venout$sims.list$muA[,i]), 0.025)
   upper.h[i,4] <- quantile(exp(venout$sims.list$muA[,i]), 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 52500)
for (i in 1:52500){
   correl.h[i,1] <- cor(venout$sims.list$lambda[i,], venout$sims.list$pY[i,1:7], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(venout$sims.list$lambda[i,], venout$sims.list$pA[i,1:7], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(venout$sims.list$lambda[i,], exp(venout$sims.list$muY[i,1:7]), use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(venout$sims.list$lambda[i,], exp(venout$sims.list$muA[i,1:7]), use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/52500
sum(correl.h[!is.na(correl.h[,2]),2]>0)/52500
sum(correl.h[!is.na(correl.h[,3]),3]>0)/52500
sum(correl.h[!is.na(correl.h[,4]),4]>0)/52500


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Yearling breeding probability", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.04 (-0.54,0.04, 0.59)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.55", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult breeding probability", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.03 (-0.58,0.04, 0.65)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.54", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Yearling productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = -0.24 (-0.68,-0.10, 0.55)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.40", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.59 (-0.23,0.44, 0.84)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.87", pos = 4, font = 3, cex = 0.8)

#  removal rate :
 par(mfrow = c(1,2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)

lower <- upper <- numeric()
T <- nyear 
for (t in 1:T){
   lower[t] <- quantile(venout$sims.list$H[,t], 0.025)
   upper[t] <- quantile(venout$sims.list$H[,t], 0.975)}
plot(y = venout$mean$H, x = (1:T)+0.5, xlim = c(1,6), type = "b", pch = 16, ylim = c(0, 1), ylab = "Removal rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T), labels = 2003:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

# Removal intensity

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyear-1, ncol = 4)

for (i in 1:(nyear-1)){
   lambda.h[i] <- mean(venout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(venout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(venout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,1] <- mean(venout$sims.list$H[,i])
   lower.h[i,1] <- quantile(venout$sims.list$H[,i], 0.025)
   upper.h[i,1] <- quantile(venout$sims.list$H[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 52500)
for (i in 1:52500){
   correl.h[i,1] <- cor(venout$sims.list$lambda[i,1:5], venout$sims.list$H[i,1:5], use = "pairwise.complete.obs")

   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)



# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/52500


# Plot Fig. 11-8

linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.4), ylab = "Population growth rate", xlab = "Removal", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")

text(x = 0.1, y = 0.7, "r = -0.51 (-0.79,-0.27,0.5)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 0.65, "P(r>0) = 0.25", pos = 4, font = 3, cex = 0.8)







##########################################
## FOUGERES POPULATION
##########################################


# WRITE THE WinBUGS MODEL

sink("vulpestot.bug")
cat("
##########################################################################
#  Integrated population model for the red fox population in Vendelais France
#
#  Age structured, female-based model (2 age classes-1-year )
#  Pre-breeding census
#
#  Combination of:
#	- Distance sampling census (2003-2010): --> state-space model
#	- productivity for harvest data(2003-2006) --> Poisson regression
#	- age-at-harvest data from trapping and hunting(2003-2006) --> yearling and adult survival:
#		model by Udevitz & Gonan (2012)
#
##########################################################################
model
{
############################################################
# 1. Priors for the parameters
############################################################

    # Survival probabilities from the age-at-harvest data
for (t in 1:(t.census-1)) {
    Sy[t] ~ dunif(0,0.8)
	}

for (t in 1:(t.census-1)) {
    Sa[t] ~ dunif(0.2,1)
	}

    # Juvenile Survival probabilities
for (t in 1:(t.census-1)) {
    Sj[t] ~ dnorm(0.377,0.2)I(0,1)    # Sj prior from Devenish Nelson average of 8 population
	}

    # Productivity
for (t in 1:t.census) {
    muY[t] ~ dnorm(1.4,0.2)I(0,2)   # constrained, because larger fecundity than exp(2) = 7 is not possible
    muA[t] ~ dnorm(1.6,0.2)I(0,2)
  }

    # Proportion of breeding female
for (t in 1:t.census) {
    pY[t] ~ dunif(0,1)
    pA[t] ~ dunif(0.5,1)
  }

    # Immigration rate
for (t in 1:t.census) {
    im[t] ~ dgamma(0.1,0.5)
	}

    # Initial population sizes
    Ny[1] ~ dnorm(160,0.01)I(0,)
    Na[1] ~ dnorm(110,0.01)I(0,)

###################################################################
# 2. Derived and fixed parameters
###################################################################
    # Annual population growth rate
for (t in 2:t.census){
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
    totaly[t] <- zy[t] / wy[t]

		za[t] <- sum(Xit[2:(nage-1),t-1]) * har[t]
		wa[t] <- har[t-1] * lmbda[t]
    totala[t] <- za[t] / wa[t]
		}

for (t in 2:nyear) {

		X[2,t] ~ dbin(Sy[t-1],totaly[t])

		X[3,t] ~ dbin(Sa[t-1],totala[t])

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

} #t



    ########################################
		# 3.2.1 Proportion of breeding female
		########################################
for (t in 1:nyear) {

 			breeds[1,t] ~ dbin(pY[t],females[1,t])
      breeds[2,t] ~ dbin(pA[t],females[1,t])

} #t


	###################################################
	# 3.3 Likelihood for population survey data
	###################################################

		#############################################################
		# 3.3.1 System process:  2 age class matrix population model
		#############################################################

	for (t in 2:t.census){

		psi1[t] <- (pY[t-1]* exp(muY[t-1])* 0.5 * Ny[t-1])+ (pA[t-1]* exp(muA[t-1])* 0.5 * Na[t-1])
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

    for (t in 1:t.census){

     Ntot[t] <- Ny[t]+Na[t]
    census[t] ~ dpois(Ntot[t])             ## female effective estimation per year given the surface of the area
                         ## Total population size observe in population survey
       } #t

}  # End Model
",fill=TRUE)
sink()


# IMPORT DATA

# 1. Census data (2003-2010)
fouD<- c(1.864,1.863,2.181,2.26,2.554,NA,NA,3.107)
fouVD<-c(0.039,0.0425,0.0467,0.08,0.074,NA,NA,0.11)

Nfou<-round(fouD*248 )
t.census<-length(Nfou)

lam<-rep(NA,t.census)
for (t in 1:t.census) {
  lam[t+1]<- Nfou[t+1]/Nfou[t]              # Eq (2) from Udevitz 2012
  }		# t

census<-c(Nfou/2)                      # work only on females
lmbda<-c(lam[1:5])

# 2. Age-at-death data (2003-2007)

tab<-read.table("DATA.txt",h=T,dec=",")
tab1<-subset(tab, tab$Mode!="D"& tab$Mode!="R"& tab$Age!="0")                   # remove biaised Digging out Data and juvenile class
FOU<-as.matrix(table(tab1$Age[tab1$GIC=="F"], tab1$Year[tab1$GIC=="F"]))        # work on Vendelais population

	Xit <- FOU[,3:7]                                                              # keep only 2003-2007 for a better quality data
	har<- colSums(Xit)
	nyear<- ncol(Xit)
  nage<-nrow(Xit)

  X <- rbind(Xit[1:2,],colSums(Xit[3:nage,]))

# 3. Data on productivity (2002-2006)
tab2<-subset(tab1, tab1$Info=="1"& tab1$GIC=="F")                               # keep Vendelais female with information on repro
fou<- subset(tab2, tab2$Year!="1"& tab2$Year!="2")                              # keep only 2003-2007 for a better quality data

a<-tapply(fou$Count[fou$Year=="3"],fou$Age2.[fou$Year=="3"],sum,na.rm=T)
b<-tapply(fou$Count[fou$Year=="4"],fou$Age2.[fou$Year=="4"],sum,na.rm=T)
d<-tapply(fou$Count[fou$Year=="5"],fou$Age2.[fou$Year=="5"],sum,na.rm=T)
e<-tapply(fou$Count[fou$Year=="6"],fou$Age2.[fou$Year=="6"],sum,na.rm=T)
f<-tapply(fou$Count[fou$Year=="7"],fou$Age2.[fou$Year=="7"],sum,na.rm=T)
newborns<-cbind(a,b,d,e,f)
			## Nb of newborns estimated from embryos count and placental scars

breeds<-table(fou$Age2.[fou$Status=="R"],fou$Year[fou$Status=="R"])
			## Matrix of number of breeding female per year and age class
			## assuming each breeding female have one brood

females<-table(fou$Age2.[fou$Info=="1"],fou$Year[fou$Info=="1"])
			## Matrix of total number of female sampled per year and age class


# DEFINE SPECIFICATIONS FOR RUNNING THE MODEL

# 1. MCMC specification   (time ~ 30mn)
thin <- 20		# thinning
chain <- 3		# number of parallel chains
iter <- 100000	# number of iterations
burn <- 50000		# number of burn-in

# 2. Define parameters to be sampled
parameters <- c("Sj","Sy","Sa","im","pY","pA","muY","muA","Ntot","lambda","Ny","Na","R","IM","H")

# 3. Create random initial values to start the chains

initstot <- function(){list(pY=c(runif(t.census,0,1)), pA=c(runif(t.census,0.5,1)), R=c(NA,rep(1,t.census-1)), JUV=c(NA,rep(1,t.census-1)), IM=c(NA,rep(1,t.census-1)),  muY=c(rnorm(t.census,1.4,0.4)), muA=c(rnorm(t.census,1.6,0.4)), Sj=c(runif((t.census-1),0,0.6)), Sy=c(runif((t.census-1),0.2,0.8)), Sa=c(runif((t.census-1),0.2,1)))}

# 4. Combine all data into a list
vulpestot <- list(census=census,t.census=t.census,lmbda=lmbda, nage=nage,nyear=nyear,newborns=newborns,breeds=breeds,females=females,Xit=Xit,X=X,har=har)

# 5. Run the model from R
fouout <- bugs(vulpestot, initstot, parameters, "vulpestot.bug", n.chains=chain, n.iter=iter, n.burn=burn, n.thin=thin, debug=T, bugs.directory = bugs.dir)

print(fouout,6)

# 6. Save the output
save(fouout, file="fouout.Rdata")



## FIGURES

#  Population size
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
year <- 2003:2010
for (i in 1:t.census){
   lower[i] <- quantile(fouout$sims.list$Ntot[,i], 0.025)
   upper[i] <- quantile(fouout$sims.list$Ntot[,i], 0.975)}
m1 <- min(c(fouout$mean$Ntot, census, lower), na.rm = T)
m2 <- max(c(fouout$mean$Ntot, census, upper), na.rm = T)
plot(0, 0, ylim = c(0, m2), xlim = c(1, t.census), ylab = "Population size", xlab = " ", col = "black", type = "l",  axes = F, frame = F)
axis(2)
axis(1, at = 1:t.census, labels = year)
polygon(x = c(1:t.census, t.census:1), y = c(lower, upper[t.census:1]), col = "grey90", border = "grey90")
points(census, type = "l", col = "grey30", lwd = 2)
points(fouout$mean$Ntot, type = "l", col = "blue", lwd = 2)
legend(x = 1, y = 100, legend = c("Counts", "Estimates"), lty = c(1, 1),lwd = c(2, 2), col = c("grey30", "blue"), bty = "n", cex = 1)


lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$Ny[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$Ny[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = fouout$mean$Ny, x = (1:T)+0.5, xlim= c(1, 10), type = "b", pch = 1, ylim = c(0, 800), ylab = "Population size", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$Na[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$Na[,t], 0.975)}
points(y=fouout$mean$Na, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 800, legend = c("Yearlings", "Adults"), lty = c(1, 1),lwd = c(2, 2), pch = c(1, 16), bty = "n", cex = 1)

#  Survival rate
par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census-1
for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$Sy[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$Sy[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = fouout$mean$Sy, x = (1:T)+0.5, xlim= c(1, 8), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$Sa[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$Sa[,t], 0.975)}
points(y=fouout$mean$Sa, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$Sj[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$Sj[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = fouout$mean$Sj, x = (1:T)+0.5, xlim= c(1, 8), type = "b", pch = 17, ylim = c(0, 1), ylab = "Annual juveile survival probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

#  Productivity rate : warning !! transform mu to have fec with fec=exp(mu)

FecY<-exp(fouout$mean$muY)
FecA<-exp(fouout$mean$muA)

par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)
lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- exp(quantile(fouout$sims.list$muY[,t], 0.025))
   upper[t] <- exp(quantile(fouout$sims.list$muY[,t], 0.975))}
par(mgp=c(3.8,1,0))
plot(y = exp(fouout$mean$muY), x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 1, ylim = c(0, 7), ylab = "Annual Productivity (offspring per female)", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- exp(quantile(fouout$sims.list$muA[,t], 0.025))
   upper[t] <- exp(quantile(fouout$sims.list$muA[,t], 0.975))}
points(y=exp(fouout$mean$muA), x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 1, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")

for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$pY[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$pY[,t], 0.975)}
par(mgp=c(3.8,1,0))
plot(y = fouout$mean$pY, x = (1:T)+0.5, xlim= c(1, 9), type = "b", pch = 1, ylim = c(0, 1), ylab = "Annual Breeding Probability", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$pA[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$pA[,t], 0.975)}
points(y=fouout$mean$pA, x = (1:T)+0.5, type = "b", pch = 16, cex = 1.5, lwd = 2)
segments((1:T)+0.5, lower, (1:T)+0.5, upper, lwd = 2)

legend(x = 1, y = 0.2, legend = c("Yearlings", "Adults"), pch = c( 1, 16), bty = "n")


#  immigration rate :

 par(mfrow = c(1, 2), cex.axis = 1.2, cex.lab = 1.2, mar = c(5, 6, 1.5, 2), las = 1)

lower <- upper <- numeric()
T <- t.census
for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$im[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$im[,t], 0.975)}
plot(y = fouout$mean$im, x = (1:T)+0.5, xlim = c(1, 9), type = "b", pch = 16, ylim = c(0, 1.1), ylab = "Immigration rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2011)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

#  growth rate :


lower <- upper <- numeric()
T <- t.census -1
for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$lambda[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$lambda[,t], 0.975)}
plot(y = fouout$mean$lambda, x = (1:T)+0.5, xlim = c(1,9), type = "b", pch = 16, ylim = c(0.5, 1.5), ylab = "Growth rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
points(y = lam[2:8], x = (1:T)+0.5, type = "b", pch = 1, cex = 1.5, lwd = 2)
axis(2)
axis(1, at = 1:(T+1), labels = 2003:2010)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

legend(x = 1, y = 0.7, legend = c("Bugs", "Distance"), pch = c( 16, 1), bty = "n")



# DESCRIPTIVE STATISTICS

# Survival  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(fouout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(fouout$sims.list$Sj[,i])
   lower.h[i,1] <- quantile(fouout$sims.list$Sj[,i], 0.025)
   upper.h[i,1] <- quantile(fouout$sims.list$Sj[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(fouout$sims.list$Sy[,i])
   lower.h[i,2] <- quantile(fouout$sims.list$Sy[,i], 0.025)
   upper.h[i,2] <- quantile(fouout$sims.list$Sy[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(fouout$sims.list$Sa[,i])
   lower.h[i,3] <- quantile(fouout$sims.list$Sa[,i], 0.025)
   upper.h[i,3] <- quantile(fouout$sims.list$Sa[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(fouout$sims.list$im[,i])
   lower.h[i,4] <- quantile(fouout$sims.list$im[,i], 0.025)
   upper.h[i,4] <- quantile(fouout$sims.list$im[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 52500)
for (i in 1:52500){
   correl.h[i,1] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$Sj[i,], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$Sy[i,], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$Sa[i,], use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$im[i,1:7], use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/52500
sum(correl.h[!is.na(correl.h[,2]),2]>0)/52500
sum(correl.h[!is.na(correl.h[,3]),3]>0)/52500
sum(correl.h[!is.na(correl.h[,4]),4]>0)/52500


# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Juvenile survival", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.75 (-0.31,0.53, 0.89)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.85", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Yearling survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.57 (-0.65,0.24, 0.87)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.65", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Adult survival", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.49 (-0.52,0.26, 0.83)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.69", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Immigration rate", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.85 (-0.48, 0.26,0.91)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.69", pos = 4, font = 3, cex = 0.8)

# Productivity  on total census

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = t.census-1, ncol = 4)

for (i in 1:(t.census-1)){
   lambda.h[i] <- mean(fouout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,1] <- mean(fouout$sims.list$pY[,i])
   lower.h[i,1] <- quantile(fouout$sims.list$pY[,i], 0.025)
   upper.h[i,1] <- quantile(fouout$sims.list$pY[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,2] <- mean(fouout$sims.list$pA[,i])
   lower.h[i,2] <- quantile(fouout$sims.list$pA[,i], 0.025)
   upper.h[i,2] <- quantile(fouout$sims.list$pA[,i], 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,3] <- mean(exp(fouout$sims.list$muY[,i]))
   lower.h[i,3] <- quantile(exp(fouout$sims.list$muY[,i]), 0.025)
   upper.h[i,3] <- quantile(exp(fouout$sims.list$muY[,i]), 0.975)
   }

for (i in 1:(t.census-1)){
   Fitted.h[i,4] <- mean(exp(fouout$sims.list$muA[,i]))
   lower.h[i,4] <- quantile(exp(fouout$sims.list$muA[,i]), 0.025)
   upper.h[i,4] <- quantile(exp(fouout$sims.list$muA[,i]), 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 52500)
for (i in 1:52500){
   correl.h[i,1] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$pY[i,1:7], use = "pairwise.complete.obs")
   correl.h[i,2] <- cor(fouout$sims.list$lambda[i,], fouout$sims.list$pA[i,1:7], use = "pairwise.complete.obs")
   correl.h[i,3] <- cor(fouout$sims.list$lambda[i,], exp(fouout$sims.list$muY[i,1:7]), use = "pairwise.complete.obs")
   correl.h[i,4] <- cor(fouout$sims.list$lambda[i,], exp(fouout$sims.list$muA[i,1:7]), use = "pairwise.complete.obs")
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
sum(correl.h[!is.na(correl.h[,1]),1]>0)/52500
sum(correl.h[!is.na(correl.h[,2]),2]>0)/52500
sum(correl.h[!is.na(correl.h[,3]),3]>0)/52500
sum(correl.h[!is.na(correl.h[,4]),4]>0)/52500



# Plot Fig. 11-8
par(mfrow = c(2, 2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)
linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab = "Yearling breeding probability", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.34 (-0.58,0.18, 0.72)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.65", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5.5, 4, 1.5, 1))
plot(y = lambda.h, Fitted.h[,2], type = "n", xlim = c(0, 1), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult breeding probability", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,2], lam.lower.h, Fitted.h[,2], lam.upper.h, col = linecol)
segments(lower.h[,2], lambda.h, upper.h[,2], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,2], pch = 19, col = "black")
text(x = 0.1, y = 1.6, "r = 0.34 (-0.29,0.26, 0.71)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 1.5, "P(r>0) = 0.77", pos = 4, font = 3, cex = 0.8)

par(mar = c(5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,3], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6), ylab = "Population growth rate", xlab =  "Yearling productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,3], lam.lower.h, Fitted.h[,3], lam.upper.h, col = linecol)
segments(lower.h[,3], lambda.h, upper.h[,3], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,3], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.22 (-0.66,0.07, 0.66)", pos = 4,  font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.56", pos = 4, font = 3, cex = 0.8)

par(mar = c( 5, 4, 2, 1))
plot(y = lambda.h, Fitted.h[,4], type = "n", xlim = c(0, 7), ylim = c(0.6, 1.6),  ylab = "", xlab = "Adult productivity", frame.plot = FALSE, pch = 19)
segments(Fitted.h[,4], lam.lower.h, Fitted.h[,4], lam.upper.h, col = linecol)
segments(lower.h[,4], lambda.h, upper.h[,4], lambda.h, col = linecol)
points(y=lambda.h, Fitted.h[,4], pch = 19, col = "black")
text(x = 0.6, y = 1.6, "r = 0.33 (-0.46,0.27, 0.75)", pos = 4, font = 3, cex = 0.8)
text(x = 0.6, y = 1.5, "P(r>0) = 0.74", pos = 4, font = 3, cex = 0.8)

#  removal rate :
 par(mfrow = c(1,2), mar = c(5.5, 4, 1.5, 1), mgp=c(3, 1, 0), las = 1, cex = 1)

lower <- upper <- numeric()
T <- nyear 
for (t in 1:T){
   lower[t] <- quantile(fouout$sims.list$H[,t], 0.025)
   upper[t] <- quantile(fouout$sims.list$H[,t], 0.975)}
plot(y = fouout$mean$H, x = (1:T)+0.5, xlim = c(1,6), type = "b", pch = 16, ylim = c(0, 1), ylab = "Removal rate", xlab = "", axes = F, cex = 1.5, frame = F, lwd = 2)
axis(2)
axis(1, at = 1:(T), labels = 2003:2007)
segments((1:T)+0.5, lower, (1:T)+0.5, upper)

# Removal intensity

lambda.h <- lam.lower.h <- lam.upper.h <- numeric()
Fitted.h <- lower.h <- upper.h <- matrix(NA, nrow = nyear-1, ncol = 4)

for (i in 1:(nyear-1)){
   lambda.h[i] <- mean(fouout$sims.list$lambda[,i])
   lam.lower.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.025)
   lam.upper.h[i] <- quantile(fouout$sims.list$lambda[,i], 0.975)
   }

for (i in 1:(nyear-1)){
   Fitted.h[i,1] <- mean(fouout$sims.list$H[,i])
   lower.h[i,1] <- quantile(fouout$sims.list$H[,i], 0.025)
   upper.h[i,1] <- quantile(fouout$sims.list$H[,i], 0.975)
   }


# Calculate some correlation coefficients
correl.h <- matrix(NA, ncol = 4, nrow = 600)
for (i in 1:600){
   correl.h[i,1] <- cor(fouout$sims.list$lambda[i,1:5], fouout$sims.list$H[i,1:5], use = "pairwise.complete.obs")

   }

# Credible intervals of correlation coefficients
quantile(correl.h[,1], c(0.05, 0.5, 0.95), na.rm = TRUE)



# Compute the posterior modes of correlation coefficients
m <- density(correl.h[,1], na.rm = TRUE)
m$x[which(m$y==max(m$y))]


# Probability that correlation coefficients (r) > 0
sum(correl.h[!is.na(correl.h[,1]),1]>0)/600


# Plot Fig. 11-8

linecol <- c("grey70")
plot(y = lambda.h, Fitted.h[,1], type = "n", xlim = c(0, 0.5), ylim = c(0.6, 1.4), ylab = "Population growth rate", xlab = "Removal", frame = FALSE, pch = 19)
segments(Fitted.h[,1], lam.lower.h, Fitted.h[,1], lam.upper.h, col = linecol)
segments(lower.h[,1], lambda.h, upper.h[,1], lambda.h, col = linecol)
points(y = lambda.h, Fitted.h[,1], pch = 19, col = "black")

text(x = 0.1, y = 0.7, "r = 0.75 (-0.77,0.29,0.87)", pos = 4, font = 3, cex = 0.8)
text(x = 0.1, y = 0.65, "P(r>0) = 0.61", pos = 4, font = 3, cex = 0.8)



####################################
## CONSTRUCT TABLE FOR ANALYSIS
####################################

d<-cbind(c(2002,2003,2004,2005,2006),
        c(outot$mean$Sj[1:4],NA),c(outot$mean$Sy[1:4],NA),c(outot$mean$Sa[1:4],NA),c(outot$mean$im[1:5]),
        c(outot$mean$pY[1:5]),c(outot$mean$pA[1:5]),c(exp(outot$mean$muY[1:5])),c(exp(outot$mean$muA[1:5])),
        c(outot$mean$Ntot[1:5]),c(outot$mean$lambda[1:5]),c(outot$mean$Ny[1:5]),c(outot$mean$Na[1:5]),
        c(NA,outot$mean$R[2:5]),c(NA,outot$mean$IM[2:5]),c(outot$mean$H[1:5]))
        
v<-cbind(c(2003,2004,2005,2006,2007),
        c(venout$mean$Sj[1:4],NA),c(venout$mean$Sy[1:4],NA),c(venout$mean$Sa[1:4],NA),c(venout$mean$im[1:5]),
        c(venout$mean$pY[1:5]),c(venout$mean$pA[1:5]),c(exp(venout$mean$muY[1:5])),c(exp(venout$mean$muA[1:5])),
        c(venout$mean$Ntot[1:5]),c(venout$mean$lambda[1:5]),c(venout$mean$Ny[1:5]),c(venout$mean$Na[1:5]),
        c(NA,venout$mean$R[2:5]),c(NA,venout$mean$IM[2:5]),c(venout$mean$H[1:5]))
        
f<-cbind(c(2003,2004,2005,2006,2007),
        c(fouout$mean$Sj[1:4],NA),c(fouout$mean$Sy[1:4],NA),c(fouout$mean$Sa[1:4],NA),c(fouout$mean$im[1:5]),
        c(fouout$mean$pY[1:5]),c(fouout$mean$pA[1:5]),c(exp(fouout$mean$muY[1:5])),c(exp(fouout$mean$muA[1:5])),
        c(fouout$mean$Ntot[1:5]),c(fouout$mean$lambda[1:5]),c(fouout$mean$Ny[1:5]),c(fouout$mean$Na[1:5]),
        c(NA,fouout$mean$R[2:5]),c(NA,fouout$mean$IM[2:5]),c(fouout$mean$H[1:5]))

total<-rbind(d,v,f)
total<-as.data.frame(total)
names(total)<-c("Year","Sj","Sy","Sa","im","pY","pA","fecY","fecA","Ntot","lambda","Ny","Na","R","IM","H")

total$GIC<-c(rep("D",5),rep("V",5),rep("F",5))

attach(total)

par(mfrow=c(2,2))
plot(lambda~H,pch=16)
summary(lm(lambda~H))
abline (lm(lambda~H))
text(x = 0.4, y = 0.9, "r? = 0.001 ", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 0.85, "p(lm) = 0.90", pos = 4, font = 3, cex = 0.8)

plot(Sa~H,pch=16)
summary(lm(Sa~H))
abline (lm(Sa~H))
text(x = 0.4, y = 0.9, "r? = -0.06 ", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 0.85, "p(lm) = 0.56", pos = 4, font = 3, cex = 0.8)

plot(Sj~H,pch=16)
summary(lm(Sj~H))
abline (lm(Sj~H))
text(x = 0.4, y = 0.9, "r? = -0.078 ", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 0.85, "p(lm) = 0.66", pos = 4, font = 3, cex = 0.8)

plot(Sj~Sa,pch=16)
summary(lm(Sj~Sa))
abline (lm(Sj~Sa))
text(x = 0.4, y = 0.9, "r? = -0.06 ", pos = 4, font = 3, cex = 0.8)
text(x = 0.4, y = 0.85, "p(lm) = 0.56", pos = 4, font = 3, cex = 0.8)

