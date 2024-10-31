########  DATA CONSTRUCTION OF LIFE TABLE ANALYSIS

### import data

setwd("C:/Users/etu-devillard/Documents/1. THESE/data/renard/Survival")

data<-read.table("FOXDATA.txt",h=T,dec=",")

names (data)    

### check data and revealing sampling biais : H1 of survival estimate = random sampling of population structure                                                             
barplot(table(data$Age[data$Site=="FDC35"],data$Mode[data$Site=="FDC35"]),beside=T,xlab="Sampling mode",ylab="Number of fox sampled")
legend(x =33, y = 2000, legend = c("Juvenile", "Yearlings", "Adults 2", "3","4","5","6+"),lty = c(1, 1), lwd = c(2, 2), col=c("grey10","grey20","grey30","grey40","grey50","grey60","grey70"), cex = 1)
title("Age structure of the FDC35 data set")
       #Unearthind is highly biased, collision have scarce data. 
       #Juvenile sampling is highly biased depending on the method used.
       
fox<-subset(data,data$Mode!="D"&data$Mode!="R"&data$Age!="0")                   # Work only on hunting and trapping to avoid biais and we removed biais juvenile count
fox35<-subset(fox,fox$Site=="FDC35")                                            # Present work on FDC 35

fox35$Site<-factor(fox35$Site)
fox35$Mode<-factor(fox35$Mode)
fox35$Age<-factor(fox35$Age)

####################################################
### MODEL SELECTION ON OVERALL AGE-AT-HARVEST DATA
####################################################

# Age structure of the data

table(fox35$Age,fox35$Year)
 barplot(table(fox35$Age),xlab="Age",ylab="Frequency",main="Age structure of age-at-harvest data from FDC35")

FOX35<-subset(fox35,fox35$Year!="2000")                                         # Remove 2000 year with few data
FOX35$Year<-factor(FOX35$Year)
FOX35$Mode<-factor(FOX35$Mode)

# Sex structure of the data
barplot(table(FOX35$Sex,FOX35$Mode),beside=T,col=c("grey","black","grey50"),xlab="Sampling mode",ylab="Number of fox sampled")     
legend(x =12, y = 800, legend = c("Female", "Male", "NI"),lty = c(1, 1), lwd = c(2, 2), col=c("grey","black","grey50"), cex = 1)
title("Sex structure of the FDC35 data set")
   #  No strong differences

# Differences of age structure between sampling method
m<-table(FOX35$Age,FOX35$Mode)
barplot(table(FOX35$Age,FOX35$Mode),beside=T,xlab="Sampling mode",ylab="Number of fox sampled")
legend(x = 1, y = 1000, legend = c("Yearlings", "Adults 2", "3","4","5","6+"),lty = c(1, 1), lwd = c(2, 2), col=c("grey10","grey20","grey30","grey40","grey50","grey60","grey70"), cex = 1)
title("Age structure of the FDC35 data set")
chisq.test(m, simulate.p.value = T, B=10000)

barplot(t(m)/colSums(m),beside=T,xlab="Sampling mode",ylab="Prop of fox sampled",legend.text=c("Hunting","Trapping"))                                            
m1<- (t(m)/colSums(m))

# Differences of age structure between site
m<-table(FOX35$Age,FOX35$GIC)
barplot(t(m)/colSums(m),beside=T,xlab="GIC",ylab="Prop of fox sampled",legend.text=c("DOM","FOU","VEN"))                                            
m1<-rbind(m[1:4,],colSums(m[5:10,]))
chisq.test(m[,c(2,3,5)], simulate.p.value = T, B=10000)
# Check the assumption of stable age structure

m<-table(FOX35$Age,FOX35$Year)
barplot(t(m)/colSums(m),beside=T,xlab="Age structure",ylab="Prop of fox sampled",legend.text=c("2001","2002","2003","2004","2005","2006","2007"))                                            
m1<-rbind(m[1:4,],colSums(m[5:10,]))
chisq.test(m1)          # unstable age structure due to 2001 and 2002 where only Domagne data
chisq.test(m1[,3:7])    # stable age structure between 2003 and 2007 in all FDC35


### lambda estimates from independant distance sampling data : N entre 2002 et 2010 [ rq: CV = ecart-type / moyenne= racine(variance)/mean

DENS<-read.table("Total.txt",h=T,dec=".")
names (DENS)                                                                 

lamD<- DENS$LAMBDA[DENS$GIC=="DOM"]
VarLamD<-DENS$VARLAM[DENS$GIC=="DOM"]

lamV<- c(NA,DENS$LAMBDA[DENS$GIC=="VEN"])
VarLamV<- c(NA,DENS$VARLAM[DENS$GIC=="VEN"])

lamF<- c(NA,DENS$LAMBDA[DENS$GIC=="FOU"])
VarLamF<- c(NA,DENS$VARLAM[DENS$GIC=="FOU"])

lambda<-colMeans(rbind(lamD,lamV,lamF),na.rm=T)
VarLam<-colMeans(rbind(VarLamD,VarLamV,VarLamF),na.rm=T)

### survival estimate between 2003 and 2006 : 2003 constraint by absence of lambda in 2002 

Tot<-table(FOX35$Age,FOX35$Year)
Xit<-Tot[,3:7]
nyear<- ncol(Xit)
nage<-nrow(Xit)
lambda<-lambda[2:7]
VarLam<-VarLam[2:7]


# table of Cit

Cit<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      Cit[i,t]<- Xit[i,t]/sum(Xit[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t


### Model: no assumptions about equality of survival rates

# table of Sit

Sit<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:(nage-1)) {         		                                       # loop
      Sit[i,t]<- Cit[i+1,t+1]*lambda[t]/Cit[i,t]                               # Eq (6) from Udevitz 2012
      }   # i
}		# t

Sit<- replace(Sit,which(is.nan(Sit)==TRUE),0)                                    # replace impossible value by 0
Sit<- replace(Sit,which(is.infinite(Sit)==TRUE),0)


nbSit<-36                                                                        # nb of survival parameters for AIC 

Var_Sit<-  matrix(0, nage, nyear-1)                                              # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit[i,t]<- (Sit[i,t]^2)*((1-Cit[i,t])/Xit[i,t]+(1-Cit[i+1,t+1])/Xit[i+1,t+1])+ VarLam[t]*(Cit[i+1,t+1]/Cit[i,t])^2                               
      }   # i
}		# t

Var_Sit<- replace(Var_Sit,which(is.nan(Var_Sit)==TRUE),0)                        # replace impossible value by 0
Var_Sit<- replace(Var_Sit,which(is.infinite(Var_Sit)==TRUE),0)

SD_Sit <- sqrt(Var_Sit)                                                          # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp<- Cit[,1:(nyear-1)]                                                      # because there is no age related group

nbCit_exp<-length(Cit_exp)                                                       # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p<-  matrix(0, nage, nyear) 
lik<-  matrix(0, nage, nyear) 
lik[,1]<- Xit[,1]*log(Cit_exp[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                        # loop
      p[i,t]<-Cit_exp[i-1,t-1]*Sit[i-1,t-1]/lambda[t-1]     
      }   # i   
      p[1,t]<-1-sum(p[2:nage,t]) 
      
      for (i in 1:nage) {                                                        # loop
      lik[i,t]<-Xit[i,t]*log(p[i,t])     
      }   # i   
      
}		# t

lik<- replace(lik,which(is.nan(lik)==TRUE),0)                                    # replace impossible value by 0
lik<- replace(lik,which(is.infinite(lik)==TRUE),0)
loglik<- colSums (lik)

# Derive factorial

fact_xit <- lfactorial(Xit)
fact_xit <- sum(fact_xit)

fact_sigma_xit <- sum(lfactorial(colSums(Xit)))

# Derive loglikelihood

loglikelihood <- fact_sigma_xit-fact_xit+sum(loglik)

# Derive Aic

AIC<- 2*(nbSit+nbCit_exp)-2*loglikelihood

## model summary
Sit
SD_Sit
Cit_exp
colSums(Cit_exp)
p
colSums(p)
lik
loglik
loglikelihood
nbSit
nbCit_exp
AIC


### Model1:  Assumes survival rates are equal for ages 6-10

# table of Sit

Sit1<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:5) {         		                                               # loop
      Sit1[i,t]<- Cit[i+1,t+1]*lambda[t]/Cit[i,t]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 6:(nage-1)) {         		                                       # loop
      Sit1[i,t]<- sum(Cit[7:nage,t+1])*lambda[t]/sum(Cit[6:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

Sit1<- replace(Sit1,which(is.nan(Sit1)==TRUE),0)                                 # replace impossible value by 0
Sit1<- replace(Sit1,which(is.infinite(Sit1)==TRUE),0)

nbSit1<-24                                                                     # nb of survival parameters for AIC 

Var_Sit1<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:5) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit1[i,t]<- (Sit1[i,t]^2)*((1-Cit[i,t])/Xit[i,t]+(1-Cit[i+1,t+1])/Xit[i+1,t+1])+ VarLam[t]*(Cit[i+1,t+1]/Cit[i,t])^2                               
      }   # i
      for (i in 6:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit1[i,t]<- (Sit1[i,t]^2)*((1-sum(Cit[6:(nage-1),t]))/sum(Xit[6:(nage-1),t])+(1-sum(Cit[7:nage,t+1]))/sum(Xit[7:nage,t+1]))+ VarLam[t]*(sum(Cit[7:nage,t+1])/sum(Cit[6:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Sit1<- replace(Var_Sit1,which(is.nan(Var_Sit1)==TRUE),0)                     # replace impossible value by 0
Var_Sit1<- replace(Var_Sit1,which(is.infinite(Var_Sit1)==TRUE),0)

SD_Sit1 <- sqrt(Var_Sit1)                                                        # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp1<-  matrix(0, nage, nyear-1)                                             # for age related group
for (t in 1:(nyear-1)) { 
      for (i in 1:5) {         		                                               # Eq (15) from Udevitz 2012
      Cit_exp1[i,t]<- Cit[i,t]                               
      }                                                                          # loop
      for (i in 6:(nage-1)) {         		                                       # Eq (15) from Udevitz 2012
      Cit_exp1[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[6:(nage-1),t])))/(sum(Xit[,t])*(sum(Xit[6:(nage-1),t])+sum(Xit[7:nage,t+1])) )
      }                                                  
}		# t
Cit_exp1[nage,]<-Cit[nage,1:(nyear-1)]
Cit_exp1[1,]<-1-colSums(Cit_exp1[2:nage,])


nbCit_exp1<-length(Cit_exp1)                                                      # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p1<-  matrix(0, nage, nyear) 
lik1<-  matrix(0, nage, nyear) 
lik1[,1]<- Xit[,1]*log(Cit_exp1[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                          # loop
      p1[i,t]<-Cit_exp1[i-1,t-1]*Sit1[i-1,t-1]/lambda[t-1]     
      }   # i   
      p1[1,t]<-1-sum(p1[2:nage,t]) 
      
      for (i in 1:nage) {                                                          # loop
      lik1[i,t]<-Xit[i,t]*log(p1[i,t])     
      }   # i   
      
}		# t

lik1<- replace(lik1,which(is.nan(lik1)==TRUE),0)                                    # replace impossible value by 0
lik1<- replace(lik1,which(is.infinite(lik1)==TRUE),0)
loglik1<- colSums (lik1)

# Derive factorial     idem for all model !!

# Derive loglikelihood

loglikelihood1 <- fact_sigma_xit-fact_xit+sum(loglik1)

# Derive Aic

AIC1<- 2*(nbSit1+nbCit_exp1)-2*loglikelihood1

## model summary
Sit1
SD_Sit1
Cit_exp1
colSums(Cit_exp1)
p1
colSums(p1)
lik1
loglik1
loglikelihood1
nbSit1
nbCit_exp1
AIC1


### Model2:  Assumes survival rates are equal for ages 3-5, 6-10

# table of Sit

Sit2<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:2) {         		                                               # loop
      Sit2[i,t]<- Cit[i+1,t+1]*lambda[t]/Cit[i,t]                                # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:5) {         		                                               # loop
      Sit2[i,t]<- sum(Cit[4:6,t+1])*lambda[t]/sum(Cit[3:5,t])                    # Eq (6) from Udevitz 2012
      }   # i
       for (i in 6:(nage-1)) {         		                                       # loop
      Sit2[i,t]<- sum(Cit[7:nage,t+1])*lambda[t]/sum(Cit[6:(nage-1),t])                   # Eq (6) from Udevitz 2012
      }   # i

}		# t

Sit2<- replace(Sit2,which(is.nan(Sit2)==TRUE),0)                                 # replace impossible value by 0
Sit2<- replace(Sit2,which(is.infinite(Sit2)==TRUE),0)

nbSit2<-16                                                                        # nb of survival parameters for AIC 

Var_Sit2<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit2[i,t]<- (Sit2[i,t]^2)*((1-Cit[i,t])/Xit[i,t]+(1-Cit[i+1,t+1])/Xit[i+1,t+1])+ VarLam[t]*(Cit[i+1,t+1]/Cit[i,t])^2                               
      }   # i
      for (i in 3:5) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit2[i,t]<- (Sit2[i,t]^2)*((1-sum(Cit[3:5,t]))/sum(Xit[3:5,t])+(1-sum(Cit[4:6,t+1]))/sum(Xit[4:6,t+1]))+ VarLam[t]*(sum(Cit[4:6,t+1])/sum(Cit[3:5,t]))^2                               
      }
      for (i in 6:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit2[i,t]<- (Sit2[i,t]^2)*((1-sum(Cit[6:(nage-1),t]))/sum(Xit[6:(nage-1),t])+(1-sum(Cit[7:nage,t+1]))/sum(Xit[7:nage,t+1]))+ VarLam[t]*(sum(Cit[7:nage,t+1])/sum(Cit[6:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Sit2<- replace(Var_Sit2,which(is.nan(Var_Sit2)==TRUE),0)                     # replace impossible value by 0
Var_Sit2<- replace(Var_Sit2,which(is.infinite(Var_Sit2)==TRUE),0)

SD_Sit2 <- sqrt(Var_Sit2)                                                        # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp2<-  matrix(0, nage, nyear-1)                                             # for age related group
for (t in 1:(nyear-1)) { 
      for (i in 1:2) {         		                                                 # Eq (15) from Udevitz 2012
      Cit_exp2[i,t]<- Cit[i,t]                               
      }                                                                           # loop
      for (i in 3:5) {         		                                                # Eq (15) from Udevitz 2012
      Cit_exp2[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[3:5,t])))/(sum(Xit[,t])*(sum(Xit[3:5,t])+sum(Xit[4:6,t+1])) )
      } 
      for (i in 6:(nage-1)) {         		                                       # Eq (15) from Udevitz 2012
      Cit_exp2[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[6:(nage-1),t])))/(sum(Xit[,t])*(sum(Xit[6:(nage-1),t])+sum(Xit[7:nage,t+1])) )
      }                                                  
}		# t
Cit_exp2[nage,]<-Cit[nage,1:(nyear-1)]
Cit_exp2[1,]<-1-colSums(Cit_exp2[2:nage,])


nbCit_exp2<-length(Cit_exp2)                                                     # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p2<-  matrix(0, nage, nyear) 
lik2<-  matrix(0, nage, nyear) 
lik2[,1]<- Xit[,1]*log(Cit_exp2[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                          # loop
      p2[i,t]<-Cit_exp2[i-1,t-1]*Sit2[i-1,t-1]/lambda[t-1]     
      }   # i   
      p2[1,t]<-1-sum(p2[2:nage,t]) 
      
      for (i in 1:nage) {                                                          # loop
      lik2[i,t]<-Xit[i,t]*log(p2[i,t])     
      }   # i   
      
}		# t

lik2<- replace(lik2,which(is.nan(lik2)==TRUE),0)                                    # replace impossible value by 0
lik2<- replace(lik2,which(is.infinite(lik2)==TRUE),0)
loglik2<- colSums (lik2)

# Derive factorial     idem for all model !!

# Derive loglikelihood

loglikelihood2 <- fact_sigma_xit-fact_xit+sum(loglik2)

# Derive Aic

AIC2<- 2*(nbSit2+nbCit_exp2)-2*loglikelihood2

 ## model summary
Sit2
SD_Sit2
Cit_exp2
colSums(Cit_exp2)
p2
colSums(p2)
lik2
loglik2
loglikelihood2
nbSit2
nbCit_exp2
AIC2


### Model3:  Assumes survival rates are equal for ages 2-5, 6-10

# table of Sit

Sit3<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:1) {         		                                               # loop
      Sit3[i,t]<- Cit[i+1,t+1]*lambda[t]/Cit[i,t]                                # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:5) {         		                                               # loop
      Sit3[i,t]<- sum(Cit[3:6,t+1])*lambda[t]/sum(Cit[2:5,t])                    # Eq (6) from Udevitz 2012
      }   # i
      for (i in 6:(nage-1)) {         		                                       # loop
      Sit3[i,t]<- sum(Cit[7:nage,t+1])*lambda[t]/sum(Cit[6:(nage-1),t])                   # Eq (6) from Udevitz 2012
      }   # i

}		# t

Sit3<- replace(Sit3,which(is.nan(Sit3)==TRUE),0)                                 # replace impossible value by 0
Sit3<- replace(Sit3,which(is.infinite(Sit3)==TRUE),0)

nbSit3<-12                                                                       # nb of survival parameters for AIC 

Var_Sit3<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit3[i,t]<- (Sit3[i,t]^2)*((1-Cit[i,t])/Xit[i,t]+(1-Cit[i+1,t+1])/Xit[i+1,t+1])+ VarLam[t]*(Cit[i+1,t+1]/Cit[i,t])^2                               
      }   # i
      for (i in 2:5) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit3[i,t]<- (Sit3[i,t]^2)*((1-sum(Cit[2:5,t]))/sum(Xit[2:5,t])+(1-sum(Cit[3:6,t+1]))/sum(Xit[3:6,t+1]))+ VarLam[t]*(sum(Cit[3:6,t+1])/sum(Cit[2:5,t]))^2                               
      }
      for (i in 6:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit3[i,t]<- (Sit3[i,t]^2)*((1-sum(Cit[6:(nage-1),t]))/sum(Xit[6:(nage-1),t])+(1-sum(Cit[7:nage,t+1]))/sum(Xit[7:nage,t+1]))+ VarLam[t]*(sum(Cit[7:nage,t+1])/sum(Cit[6:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Sit3<- replace(Var_Sit3,which(is.nan(Var_Sit3)==TRUE),0)                     # replace impossible value by 0
Var_Sit3<- replace(Var_Sit3,which(is.infinite(Var_Sit3)==TRUE),0)

SD_Sit3 <- sqrt(Var_Sit3)                                                        # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp3<-  matrix(0, nage, nyear-1)                                                                   
for (t in 1:(nyear-1)) { 
      for (i in 1:1) {         		                                               # Eq (15) from Udevitz 2012
      Cit_exp3[i,t]<- Cit[i,t]                               
      }                                                                           # loop
      for (i in 2:5) {         		                                                # Eq (15) from Udevitz 2012
      Cit_exp3[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[2:5,t])))/(sum(Xit[,t])*(sum(Xit[2:5,t])+sum(Xit[3:6,t+1])) )
      } 
      for (i in 6:(nage-1)) {         		                                       # Eq (15) from Udevitz 2012
      Cit_exp3[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[6:(nage-1),t])))/(sum(Xit[,t])*(sum(Xit[6:(nage-1),t])+sum(Xit[7:nage,t+1])) )
      }                                                  
}		# t
Cit_exp3[nage,]<-Cit[nage,1:(nyear-1)]
Cit_exp3[1,]<-1-colSums(Cit_exp3[2:nage,])


nbCit_exp3<-length(Cit_exp3)                                                      # nb of survival parameters for AIC 



### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p3<-  matrix(0, nage, nyear) 
lik3<-  matrix(0, nage, nyear) 
lik3[,1]<- Xit[,1]*log(Cit_exp3[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                          # loop
      p3[i,t]<-Cit_exp3[i-1,t-1]*Sit3[i-1,t-1]/lambda[t-1]     
      }   # i   
      p3[1,t]<-1-sum(p3[2:nage,t]) 
      
      for (i in 1:nage) {                                                          # loop
      lik3[i,t]<-Xit[i,t]*log(p3[i,t])     
      }   # i   
      
}		# t

lik3<- replace(lik3,which(is.nan(lik3)==TRUE),0)                                    # replace impossible value by 0
lik3<- replace(lik3,which(is.infinite(lik3)==TRUE),0)
loglik3<- colSums (lik3)

# Derive factorial     idem for all model !!

# Derive loglikelihood

loglikelihood3 <- fact_sigma_xit-fact_xit+sum(loglik3)

# Derive Aic

AIC3<- 2*(nbSit3+nbCit_exp3)-2*loglikelihood3

 ## model summary
Sit3
SD_Sit3
Cit_exp3
colSums(Cit_exp3)
p3
colSums(p3)
lik3
loglik3
loglikelihood3
nbSit3
nbCit_exp3
AIC3


   
   
### Model4:  Assumes survival rates are equal for ages 2-10

# table of Sit

Sit4<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      Sit4[1,t]<- Cit[2,t+1]*lambda[t]/Cit[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      Sit4[i,t]<- sum(Cit[3:nage,t+1])*lambda[t]/sum(Cit[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

Sit4<- replace(Sit4,which(is.nan(Sit4)==TRUE),0)                                 # replace impossible value by 0
Sit4<- replace(Sit4,which(is.infinite(Sit4)==TRUE),0)

nbSit4<-8                                                                       # nb of survival parameters for AIC 

Var_Sit4<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_Sit4[1,t]<- (Sit4[1,t]^2)*((1-Cit[1,t])/Xit[1,t]+(1-Cit[2,t+1])/Xit[2,t+1])+ VarLam[t]*(Cit[2,t+1]/Cit[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit4[i,t]<- (Sit4[i,t]^2)*((1-sum(Cit[2:(nage-1),t]))/sum(Xit[2:(nage-1),t])+(1-sum(Cit[3:nage,t+1]))/sum(Xit[3:nage,t+1]))+ VarLam[t]*(sum(Cit[3:nage,t+1])/sum(Cit[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Sit4<- replace(Var_Sit4,which(is.nan(Var_Sit4)==TRUE),0)                     # replace impossible value by 0
Var_Sit4<- replace(Var_Sit4,which(is.infinite(Var_Sit4)==TRUE),0)

SD_Sit4 <- sqrt(Var_Sit4)                                                        # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp4<-  matrix(0, nage, nyear-1)                                                                    
for (t in 1:(nyear-1)) { 
                                                                                 # loop
      for (i in 2:(nage-1)) {         		                                       # Eq (15) from Udevitz 2012
      Cit_exp4[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[2:(nage-1),t])))/(sum(Xit[,t])*(sum(Xit[2:(nage-1),t])+sum(Xit[3:nage,t+1])) )
      }                                                  
}		# t
Cit_exp4[nage,]<-Cit[nage,1:(nyear-1)]
Cit_exp4[1,]<-1-colSums(Cit_exp4[2:nage,])


nbCit_exp4<-length(Cit_exp4)                                                      # nb of survival parameters for AIC 
                                            


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p4<-  matrix(0, nage, nyear) 
lik4<-  matrix(0, nage, nyear) 
lik4[,1]<- Xit[,1]*log(Cit_exp4[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                          # loop
      p4[i,t]<-Cit_exp4[i-1,t-1]*Sit4[i-1,t-1]/lambda[t-1]     
      }   # i   
      p4[1,t]<-1-sum(p4[2:nage,t]) 
      
      for (i in 1:nage) {                                                          # loop
      lik4[i,t]<-Xit[i,t]*log(p4[i,t])     
      }   # i   
      
}		# t

lik4<- replace(lik4,which(is.nan(lik4)==TRUE),0)                                   # replace impossible value by 0
lik4<- replace(lik4,which(is.infinite(lik4)==TRUE),0)
loglik4<- colSums (lik4)

# Derive factorial     idem for all model !!

# Derive loglikelihood

loglikelihood4 <- fact_sigma_xit-fact_xit+sum(loglik4)

# Derive Aic

AIC4<- 2*(nbSit4+nbCit_exp4)-2*loglikelihood4

 ## model summary
Sit4
SD_Sit4
Cit_exp4
colSums(Cit_exp4)
p4
colSums(p4)
lik4
loglik4
loglikelihood4
nbSit4
nbCit_exp4
AIC4
 
 
### Model5:  Assumes survival rates are equal for ages 1-8

# table of Sit

Sit5<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
       
      for (i in 1:(nage-1)) {         		                                       # loop
      Sit5[i,t]<- sum(Cit[2:nage,t+1])*lambda[t]/sum(Cit[1:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

Sit5<- replace(Sit5,which(is.nan(Sit5)==TRUE),0)                                 # replace impossible value by 0
Sit5<- replace(Sit5,which(is.infinite(Sit5)==TRUE),0)

nbSit5<-4                                                                        # nb of survival parameters for AIC 

Var_Sit5<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
     
      for (i in 1:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit5[i,t]<- (Sit5[i,t]^2)*((1-sum(Cit[1:(nage-1),t]))/sum(Xit[1:(nage-1),t])+(1-sum(Cit[2:nage,t+1]))/sum(Xit[2:nage,t+1]))+ VarLam[t]*(sum(Cit[2:nage,t+1])/sum(Cit[1:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Sit5<- replace(Var_Sit5,which(is.nan(Var_Sit5)==TRUE),0)                     # replace impossible value by 0
Var_Sit5<- replace(Var_Sit5,which(is.infinite(Var_Sit5)==TRUE),0)

SD_Sit5 <- sqrt(Var_Sit5)                                                        # matrix of standard deviation of survival estimate 

# table of C*it

Cit_exp5<-  matrix(0, nage, nyear-1)                                                                    
for (t in 1:(nyear-1)) { 
                                                                                 # loop
      for (i in 2:(nage-1)) {         		                                       # Eq (15) from Udevitz 2012
      Cit_exp5[i,t]<- ((Xit[i,t]+Xit[i+1,t+1])*(sum(Xit[1:(nage-1),t])))/(sum(Xit[,t])*(sum(Xit[1:(nage-1),t])+sum(Xit[2:nage,t+1])) )
      }                                                  
}		# t
Cit_exp5[nage,]<-Cit[nage,1:(nyear-1)]
Cit_exp5[1,]<-1-colSums(Cit_exp5[2:nage,])


nbCit_exp5<-length(Cit_exp5)                                                     # nb of survival parameters for AIC 
                                            


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

p5<-  matrix(0, nage, nyear) 
lik5<-  matrix(0, nage, nyear) 
lik5[,1]<- Xit[,1]*log(Cit_exp5[,1])
for (t in 2:nyear) { 
      for (i in 2:nage) {                                                          # loop
      p5[i,t]<-Cit_exp5[i-1,t-1]*Sit5[i-1,t-1]/lambda[t-1]     
      }   # i   
      p5[1,t]<-1-sum(p5[2:nage,t]) 
      
      for (i in 1:nage) {                                                          # loop
      lik5[i,t]<-Xit[i,t]*log(p5[i,t])     
      }   # i   
      
}		# t

lik5<- replace(lik5,which(is.nan(lik5)==TRUE),0)                                 # replace impossible value by 0
lik5<- replace(lik5,which(is.infinite(lik5)==TRUE),0)
loglik5<- colSums (lik5)

# Derive factorial     idem for all model !!

# Derive loglikelihood

loglikelihood5 <- fact_sigma_xit-fact_xit+sum(loglik5)

# Derive Aic

AIC5<- 2*(nbSit5+nbCit_exp5)-2*loglikelihood5

 ## model summary
Sit5
SD_Sit5
Cit_exp5
colSums(Cit_exp5)
p5
colSums(p5)
lik5
loglik5
loglikelihood5
nbSit5
nbCit_exp5
AIC5
 
 
 
##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: Model6: no assumptions about equality of survival rates

# VErify assumption of stable age structure


# table of Cit

Cit6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      Cit6[i]<- sum(Xit[i,])/sum(Xit)                                            # Eq (12) from Udevitz 2012
      }   # i

# table of Sit

Sit6<- matrix(0, nage)  

      for (i in 1:(nage-1)) {         		                                       # loop
      Sit6[i]<- Cit6[i+1]*mean(lambda[1:(nyear-1)])/Cit6[i]                              # Eq (6) from Udevitz 2012
      }   # i

Sit6<- replace(Sit6,which(is.nan(Sit6)==TRUE),0)                                 # replace impossible value by 0
Sit6<- replace(Sit6,which(is.infinite(Sit6)==TRUE),0)


nbSit6<-nage-1                                                                   # nb of survival parameters for AIC 

Var_Sit6<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit6[i]<- (Sit6[i]^2)*((1-Cit6[i])/Xit[i]+(1-Cit6[i+1])/Xit[i+1])+ mean(VarLam[1:(nyear-1)])*(Cit6[i+1]/Cit6[i])^2                               
      }   # i


Var_Sit6<- replace(Var_Sit6,which(is.nan(Var_Sit6)==TRUE),0)                     # replace impossible value by 0
Var_Sit6<- replace(Var_Sit6,which(is.infinite(Var_Sit6)==TRUE),0)

SD_Sit6 <- sqrt(Var_Sit6)                                                        # matrix of standard deviation of survival estimate 

# table of Lit

Lit6<- matrix(0, nage)  
Lit6[1]<-1
      for (i in 2:(nage-1)) {         		                                       # loop
      Lit6[i]<- prod(Sit6[1:i-1])                                                # Eq (12) from Udevitz 2012
      }   # i

nbLit6<-nage-1                                                                   # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

mlb<-rep(0, nage)

for (i in 1:nage) {                                                              # loop
      mlb[i]<-mean(lambda[1:(nyear-1)])^i    
      }   # i   
      
p6<-  rep(0, nage) 
lik6<-  matrix(0, nage, nyear) 

      for (i in 2:nage) {                                                        # loop
      p6[i]<-Lit6[i]/(mlb[i]*sum(Lit6/mlb))    
      }   # i   
      
      p6[1]<-1-sum(p6[2:nage]) 
      
for (t in 1:nyear) {       
      for (i in 1:nage) {                                                        # loop
      lik6[i,t]<-Xit[i,t]*log(p6[i])     
      }   # i   
  }   # t       

lik6<- replace(lik6,which(is.nan(lik6)==TRUE),0)                                 # replace impossible value by 0
lik6<- replace(lik6,which(is.infinite(lik6)==TRUE),0)
loglik6<- colSums (lik6)

# Derive loglikelihood

loglikelihood6 <- fact_sigma_xit-fact_xit+sum(loglik6)

# Derive Aic

AIC6<- 2*(nbSit6+nbLit6)-2*loglikelihood6

## model summary
Sit6
SD_Sit6
Lit6
Cit6
sum(Cit6)
p6
sum(p6)
lik6
loglik6
loglikelihood6
nbSit6
nbLit6
AIC6  

 ### Model 7 Time invariant : assumptions of equal survival between 2-8
# table of Sit

Sit7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      Sit7[i]<- Cit6[i+1]*mean(lambda[1:(nyear-1)])/Cit6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      Sit7[i]<- sum(Cit6[3:nage])*mean(lambda[1:(nyear-1)])/sum(Cit6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
Sit7<- replace(Sit7,which(is.nan(Sit7)==TRUE),0)                                 # replace impossible value by 0
Sit7<- replace(Sit7,which(is.infinite(Sit7)==TRUE),0)


nbSit7<-2                                                                        # nb of survival parameters for AIC 

Var_Sit7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit7[i]<- (Sit7[i]^2)*((1-Cit6[i])/Xit[i]+(1-Cit6[i+1])/Xit[i+1])+ mean(VarLam[1:(nyear-1)])*(Cit6[i+1]/Cit6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit7[i]<- (Sit7[i]^2)*((1-sum(Cit6[2:(nage-1)]))/sum(Xit[2:(nage-1)])+(1-sum(Cit6[3:nage]))/sum(Xit[3:nage]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[3:nage])/sum(Cit6[2:(nage-1)]))^2                               
      }   # i

Var_Sit7<- replace(Var_Sit7,which(is.nan(Var_Sit7)==TRUE),0)                     # replace impossible value by 0
Var_Sit7<- replace(Var_Sit7,which(is.infinite(Var_Sit7)==TRUE),0)

SD_Sit7 <- sqrt(Var_Sit7)                                                        # matrix of standard deviation of survival estimate 

# table of Lit

Lit7<- matrix(0, nage)  
Lit7[1]<-1
      for (i in 2:(nage-1)) {         		                                       # loop
      Lit7[i]<- prod(Sit7[1:i-1])                                                # Eq (12) from Udevitz 2012
      }   # i

nbLit7<-nage-1                                                                   # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

mlb<-rep(0, nage)

for (i in 1:nage) {                                                              # loop
      mlb[i]<-mean(lambda[1:(nyear-1)])^i    
      }   # i   
      
p7<-  rep(0, nage) 
lik7<-  matrix(0, nage, nyear) 

      for (i in 2:nage) {                                                        # loop
      p7[i]<-Lit7[i]/(mlb[i]*sum(Lit7/mlb))    
      }   # i   
      
      p7[1]<-1-sum(p7[2:nage]) 
      
for (t in 1:nyear) {       
      for (i in 1:nage) {                                                        # loop
      lik7[i,t]<-Xit[i,t]*log(p7[i])     
      }   # i   
  }   # t       

lik7<- replace(lik7,which(is.nan(lik7)==TRUE),0)                                 # replace impossible value by 0
lik7<- replace(lik7,which(is.infinite(lik7)==TRUE),0)
loglik7<- colSums (lik7)

# Derive loglikelihood

loglikelihood7 <- fact_sigma_xit-fact_xit+sum(loglik7)

# Derive Aic

AIC7<- 2*(nbSit7+nbLit7)-2*loglikelihood7

## model summary
Sit7
SD_Sit7
Lit7
Cit6
sum(Cit6)
p7
sum(p7)
lik7
loglik7
loglikelihood7
nbSit7
nbLit7
AIC7  
 
 ### Model 8 Time invariant : assumptions of equal survival between 1-8
# table of Sit

Sit8<- matrix(0, nage)  

      for (i in 1:(nage-1)) {         		                                       # loop
      Sit8[i]<- sum(Cit6[2:nage])*mean(lambda[1:(nyear-1)])/sum(Cit6[1:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
Sit8<- replace(Sit8,which(is.nan(Sit8)==TRUE),0)                                 # replace impossible value by 0
Sit8<- replace(Sit8,which(is.infinite(Sit8)==TRUE),0)


nbSit8<-1                                                                        # nb of survival parameters for AIC 

Var_Sit8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit8[i]<- (Sit8[i]^2)*((1-sum(Cit6[1:(nage-1)]))/sum(Xit[1:(nage-1)])+(1-sum(Cit6[2:nage]))/sum(Xit[2:nage]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[2:nage])/sum(Cit6[1:(nage-1)]))^2                               
      }   # i

Var_Sit8<- replace(Var_Sit8,which(is.nan(Var_Sit8)==TRUE),0)                     # replace impossible value by 0
Var_Sit8<- replace(Var_Sit8,which(is.infinite(Var_Sit8)==TRUE),0)

SD_Sit8 <- sqrt(Var_Sit8)                                                        # matrix of standard deviation of survival estimate 

# table of Lit

Lit8<- matrix(0, nage)  
Lit8[1]<-1
      for (i in 2:(nage-1)) {         		                                       # loop
      Lit8[i]<- prod(Sit8[1:i-1])                                                # Eq (12) from Udevitz 2012
      }   # i

nbLit8<-nage-1                                                                   # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

mlb<-rep(0, nage)

for (i in 1:nage) {                                                              # loop
      mlb[i]<-mean(lambda[1:(nyear-1)])^i    
      }   # i   
      
p8<-  rep(0, nage) 
lik8<-  matrix(0, nage, nyear) 

      for (i in 2:nage) {                                                        # loop
      p8[i]<-Lit8[i]/(mlb[i]*sum(Lit8/mlb))    
      }   # i   
      
      p8[1]<-1-sum(p8[2:nage]) 
      
for (t in 1:nyear) {       
      for (i in 1:nage) {                                                        # loop
      lik8[i,t]<-Xit[i,t]*log(p8[i])     
      }   # i   
  }   # t       

lik8<- replace(lik8,which(is.nan(lik8)==TRUE),0)                                 # replace impossible value by 0
lik8<- replace(lik8,which(is.infinite(lik8)==TRUE),0)
loglik8<- colSums (lik8)

# Derive loglikelihood

loglikelihood8 <- fact_sigma_xit-fact_xit+sum(loglik8)

# Derive Aic

AIC8<- 2*(nbSit8+nbLit8)-2*loglikelihood8

## model summary
Sit8
SD_Sit8
Lit8
Cit6
sum(Cit6)
p8
sum(p8)
lik8
loglik8
loglikelihood8
nbSit8
nbLit8
AIC8   

 ### Model 9 Time invariant : assumptions of equal survival between 1;2-5;6-10
# table of Sit

Sit9<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      Sit9[i]<- Cit6[i+1]*mean(lambda[1:(nyear-1)])/Cit6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:5) {         		                                       # loop
      Sit9[i]<- sum(Cit6[3:6])*mean(lambda[1:(nyear-1)])/sum(Cit6[2:5] )       # Eq (6) from Udevitz 2012
      }   # i
      for (i in 6:(nage-1)) {         		                                       # loop
      Sit9[i]<- sum(Cit6[7:nage])*mean(lambda[1:(nyear-1)])/sum(Cit6[6:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
Sit9<- replace(Sit9,which(is.nan(Sit9)==TRUE),0)                                 # replace impossible value by 0
Sit9<- replace(Sit9,which(is.infinite(Sit9)==TRUE),0)


nbSit9<-3                                                                        # nb of survival parameters for AIC 

Var_Sit9<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit9[i]<- (Sit9[i]^2)*((1-Cit6[i])/Xit[i]+(1-Cit6[i+1])/Xit[i+1])+ mean(VarLam[1:(nyear-1)])*(Cit6[i+1]/Cit6[i])^2                               
      }   # i
      for (i in 2:5) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit9[i]<- (Sit9[i]^2)*((1-sum(Cit6[2:5]))/sum(Xit[2:5])+(1-sum(Cit6[3:6]))/sum(Xit[3:6]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[3:6])/sum(Cit6[2:5]))^2                               
      }   # i
      for (i in 6:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit9[i]<- (Sit9[i]^2)*((1-sum(Cit6[6:(nage-1)]))/sum(Xit[6:(nage-1)])+(1-sum(Cit6[7:nage]))/sum(Xit[7:nage]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[7:nage])/sum(Cit6[6:(nage-1)]))^2                               
      }   # i

Var_Sit9<- replace(Var_Sit9,which(is.nan(Var_Sit9)==TRUE),0)                     # replace impossible value by 0
Var_Sit9<- replace(Var_Sit9,which(is.infinite(Var_Sit9)==TRUE),0)

SD_Sit9 <- sqrt(Var_Sit9)                                                        # matrix of standard deviation of survival estimate 

# table of Lit

Lit9<- matrix(0, nage)  
Lit9[1]<-1
      for (i in 2:(nage-1)) {         		                                       # loop
      Lit9[i]<- prod(Sit9[1:i-1])                                                # Eq (12) from Udevitz 2012
      }   # i

nbLit9<-nage-1                                                                   # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 1998 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

mlb<-rep(0, nage)

for (i in 1:nage) {                                                              # loop
      mlb[i]<-mean(lambda[1:(nyear-1)])^i    
      }   # i   
      
p9<-  rep(0, nage) 
lik9<-  matrix(0, nage, nyear) 

      for (i in 2:nage) {                                                        # loop
      p9[i]<-Lit9[i]/(mlb[i]*sum(Lit9/mlb))    
      }   # i   
      
      p9[1]<-1-sum(p9[2:nage]) 
      
for (t in 1:nyear) {       
      for (i in 1:nage) {                                                        # loop
      lik9[i,t]<-Xit[i,t]*log(p9[i])     
      }   # i   
  }   # t       

lik9<- replace(lik9,which(is.nan(lik9)==TRUE),0)                                 # replace impossible value by 0
lik9<- replace(lik9,which(is.infinite(lik9)==TRUE),0)
loglik9<- colSums (lik9)

# Derive loglikelihood

loglikelihood9 <- fact_sigma_xit-fact_xit+sum(loglik9)

# Derive Aic

AIC9<- 2*(nbSit9+nbLit9)-2*loglikelihood9

## model summary
Sit9
SD_Sit9
Lit9
Cit6
sum(Cit6)
p9
sum(p9)
lik9
loglik9
loglikelihood9
nbSit9
nbLit9
AIC9  


 ### Model 10 Time invariant : assumptions of equal survival between 1;2-5;6-10
# table of Sit

Sit10<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      Sit10[i]<- Cit6[i+1]*mean(lambda[1:(nyear-1)])/Cit6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:5) {         		                                       # loop
      Sit10[i]<- sum(Cit6[4:6])*mean(lambda[1:(nyear-1)])/sum(Cit6[3:5] )       # Eq (6) from Udevitz 2012
      }   # i
      for (i in 6:(nage-1)) {         		                                       # loop
      Sit10[i]<- sum(Cit6[7:nage])*mean(lambda[1:(nyear-1)])/sum(Cit6[6:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
Sit10<- replace(Sit10,which(is.nan(Sit10)==TRUE),0)                                 # replace impossible value by 0
Sit10<- replace(Sit10,which(is.infinite(Sit10)==TRUE),0)


nbSit10<-4                                                                        # nb of survival parameters for AIC 

Var_Sit10<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_Sit10[i]<- (Sit10[i]^2)*((1-Cit6[i])/Xit[i]+(1-Cit6[i+1])/Xit[i+1])+ mean(VarLam[1:(nyear-1)])*(Cit6[i+1]/Cit6[i])^2                               
      }   # i
      for (i in 3:5) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit10[i]<- (Sit10[i]^2)*((1-sum(Cit6[3:5]))/sum(Xit[3:5])+(1-sum(Cit6[4:6]))/sum(Xit[4:6]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[4:6])/sum(Cit6[3:5]))^2                               
      }   # i
      for (i in 6:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Sit10[i]<- (Sit10[i]^2)*((1-sum(Cit6[6:(nage-1)]))/sum(Xit[6:(nage-1)])+(1-sum(Cit6[7:nage]))/sum(Xit[7:nage]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[7:nage])/sum(Cit6[6:(nage-1)]))^2                               
      }   # i

Var_Sit10<- replace(Var_Sit10,which(is.nan(Var_Sit10)==TRUE),0)                     # replace impossible value by 0
Var_Sit10<- replace(Var_Sit10,which(is.infinite(Var_Sit10)==TRUE),0)

SD_Sit10 <- sqrt(Var_Sit10)                                                        # matrix of standard deviation of survival estimate 

# table of Lit

Lit10<- matrix(0, nage)  
Lit10[1]<-1
      for (i in 2:(nage-1)) {         		                                       # loop
      Lit10[i]<- prod(Sit10[1:i-1])                                                # Eq (12) from Udevitz 2012
      }   # i

nbLit10<-nage-1                                                                   # nb of survival parameters for AIC 


### Likelihood estimate for model selection from Eq (4) and (6) of Udevitz 110108 and Eq (1) and (5) from Udevitz 2012

# Derive loglik(t)

mlb<-rep(0, nage)

for (i in 1:nage) {                                                              # loop
      mlb[i]<-mean(lambda[1:(nyear-1)])^i    
      }   # i   
      
p10<-  rep(0, nage) 
lik10<-  matrix(0, nage, nyear) 

      for (i in 2:nage) {                                                        # loop
      p10[i]<-Lit10[i]/(mlb[i]*sum(Lit10/mlb))    
      }   # i   
      
      p10[1]<-1-sum(p10[2:nage]) 
      
for (t in 1:nyear) {       
      for (i in 1:nage) {                                                        # loop
      lik10[i,t]<-Xit[i,t]*log(p10[i])     
      }   # i   
  }   # t       

lik10<- replace(lik10,which(is.nan(lik10)==TRUE),0)                                 # replace impossible value by 0
lik10<- replace(lik10,which(is.infinite(lik10)==TRUE),0)
loglik10<- colSums (lik10)

# Derive loglikelihood

loglikelihood10 <- fact_sigma_xit-fact_xit+sum(loglik10)

# Derive Aic

AIC10<- 2*(nbSit10+nbLit10)-2*loglikelihood10

## model summary
Sit10
SD_Sit10
Lit10
Cit6
sum(Cit6)
p10
sum(p10)
lik10
loglik10
loglikelihood10
nbSit10
nbLit10
AIC10  


### MODEL SELECTION

TotRes<-cbind(c("all","6-10","3-5;6-10","2-5;6-10","2-10","1-10","all; time inv ","2-10; time inv","1-10; time inv","2-5;6-10; time inv","3-5;6-10; time inv"),
              c(nbSit,nbSit1,nbSit2,nbSit3,nbSit4,nbSit5,nbSit6,nbSit7,nbSit8,nbSit9,nbSit10),
              c(nbCit_exp,nbCit_exp1,nbCit_exp2,nbCit_exp3,nbCit_exp4,nbCit_exp5,nbLit6,nbLit7,nbLit8,nbLit9,nbLit10),
              c(loglikelihood, loglikelihood1, loglikelihood2, loglikelihood3, loglikelihood4, loglikelihood5, loglikelihood6, loglikelihood7, loglikelihood8, loglikelihood9, loglikelihood10),
              c(AIC,AIC1,AIC2,AIC3,AIC4,AIC5,AIC6,AIC7,AIC8,AIC9,AIC10))
TotRes<-as.data.frame(TotRes)
names(TotRes )<-c("Model","nb Si","nb Ci","logL", "AIC")
TotRes

Sit5   #best model for annual estimation = one survival rate for all adult class
Sit7   #best model for constant estimation = yearling and adult estimation

# plot results

plot(c(2003,2004,2005,2006),t(Sit5[1,]),type="b",pch=16,ylim=c(0,1),col=c("black"),lwd=2,cex=2,ylab="Survival probalility with SD",xlab="Year")
  arrows(c(2003,2004,2005,2006),Sit5[1,]-SD_Sit5[1,],c(2003,2004,2005,2006), Sit5[1,]+SD_Sit5[1,], angle=90, col="black",code=3,length=0.1)
 
title("Survival Probability in FDC35")
  legend(x = 2005, y = 1, legend = c("All data", "Trapping", "Hunting"), pch=16:18,lty = c(1, 1), lwd = c(2, 2), col=c("black","grey50","grey50"), bty = "n", cex = 1)


############################################
#### Sampling method influence on estimation
###########################################"

### survival estimate between 2002 and 2007

Hunt<-subset(FOX35,FOX35$Mode=="C")          
Ch<-table(Hunt$Age,Hunt$Year)
Trap<-subset(FOX35,FOX35$Mode=="P")            
Tr<-table(Trap$Age,Trap$Year)

Xit<-Ch[,3:7]                               ## change it to choice the sampling method !
nyear<- ncol(Xit)
nage<-nrow(Xit)

# table of Cit

Cit<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      Cit[i,t]<- Xit[i,t]/sum(Xit[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### Model5:  Assumes survival rates are equal for ages 1-8

# table of Chit

Chit5<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
       
      for (i in 1:(nage-1)) {         		                                       # loop
      Chit5[i,t]<- sum(Cit[2:nage,t+1])*lambda[t]/sum(Cit[1:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

Chit5<- replace(Chit5,which(is.nan(Chit5)==TRUE),0)                                 # replace impossible value by 0
Chit5<- replace(Chit5,which(is.infinite(Chit5)==TRUE),0)

Var_Chit5<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
     
      for (i in 1:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Chit5[i,t]<- (Chit5[i,t]^2)*((1-sum(Cit[1:(nage-1),t]))/sum(Xit[1:(nage-1),t])+(1-sum(Cit[2:nage,t+1]))/sum(Xit[2:nage,t+1]))+ VarLam[t]*(sum(Cit[2:nage,t+1])/sum(Cit[1:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_Chit5<- replace(Var_Chit5,which(is.nan(Var_Chit5)==TRUE),0)                     # replace impossible value by 0
Var_Chit5<- replace(Var_Chit5,which(is.infinite(Var_Chit5)==TRUE),0)

SD_Chit5 <- sqrt(Var_Chit5)                                                        # matrix of standard deviation of survival estimate 

### Model 7 Time invariant : assumptions of equal survival between 2-8
# table of Chit

Chit7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      Chit7[i]<- Cit6[i+1]*mean(lambda[1:(nyear-1)])/Cit6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      Chit7[i]<- sum(Cit6[3:nage])*mean(lambda[1:(nyear-1)])/sum(Cit6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
Chit7<- replace(Chit7,which(is.nan(Chit7)==TRUE),0)                                 # replace impossible value by 0
Chit7<- replace(Chit7,which(is.infinite(Chit7)==TRUE),0)


Var_Chit7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_Chit7[i]<- (Chit7[i]^2)*((1-Cit6[i])/Xit[i]+(1-Cit6[i+1])/Xit[i+1])+ mean(VarLam[1:(nyear-1)])*(Cit6[i+1]/Cit6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_Chit7[i]<- (Chit7[i]^2)*((1-sum(Cit6[2:(nage-1)]))/sum(Xit[2:(nage-1)])+(1-sum(Cit6[3:nage]))/sum(Xit[3:nage]))+ mean(VarLam[1:(nyear-1)])*(sum(Cit6[3:nage])/sum(Cit6[2:(nage-1)]))^2                               
      }   # i

Var_Chit7<- replace(Var_Chit7,which(is.nan(Var_Chit7)==TRUE),0)                     # replace impossible value by 0
Var_Chit7<- replace(Var_Chit7,which(is.infinite(Var_Chit7)==TRUE),0)

SD_Chit7 <- sqrt(Var_Chit7)                                                        # matrix of standard deviation of survival estimate 

Chit5
Chit7

# plot results

points(c(2003,2004,2005,2006),t(Chit5[1,]),type="b",pch=18,ylim=c(0,1),col=c("grey50"),lwd=2,cex=2,ylab="Survival probalility with SD",xlab="Year")
  arrows(c(2003,2004,2005,2006),Chit5[1,]-SD_Chit5[1,],c(2003,2004,2005,2006), Chit5[1,]+SD_Chit5[1,], angle=90, col="grey50",code=3,length=0.1)
 


#############################
#### Site by Site estimation
#############################

names(FOX35)

### Construct life table of standing age structure

dom<-as.matrix(table(FOX35$Age[FOX35$GIC=="D"],FOX35$Year[FOX35$GIC=="D"]))
fou<-as.matrix(table(FOX35$Age[FOX35$GIC=="F"],FOX35$Year[FOX35$GIC=="F"]))
ven<-as.matrix(table(FOX35$Age[FOX35$GIC=="V"],FOX35$Year[FOX35$GIC=="V"]))

###########################################
###  DOMAGNE GIC   entre 2002 et 2007
###########################################

dom
lamD
VarLamD

### other data
DensD<- DENS$DENSITY[DENS$GIC=="DOM"]
CVDD<-DENS$DCV[DENS$GIC=="DOM"]/100


LSD<-DENS[1:6,36:39]                          # Litter size estimate
rownames(LSD)<-c("2002","2003","2004","2005","2006","2007")

PBD<-DENS[1:6,30:32]                        # Probability of breeding
rownames(PBD)<-c("2002","2003","2004","2005","2006","2007")


### ESTIMATE X0F, number of newborns per year based on the number of female caugth
#################

domF<-as.matrix(table(FOX35$Age[FOX35$GIC=="D"&FOX35$Sex=="F"],FOX35$Year[FOX35$GIC=="D"&FOX35$Sex=="F"]))      #nb of female caugth
domF<-domF[,2:7]

X0D<-rep(NA,6)

for (t in 1:length(X0D)) {
  X0D1<- domF[1,t]*PBD[t,1]*LSD[t,1] 
  X0D2<- domF[2,t]*PBD[t,2]*LSD[t,2] 
  X0D35<- sum(domF[3:5,t])*PBD[t,2]*LSD[t,3] 
  X0D610<-sum(domF[6:10,t])*PBD[t,3]*LSD[t,4] 
  
  X0D[t]<- X0D1+X0D2+X0D35+X0D610
      }   # t
      
X0D<-round(X0D)

### survival estimate between 2002 and 2007

XitD<-rbind(X0D,dom[,2:7] )                                                        # add X0 to life table !
nyear<- ncol(XitD)
nage<-nrow(XitD)

X0D/colSums(XitD)
# table of Cit

CitD<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitD[i,t]<- XitD[i,t]/sum(XitD[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

# Derive factorial

fact_xitD <- lfactorial(XitD)
fact_xitD <- sum(fact_xitD)

fact_sigma_xitD <- sum(lfactorial(colSums(XitD)))


### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitD<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitD[1,t]<- CitD[2,t+1]*lamD[t]/CitD[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitD[i,t]<- sum(CitD[3:nage,t+1])*lamD[t]/sum(CitD[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitD<- replace(SitD,which(is.nan(SitD)==TRUE),0)                                 # replace impossible value by 0
SitD<- replace(SitD,which(is.infinite(SitD)==TRUE),0)

nbSitD<-10  
                                                                   # nb of survival parameters for AIC 
Var_SitD<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitD[1,t]<- (SitD[1,t]^2)*((1-CitD[1,t])/XitD[1,t]+(1-CitD[2,t+1])/XitD[2,t+1])+ VarLamD[t]*(CitD[2,t+1]/CitD[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitD[i,t]<- (SitD[i,t]^2)*((1-sum(CitD[2:(nage-1),t]))/sum(XitD[2:(nage-1),t])+(1-sum(CitD[3:nage,t+1]))/sum(XitD[3:nage,t+1]))+ VarLamD[t]*(sum(CitD[3:nage,t+1])/sum(CitD[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitD<- replace(Var_SitD,which(is.nan(Var_SitD)==TRUE),0)                     # replace impossible value by 0
Var_SitD<- replace(Var_SitD,which(is.infinite(Var_SitD)==TRUE),0)

SD_SitD <- sqrt(Var_SitD)                                                        # matrix of standard deviation of survival estimate 


  
##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitD6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitD6[i]<- sum(XitD[i,])/sum(XitD)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitD7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitD7[i]<- CitD6[i+1]*mean(lamD[1:(nyear-1)])/CitD6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitD7[i]<- sum(CitD6[3:nage])*mean(lamD[1:(nyear-1)])/sum(CitD6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitD7<- replace(SitD7,which(is.nan(SitD7)==TRUE),0)                                 # replace impossible value by 0
SitD7<- replace(SitD7,which(is.infinite(SitD7)==TRUE),0)

nbSitD7<-2

Var_SitD7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitD7[i]<- (SitD7[i]^2)*((1-CitD6[i])/XitD[i]+(1-CitD6[i+1])/XitD[i+1])+ mean(VarLamD[1:(nyear-1)])*(CitD6[i+1]/CitD6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitD7[i]<- (SitD7[i]^2)*((1-sum(CitD6[2:(nage-1)]))/sum(XitD[2:(nage-1)])+(1-sum(CitD6[3:nage]))/sum(XitD[3:nage]))+ mean(VarLamD[1:(nyear-1)])*(sum(CitD6[3:nage])/sum(CitD6[2:(nage-1)]))^2                               
      }   # i

Var_SitD7<- replace(Var_SitD7,which(is.nan(Var_SitD7)==TRUE),0)                     # replace impossible value by 0
Var_SitD7<- replace(Var_SitD7,which(is.infinite(Var_SitD7)==TRUE),0)

SD_SitD7 <- sqrt(Var_SitD7)                                                        # matrix of standard deviation of survival estimate 


##### plot results

DomSit<- rbind(SitD[1:2,],SD_SitD[1:2,])
rownames(DomSit)<-c("JuvenileS0F","Adults", "SD_S0","SD_SA")
DS<-cbind(DomSit,c(SitD7[1:2],SD_SitD7[1:2]))
colnames(DS)<-c("2002","2003","2004","2005","2006","Average")


### ESTIMATE X0T, number of newborns per year based on the total number of removals and assuming a 1:1 sex ratio
#################

domF<-dom[,2:7]*0.5     #nb of female caugth = 0.5* total number of removals


X0D<-rep(NA,6)

for (t in 1:length(X0D)) {
  X0D1<- domF[1,t]*PBD[t,1]*LSD[t,1] 
  X0D2<- domF[2,t]*PBD[t,2]*LSD[t,2] 
  X0D35<- sum(domF[3:5,t])*PBD[t,2]*LSD[t,3] 
  X0D610<-sum(domF[6:10,t])*PBD[t,3]*LSD[t,4] 
  
  X0D[t]<- X0D1+X0D2+X0D35+X0D610
      }   # t
      
X0D<-round(X0D)

### survival estimate between 2002 and 2007

XitD<-rbind(X0D,dom[,2:7] )                                                        # add X0 to life table !
nyear<- ncol(XitD)
nage<-nrow(XitD)

X0D/colSums(XitD)

# table of Cit
CitD<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitD[i,t]<- XitD[i,t]/sum(XitD[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitD<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitD[1,t]<- CitD[2,t+1]*lamD[t]/CitD[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitD[i,t]<- sum(CitD[3:nage,t+1])*lamD[t]/sum(CitD[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitD<- replace(SitD,which(is.nan(SitD)==TRUE),0)                                 # replace impossible value by 0
SitD<- replace(SitD,which(is.infinite(SitD)==TRUE),0)

nbSitD<-10  
                                                                   # nb of survival parameters for AIC 
Var_SitD<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitD[1,t]<- (SitD[1,t]^2)*((1-CitD[1,t])/XitD[1,t]+(1-CitD[2,t+1])/XitD[2,t+1])+ VarLamD[t]*(CitD[2,t+1]/CitD[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitD[i,t]<- (SitD[i,t]^2)*((1-sum(CitD[2:(nage-1),t]))/sum(XitD[2:(nage-1),t])+(1-sum(CitD[3:nage,t+1]))/sum(XitD[3:nage,t+1]))+ VarLamD[t]*(sum(CitD[3:nage,t+1])/sum(CitD[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitD<- replace(Var_SitD,which(is.nan(Var_SitD)==TRUE),0)                     # replace impossible value by 0
Var_SitD<- replace(Var_SitD,which(is.infinite(Var_SitD)==TRUE),0)

SD_SitD <- sqrt(Var_SitD)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitD6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitD6[i]<- sum(XitD[i,])/sum(XitD)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitD7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitD7[i]<- CitD6[i+1]*mean(lamD[1:(nyear-1)])/CitD6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitD7[i]<- sum(CitD6[3:nage])*mean(lamD[1:(nyear-1)])/sum(CitD6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitD7<- replace(SitD7,which(is.nan(SitD7)==TRUE),0)                                 # replace impossible value by 0
SitD7<- replace(SitD7,which(is.infinite(SitD7)==TRUE),0)

nbSitD7<-2

Var_SitD7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitD7[i]<- (SitD7[i]^2)*((1-CitD6[i])/XitD[i]+(1-CitD6[i+1])/XitD[i+1])+ mean(VarLamD[1:(nyear-1)])*(CitD6[i+1]/CitD6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitD7[i]<- (SitD7[i]^2)*((1-sum(CitD6[2:(nage-1)]))/sum(XitD[2:(nage-1)])+(1-sum(CitD6[3:nage]))/sum(XitD[3:nage]))+ mean(VarLamD[1:(nyear-1)])*(sum(CitD6[3:nage])/sum(CitD6[2:(nage-1)]))^2                               
      }   # i

Var_SitD7<- replace(Var_SitD7,which(is.nan(Var_SitD7)==TRUE),0)                     # replace impossible value by 0
Var_SitD7<- replace(Var_SitD7,which(is.infinite(Var_SitD7)==TRUE),0)

SD_SitD7 <- sqrt(Var_SitD7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitD8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitD8[i]<- CitD6[i+1]*mean(lamD[1:(nyear-1)])/CitD6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitD8[i]<- sum(CitD6[4:nage])*mean(lamD[1:(nyear-1)])/sum(CitD6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitD8<- replace(SitD8,which(is.nan(SitD8)==TRUE),0)                                 # replace impossible value by 0
SitD8<- replace(SitD8,which(is.infinite(SitD8)==TRUE),0)

nbSitD8<-2

Var_SitD8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitD8[i]<- (SitD8[i]^2)*((1-CitD6[i])/XitD[i]+(1-CitD6[i+1])/XitD[i+1])+ mean(VarLamD[1:(nyear-1)])*(CitD6[i+1]/CitD6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitD8[i]<- (SitD8[i]^2)*((1-sum(CitD6[3:(nage-1)]))/sum(XitD[3:(nage-1)])+(1-sum(CitD6[4:nage]))/sum(XitD[4:nage]))+ mean(VarLamD[1:(nyear-1)])*(sum(CitD6[4:nage])/sum(CitD6[3:(nage-1)]))^2                               
      }   # i

Var_SitD8<- replace(Var_SitD8,which(is.nan(Var_SitD8)==TRUE),0)                     # replace impossible value by 0
Var_SitD8<- replace(Var_SitD8,which(is.infinite(Var_SitD8)==TRUE),0)

SD_SitD8 <- sqrt(Var_SitD8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

DomSitT<- rbind(SitD[1:2,],SD_SitD[1:2,])
rownames(DomSitT)<-c("JuvenileS0F","Adults", "SD_S0","SD_SA")
DST<-cbind(DomSitT,c(SitD7[1:2],SD_SitD7[1:2]))
colnames(DST)<-c("2002","2003","2004","2005","2006","Average")


plot(c(2002,2003,2004,2005,2006),DST[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2002,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2002.1,2003.1,2004.1,2005.1,2006.1),DST[2,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2002,2003,2004,2005,2006),DST[1,1:5]-DST[3,1:5],c(2002,2003,2004,2005,2006), DST[1,1:5]+DST[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2002.1,2003.1,2004.1,2005.1,2006.1),DST[2,1:5]-DST[4,1:5],c(2002.1,2003.1,2004.1,2005.1,2006.1), DST[2,1:5]+DST[4,1:5], angle=90, col="black",code=3,length=0)

  
  legend(x =2002, y = 1.05, legend = c("Juveniles : 0.28 +/- 0.05", "Adults : 0.51 +/- 0.10"), pch=16:18,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1,"Domagne",font=2)

# Plot difference between juvenile survival

plot(c(2002,2003,2004,2005,2006),DS[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2002,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2002.1,2003.1,2004.1,2005.1,2006.1),DST[1,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2002,2003,2004,2005,2006),DS[1,1:5]-DS[3,1:5],c(2002,2003,2004,2005,2006), DS[1,1:5]+DS[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2002.1,2003.1,2004.1,2005.1,2006.1),DST[1,1:5]-DST[3,1:5],c(2002.1,2003.1,2004.1,2005.1,2006.1), DST[1,1:5]+DST[3,1:5], angle=90, col="black",code=3,length=0)

  
  legend(x =2002, y = 1.05, legend = c("S0F: 0.29 +/- 0.06","S0T : 0.28 +/- 0.05"), pch=16:17,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1,"Domagne",font=2)

  #three age class model
  
TCMD<- cbind(SitD8[1:3],SD_SitD8[1:3])       
TCMD
colnames(TCMD)<-c("Phi","SD_Phi")
rownames(TCMD)<-c("SJ","SY","SA")

plot(c(1,2,3),TCMD[,1],pch=16,cex=1.5,ylim=c(0,1),xaxt="n",xlim=c(0.5,3.5),ylab="Survival probalility with SD",xlab="")
    arrows(c(1,2,3),TCMD[,1]-TCMD[,2],c(1,2,3),TCMD[,1]+TCMD[,2], angle=90, col="black",code=3,length=0)
axis(1, at = 1, lab = expression("Juvenile"))
axis(1, at = 2, lab = expression("Yearling")) 
axis(1, at = 3, lab = expression("Adult")) 
text (x=2, y=1,"Domagne",font=2)  


###########################################
###  Vendelais GIC   entre 2003 et 2007
###########################################

ven <-ven [,3:7]
lamV <-lamV[2:6]
VarLamV <-VarLamV[2:6]


### other data
DensV<- DENS$DENSITY[DENS$GIC=="VEN"]
CVDV<-DENS$DCV[DENS$GIC=="VEN"]/100

LSV<-DENS[12:16,36:39]                        # Litter size estimate
rownames(LSV)<-c("2003","2004","2005","2006","2007")

PBV<-DENS[12:16,30:32]                    # Probability of breeding
rownames(PBV)<-c("2003","2004","2005","2006","2007")

### ESTIMATE X0F, number of newborns per year based on the number of female caugth
#################

venF<-as.matrix(table(FOX35$Age[FOX35$GIC=="V"&FOX35$Sex=="F"],FOX35$Year[FOX35$GIC=="V"&FOX35$Sex=="F"]))      #nb of female caugth
venF<-venF[,3:7]

X0V<-rep(NA,5)

for (t in 1:length(X0V)) {
  X0V1<- venF[1,t]*PBV[t,1]*LSV[t,1] 
  X0V2<- venF[2,t]*PBV[t,2]*LSV[t,2] 
  X0V35<- sum(venF[3:5,t])*PBV[t,2]*LSV[t,3] 
  X0V610<-sum(venF[6:10,t])*PBV[t,3]*LSV[t,4] 
  
  X0V[t]<- X0V1+X0V2+X0V35+X0V610
      }   # t
      
X0V<-round(X0V)

### survival estimate between 2002 and 2007

XitV<-rbind(X0V,ven )                                                        # add X0 to life table !
nyear<- ncol(XitV)
nage<-nrow(XitV)

X0V/colSums(XitV)

# table of Cit
CitV<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitV[i,t]<- XitV[i,t]/sum(XitV[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitV<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitV[1,t]<- CitV[2,t+1]*lamV[t]/CitV[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitV[i,t]<- sum(CitV[3:nage,t+1])*lamV[t]/sum(CitV[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitV<- replace(SitV,which(is.nan(SitV)==TRUE),0)                                 # replace impossible value by 0
SitV<- replace(SitV,which(is.infinite(SitV)==TRUE),0)

Var_SitV<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitV[1,t]<- (SitV[1,t]^2)*((1-CitV[1,t])/XitV[1,t]+(1-CitV[2,t+1])/XitV[2,t+1])+ VarLamV[t]*(CitV[2,t+1]/CitV[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitV[i,t]<- (SitV[i,t]^2)*((1-sum(CitV[2:(nage-1),t]))/sum(XitV[2:(nage-1),t])+(1-sum(CitV[3:nage,t+1]))/sum(XitV[3:nage,t+1]))+ VarLamV[t]*(sum(CitV[3:nage,t+1])/sum(CitV[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitV<- replace(Var_SitV,which(is.nan(Var_SitV)==TRUE),0)                     # replace impossible value by 0
Var_SitV<- replace(Var_SitV,which(is.infinite(Var_SitV)==TRUE),0)

SD_SitV <- sqrt(Var_SitV)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitV6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitV6[i]<- sum(XitV[i,])/sum(XitV)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitV7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitV7[i]<- CitV6[i+1]*mean(lamV[1:(nyear-1)])/CitV6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitV7[i]<- sum(CitV6[3:nage])*mean(lamV[1:(nyear-1)])/sum(CitV6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitV7<- replace(SitV7,which(is.nan(SitV7)==TRUE),0)                                 # replace impossible value by 0
SitV7<- replace(SitV7,which(is.infinite(SitV7)==TRUE),0)


Var_SitV7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitV7[i]<- (SitV7[i]^2)*((1-CitV6[i])/XitV[i]+(1-CitV6[i+1])/XitV[i+1])+ mean(VarLamV[1:(nyear-1)])*(CitV6[i+1]/CitV6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitV7[i]<- (SitV7[i]^2)*((1-sum(CitV6[2:(nage-1)]))/sum(XitV[2:(nage-1)])+(1-sum(CitV6[3:nage]))/sum(XitV[3:nage]))+ mean(VarLamV[1:(nyear-1)])*(sum(CitV6[3:nage])/sum(CitV6[2:(nage-1)]))^2                               
      }   # i

Var_SitV7<- replace(Var_SitV7,which(is.nan(Var_SitV7)==TRUE),0)                     # replace impossible value by 0
Var_SitV7<- replace(Var_SitV7,which(is.infinite(Var_SitV7)==TRUE),0)

SD_SitV7 <- sqrt(Var_SitV7)                                                        # matrix of standard deviation of survival estimate 


##### plot results

VenSit<- rbind( SitV[1:2,],SD_SitV[1:2,])
rownames(VenSit)<-c("JuvenileS0F","Adults","SD_S0", "SD_SA")

VS<-cbind(VenSit,c(SitV7[1:2],SD_SitV7[1:2]))
colnames(VS)<-c("2003","2004","2005","2006","Average")


### ESTIMATE X0T, number of newborns per year based on the total number of removals and assuming a 1:1 sex ratio
#################

venF<- ven*0.5     #nb of female caugth = 0.5* total number of removals

X0V<-rep(NA,5)

for (t in 1:length(X0V)) {
  X0V1<- venF[1,t]*PBV[t,1]*LSV[t,1] 
  X0V2<- venF[2,t]*PBV[t,2]*LSV[t,2] 
  X0V35<- sum(venF[3:5,t])*PBV[t,2]*LSV[t,3] 
  X0V610<-sum(venF[6:10,t])*PBV[t,3]*LSV[t,4] 
  
  X0V[t]<- X0V1+X0V2+X0V35+X0V610
      }   # t
      
X0V<-round(X0V)

### survival estimate between 2002 and 2007

XitV<-rbind(X0V,ven )                                                        # add X0 to life table !
nyear<- ncol(XitV)
nage<-nrow(XitV)

X0V/colSums(XitV)
# table of Cit

# table of Cit
CitV<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitV[i,t]<- XitV[i,t]/sum(XitV[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitV<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitV[1,t]<- CitV[2,t+1]*lamV[t]/CitV[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitV[i,t]<- sum(CitV[3:nage,t+1])*lamV[t]/sum(CitV[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitV<- replace(SitV,which(is.nan(SitV)==TRUE),0)                                 # replace impossible value by 0
SitV<- replace(SitV,which(is.infinite(SitV)==TRUE),0)

Var_SitV<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitV[1,t]<- (SitV[1,t]^2)*((1-CitV[1,t])/XitV[1,t]+(1-CitV[2,t+1])/XitV[2,t+1])+ VarLamV[t]*(CitV[2,t+1]/CitV[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitV[i,t]<- (SitV[i,t]^2)*((1-sum(CitV[2:(nage-1),t]))/sum(XitV[2:(nage-1),t])+(1-sum(CitV[3:nage,t+1]))/sum(XitV[3:nage,t+1]))+ VarLamV[t]*(sum(CitV[3:nage,t+1])/sum(CitV[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitV<- replace(Var_SitV,which(is.nan(Var_SitV)==TRUE),0)                     # replace impossible value by 0
Var_SitV<- replace(Var_SitV,which(is.infinite(Var_SitV)==TRUE),0)

SD_SitV <- sqrt(Var_SitV)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitV6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitV6[i]<- sum(XitV[i,])/sum(XitV)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitV7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitV7[i]<- CitV6[i+1]*mean(lamV[1:(nyear-1)])/CitV6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitV7[i]<- sum(CitV6[3:nage])*mean(lamV[1:(nyear-1)])/sum(CitV6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitV7<- replace(SitV7,which(is.nan(SitV7)==TRUE),0)                                 # replace impossible value by 0
SitV7<- replace(SitV7,which(is.infinite(SitV7)==TRUE),0)


Var_SitV7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitV7[i]<- (SitV7[i]^2)*((1-CitV6[i])/XitV[i]+(1-CitV6[i+1])/XitV[i+1])+ mean(VarLamV[1:(nyear-1)])*(CitV6[i+1]/CitV6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitV7[i]<- (SitV7[i]^2)*((1-sum(CitV6[2:(nage-1)]))/sum(XitV[2:(nage-1)])+(1-sum(CitV6[3:nage]))/sum(XitV[3:nage]))+ mean(VarLamV[1:(nyear-1)])*(sum(CitV6[3:nage])/sum(CitV6[2:(nage-1)]))^2                               
      }   # i

Var_SitV7<- replace(Var_SitV7,which(is.nan(Var_SitV7)==TRUE),0)                     # replace impossible value by 0
Var_SitV7<- replace(Var_SitV7,which(is.infinite(Var_SitV7)==TRUE),0)

SD_SitV7 <- sqrt(Var_SitV7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitV8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitV8[i]<- CitV6[i+1]*mean(lamV[1:(nyear-1)])/CitV6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitV8[i]<- sum(CitV6[4:nage])*mean(lamV[1:(nyear-1)])/sum(CitV6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitV8<- replace(SitV8,which(is.nan(SitV8)==TRUE),0)                                 # replace impossible value by 0
SitV8<- replace(SitV8,which(is.infinite(SitV8)==TRUE),0)


Var_SitV8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitV8[i]<- (SitV8[i]^2)*((1-CitV6[i])/XitV[i]+(1-CitV6[i+1])/XitV[i+1])+ mean(VarLamV[1:(nyear-1)])*(CitV6[i+1]/CitV6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitV8[i]<- (SitV8[i]^2)*((1-sum(CitV6[3:(nage-1)]))/sum(XitV[3:(nage-1)])+(1-sum(CitV6[4:nage]))/sum(XitV[4:nage]))+ mean(VarLamV[1:(nyear-1)])*(sum(CitV6[4:nage])/sum(CitV6[3:(nage-1)]))^2                               
      }   # i

Var_SitV8<- replace(Var_SitV8,which(is.nan(Var_SitV8)==TRUE),0)                     # replace impossible value by 0
Var_SitV8<- replace(Var_SitV8,which(is.infinite(Var_SitV8)==TRUE),0)

SD_SitV8 <- sqrt(Var_SitV8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

VenSitT<- rbind( SitV[1:2,],SD_SitV[1:2,])
rownames(VenSitT)<-c("JuvenileS0F","Adults","SD_S0", "SD_SA")

VST<-cbind(VenSitT,c(SitV7[1:2],SD_SitV7[1:2]))
colnames(VST)<-c("2003","2004","2005","2006","Average")

plot(c(2003,2004,2005,2006),VST[1,1:4],type="b",pch=16,ylim=c(0,1.5),xaxt="n",xlim=c(2003,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2003.1,2004.1,2005.1,2006.1),VST[2,1:4],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2003,2004,2005,2006),VST[1,1:4]-VST[3,1:4],c(2003,2004,2005,2006), VST[1,1:4]+VST[3,1:4], angle=90, col="black",code=3,length=0)
  arrows(c(2003.1,2004.1,2005.1,2006.1),VST[2,1:4]-VST[4,1:4],c(2003.1,2004.1,2005.1,2006.1), VST[2,1:4]+VST[4,1:4], angle=90, col="black",code=3,length=0)

axis(1, at = 2003, lab = expression("2003"))
axis(1, at = 2004, lab = expression("2004")) 
axis(1, at = 2005, lab = expression("2005")) 
axis(1, at = 2006, lab = expression("2006")) 

  legend(x =2003, y = 1.55, legend = c("Juveniles : 0.27 +/- 0.07", "Adults : 0.42 +/- 0.12"), pch=16:18,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1.4,"Vendelais",font=2)

# Plot difference between juvenile survival

plot(c(2003,2004,2005,2006),VS[1,1:4],type="b",pch=16,ylim=c(0,1),xaxt="n",xlim=c(2003,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2003.1,2004.1,2005.1,2006.1),VST[1,1:4],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2003,2004,2005,2006),VS[1,1:4]-VS[3,1:4],c(2003,2004,2005,2006), VS[1,1:4]+VS[3,1:4], angle=90, col="black",code=3,length=0)
  arrows(c(2003.1,2004.1,2005.1,2006.1),VST[1,1:4]-VST[3,1:4],c(2003.1,2004.1,2005.1,2006.1), VST[1,1:4]+VST[3,1:4], angle=90, col="black",code=3,length=0)
axis(1, at = 2003, lab = expression("2003"))
axis(1, at = 2004, lab = expression("2004")) 
axis(1, at = 2005, lab = expression("2005")) 
axis(1, at = 2006, lab = expression("2006")) 

  legend(x =2003, y = 1.05, legend = c("S0F: 0.28 +/- 0.07","S0T : 0.27 +/- 0.07"), pch=16:17,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1,"Vendelais",font=2)

        # three age class constant model
   
TCMV<- cbind(SitV8[1:3],SD_SitV8[1:3])       
TCMV
colnames(TCMV)<-c("Phi","SD_Phi")
rownames(TCMV)<-c("SJ","SY","SA")

plot(c(1,2,3),TCMV[,1],pch=16,cex=1.5,ylim=c(0,1),xaxt="n",xlim=c(0.5,3.5),ylab="Survival probalility with SD",xlab="")
    arrows(c(1,2,3),TCMV[,1]-TCMV[,2],c(1,2,3),TCMV[,1]+TCMV[,2], angle=90, col="black",code=3,length=0)
axis(1, at = 1, lab = expression("Juvenile"))
axis(1, at = 2, lab = expression("Yearling")) 
axis(1, at = 3, lab = expression("Adult")) 
text (x=2, y=1,"Vendelais",font=2)  


###########################################
###  FOUGERES GIC   entre 2002 et 2007
###########################################

fou  <- fou[,3:7]
lamF  <- lamF[2:6]
VarLamF  <- VarLamF[2:6]


###other data
DensF<- DENS$DENSITY[DENS$GIC=="FOU"]
CVDF<-DENS$DCV[DENS$GIC=="FOU"]/100

LSF<-DENS[7:11,36:39]                           # Litter size estimate
rownames(LSF)<-c("2003","2004","2005","2006","2007")

PBF<-DENS[7:11,30:32]                      # Probability of breeding
rownames(PBF)<-c("2003","2004","2005","2006","2007")

### ESTIMATE X0F, number of newborns per year based on the number of female caugth
#################

fouF<-as.matrix(table(FOX35$Age[FOX35$GIC=="F"&FOX35$Sex=="F"],FOX35$Year[FOX35$GIC=="F"&FOX35$Sex=="F"]))      #nb of female caugth
fouF<-fouF[,2:7]

X0F<-rep(NA,5)

for (t in 1:length(X0F)) {
  X0F1<- fouF[1,t]*PBF[t,1]*LSF[t,1] 
  X0F2<- fouF[2,t]*PBF[t,2]*LSF[t,2] 
  X0F35<- sum(fouF[3:5,t])*PBF[t,2]*LSF[t,3] 
  X0F610<-sum(fouF[6:10,t])*PBF[t,3]*LSF[t,4] 
  
  X0F[t]<- X0F1+X0F2+X0F35+X0F610
      }   # t
      
X0F<-round(X0F)
### survival estimate between 2002 and 2007

XitF<-rbind(X0F,fou)                                                        # add X0 to life table !
nyear<- ncol(XitF)
nage<-nrow(XitF)

X0F/colSums(XitF)

# table of Cit
CitF<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitF[i,t]<- XitF[i,t]/sum(XitF[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t


### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitF<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitF[1,t]<- CitF[2,t+1]*lamF[t]/CitF[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitF[i,t]<- sum(CitF[3:nage,t+1])*lamF[t]/sum(CitF[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitF<- replace(SitF,which(is.nan(SitF)==TRUE),0)                                 # replace impossible value by 0
SitF<- replace(SitF,which(is.infinite(SitF)==TRUE),0)

Var_SitF<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitF[1,t]<- (SitF[1,t]^2)*((1-CitF[1,t])/XitF[1,t]+(1-CitF[2,t+1])/XitF[2,t+1])+ VarLamF[t]*(CitF[2,t+1]/CitF[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitF[i,t]<- (SitF[i,t]^2)*((1-sum(CitF[2:(nage-1),t]))/sum(XitF[2:(nage-1),t])+(1-sum(CitF[3:nage,t+1]))/sum(XitF[3:nage,t+1]))+ VarLamF[t]*(sum(CitF[3:nage,t+1])/sum(CitF[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitF<- replace(Var_SitF,which(is.nan(Var_SitF)==TRUE),0)                     # replace impossible value by 0
Var_SitF<- replace(Var_SitF,which(is.infinite(Var_SitF)==TRUE),0)

SD_SitF <- sqrt(Var_SitF)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitF6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitF6[i]<- sum(XitF[i,])/sum(XitF)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitF7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitF7[i]<- CitF6[i+1]*mean(lamF[1:(nyear-1)])/CitF6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitF7[i]<- sum(CitF6[3:nage])*mean(lamF[1:(nyear-1)])/sum(CitF6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitF7<- replace(SitF7,which(is.nan(SitF7)==TRUE),0)                                 # replace impossible value by 0
SitF7<- replace(SitF7,which(is.infinite(SitF7)==TRUE),0)


Var_SitF7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitF7[i]<- (SitF7[i]^2)*((1-CitF6[i])/XitF[i]+(1-CitF6[i+1])/XitF[i+1])+ mean(VarLamF[1:(nyear-1)])*(CitF6[i+1]/CitF6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitF7[i]<- (SitF7[i]^2)*((1-sum(CitF6[2:(nage-1)]))/sum(XitF[2:(nage-1)])+(1-sum(CitF6[3:nage]))/sum(XitF[3:nage]))+ mean(VarLamF[1:(nyear-1)])*(sum(CitF6[3:nage])/sum(CitF6[2:(nage-1)]))^2                               
      }   # i

Var_SitF7<- replace(Var_SitF7,which(is.nan(Var_SitF7)==TRUE),0)                     # replace impossible value by 0
Var_SitF7<- replace(Var_SitF7,which(is.infinite(Var_SitF7)==TRUE),0)

SD_SitF7 <- sqrt(Var_SitF7)                                                        # matrix of standard deviation of survival estimate 


##### plot results

FouSit<- rbind( SitF[1:2,],SD_SitF[1:2,])
rownames(FouSit)<-c("JuvenileS0F","Adults", "SD_S0","SD_SA")

FS<-cbind(FouSit,c(SitF7[1:2],SD_SitF7[1:2]))
colnames(FS)<-c("2003","2004","2005","2006","Average")


 ### ESTIMATE X0T, number of newborns per year based on the number of total removal assuming a 1:1 sex ratio
#################

fouF<- fou*0.5

X0F<-rep(NA,5)

for (t in 1:length(X0F)) {
  X0F1<- fouF[1,t]*PBF[t,1]*LSF[t,1] 
  X0F2<- fouF[2,t]*PBF[t,2]*LSF[t,2] 
  X0F35<- sum(fouF[3:5,t])*PBF[t,2]*LSF[t,3] 
  X0F610<-sum(fouF[6:10,t])*PBF[t,3]*LSF[t,4] 
  
  X0F[t]<- X0F1+X0F2+X0F35+X0F610
      }   # t
      
X0F<-round(X0F)
### survival estimate between 2002 and 2007

XitF<-rbind(X0F,fou)                                                        # add X0 to life table !
nyear<- ncol(XitF)
nage<-nrow(XitF)

X0F/colSums(XitF)

# table of Cit
CitF<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitF[i,t]<- XitF[i,t]/sum(XitF[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t


### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitF<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitF[1,t]<- CitF[2,t+1]*lamF[t]/CitF[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitF[i,t]<- sum(CitF[3:nage,t+1])*lamF[t]/sum(CitF[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitF<- replace(SitF,which(is.nan(SitF)==TRUE),0)                                 # replace impossible value by 0
SitF<- replace(SitF,which(is.infinite(SitF)==TRUE),0)

Var_SitF<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitF[1,t]<- (SitF[1,t]^2)*((1-CitF[1,t])/XitF[1,t]+(1-CitF[2,t+1])/XitF[2,t+1])+ VarLamF[t]*(CitF[2,t+1]/CitF[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitF[i,t]<- (SitF[i,t]^2)*((1-sum(CitF[2:(nage-1),t]))/sum(XitF[2:(nage-1),t])+(1-sum(CitF[3:nage,t+1]))/sum(XitF[3:nage,t+1]))+ VarLamF[t]*(sum(CitF[3:nage,t+1])/sum(CitF[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitF<- replace(Var_SitF,which(is.nan(Var_SitF)==TRUE),0)                     # replace impossible value by 0
Var_SitF<- replace(Var_SitF,which(is.infinite(Var_SitF)==TRUE),0)

SD_SitF <- sqrt(Var_SitF)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitF6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitF6[i]<- sum(XitF[i,])/sum(XitF)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitF7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitF7[i]<- CitF6[i+1]*mean(lamF[1:(nyear-1)])/CitF6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitF7[i]<- sum(CitF6[3:nage])*mean(lamF[1:(nyear-1)])/sum(CitF6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitF7<- replace(SitF7,which(is.nan(SitF7)==TRUE),0)                                 # replace impossible value by 0
SitF7<- replace(SitF7,which(is.infinite(SitF7)==TRUE),0)


Var_SitF7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitF7[i]<- (SitF7[i]^2)*((1-CitF6[i])/XitF[i]+(1-CitF6[i+1])/XitF[i+1])+ mean(VarLamF[1:(nyear-1)])*(CitF6[i+1]/CitF6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitF7[i]<- (SitF7[i]^2)*((1-sum(CitF6[2:(nage-1)]))/sum(XitF[2:(nage-1)])+(1-sum(CitF6[3:nage]))/sum(XitF[3:nage]))+ mean(VarLamF[1:(nyear-1)])*(sum(CitF6[3:nage])/sum(CitF6[2:(nage-1)]))^2                               
      }   # i

Var_SitF7<- replace(Var_SitF7,which(is.nan(Var_SitF7)==TRUE),0)                     # replace impossible value by 0
Var_SitF7<- replace(Var_SitF7,which(is.infinite(Var_SitF7)==TRUE),0)

SD_SitF7 <- sqrt(Var_SitF7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitF8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitF8[i]<- CitF6[i+1]*mean(lamF[1:(nyear-1)])/CitF6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitF8[i]<- sum(CitF6[4:nage])*mean(lamF[1:(nyear-1)])/sum(CitF6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitF8<- replace(SitF8,which(is.nan(SitF8)==TRUE),0)                                 # replace impossible value by 0
SitF8<- replace(SitF8,which(is.infinite(SitF8)==TRUE),0)

Var_SitF8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitF8[i]<- (SitF8[i]^2)*((1-CitF6[i])/XitF[i]+(1-CitF6[i+1])/XitF[i+1])+ mean(VarLamF[1:(nyear-1)])*(CitF6[i+1]/CitF6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitF8[i]<- (SitF8[i]^2)*((1-sum(CitF6[3:(nage-1)]))/sum(XitF[3:(nage-1)])+(1-sum(CitF6[4:nage]))/sum(XitF[4:nage]))+ mean(VarLamF[1:(nyear-1)])*(sum(CitF6[4:nage])/sum(CitF6[3:(nage-1)]))^2                               
      }   # i

Var_SitF8<- replace(Var_SitF8,which(is.nan(Var_SitF8)==TRUE),0)                     # replace impossible value by 0
Var_SitF8<- replace(Var_SitF8,which(is.infinite(Var_SitF8)==TRUE),0)

SD_SitF8 <- sqrt(Var_SitF8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

FouSitT<- rbind( SitF[1:2,],SD_SitF[1:2,])
rownames(FouSitT)<-c("JuvenileS0F","Adults", "SD_S0","SD_SA")

FST<-cbind(FouSitT,c(SitF7[1:2],SD_SitF7[1:2]))
colnames(FST)<-c("2003","2004","2005","2006","Average")

plot(c(2003,2004,2005,2006),FST[1,1:4],type="b",pch=16,ylim=c(0,1.2),xaxt="n",xlim=c(2003,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2003.1,2004.1,2005.1,2006.1),FST[2,1:4],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2003,2004,2005,2006),FST[1,1:4]-FST[3,1:4],c(2003,2004,2005,2006), FST[1,1:4]+FST[3,1:4], angle=90, col="black",code=3,length=0)
  arrows(c(2003.1,2004.1,2005.1,2006.1),FST[2,1:4]-FST[4,1:4],c(2003.1,2004.1,2005.1,2006.1), FST[2,1:4]+FST[4,1:4], angle=90, col="black",code=3,length=0)
axis(1, at = 2003, lab = expression("2003"))
axis(1, at = 2004, lab = expression("2004")) 
axis(1, at = 2005, lab = expression("2005")) 
axis(1, at = 2006, lab = expression("2006")) 

  legend(x =2003, y = 1.25, legend = c("Juveniles : 0.28 +/- 0.06", "Adults : 0.47 +/- 0.11"), pch=16:18,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1.1,"Fougères",font=2)

# Plot difference between juvenile survival

plot(c(2003,2004,2005,2006),FS[1,1:4],type="b",pch=16,ylim=c(0,1),xaxt="n",xlim=c(2003,2006.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2003.1,2004.1,2005.1,2006.1),FST[1,1:4],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2003,2004,2005,2006),FS[1,1:4]-FS[3,1:4],c(2003,2004,2005,2006), FS[1,1:4]+FS[3,1:4], angle=90, col="black",code=3,length=0)
  arrows(c(2003.1,2004.1,2005.1,2006.1),FST[1,1:4]-FST[3,1:4],c(2003.1,2004.1,2005.1,2006.1), FST[1,1:4]+FST[3,1:4], angle=90, col="black",code=3,length=0)
axis(1, at = 2003, lab = expression("2003"))
axis(1, at = 2004, lab = expression("2004")) 
axis(1, at = 2005, lab = expression("2005")) 
axis(1, at = 2006, lab = expression("2006")) 

  legend(x =2003, y = 1.05, legend = c("S0F: 0.24 +/- 0.05","S0T : 0.28 +/- 0.06"), pch=16:17,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2005.8, y=1,"Fougères",font=2)
  
      # three age class constant model
   
TCMF<- cbind(SitF8[1:3],SD_SitF8[1:3])       
TCMF
colnames(TCMF)<-c("Phi","SD_Phi")
rownames(TCMF)<-c("SJ","SY","SA")

plot(c(1,2,3),TCMF[,1],pch=16,cex=1.5,ylim=c(0,1),xaxt="n",xlim=c(0.5,3.5),ylab="Survival probalility with SD",xlab="")
    arrows(c(1,2,3),TCMF[,1]-TCMF[,2],c(1,2,3),TCMF[,1]+TCMF[,2], angle=90, col="black",code=3,length=0)
axis(1, at = 1, lab = expression("Juvenile"))
axis(1, at = 2, lab = expression("Yearling")) 
axis(1, at = 3, lab = expression("Adult")) 
text (x=2, y=1,"Fougères",font=2)  

  
### ADD SURVIVAL RATE TO THE GLOBAL TABLE

DEND <-cbind(DS[1,],t(DST))
DENF <-cbind(FS[1,],t(FST))
DENV <-cbind(VS[1,],t(VST))

DEN <-rbind(DEND,DENF,DENV)
colnames(DEN)<-c("SurYear","SJF","SJT","SA","SD_SJ","SD_SA")          

TotalVR<-cbind(DENS,DEN)                                                 