########  DATA CONSTRUCTION OF LIFE TABLE ANALYSIS

### import data

setwd("C:/Users/etu-devillard/Documents/1. THESE/data/renard/Survival")

data<-read.table("FOXDATA.txt",h=T,dec=",")

names (data)                                                                 

### check data and revealing sampling biais : H1 of survival estimate = random sampling of population structure                                                             
barplot(table(data$Age[data$Site=="FDC10"],data$Mode[data$Site=="FDC10"]),beside=T,xlab="Sampling mode",ylab="Number of fox sampled")
legend(x =5, y = 300, legend = c("Juvenile", "Yearlings", "Adults 2", "3","4","5","6+"),lty = c(1, 1), lwd = c(2, 2), col=c("grey10","grey20","grey30","grey40","grey50","grey60","grey70"), cex = 1)
title("Age structure of the FDC10 data set")
                #Unearthind is highly biased, collision have scarce data. 
                #Juvenile sampling is highly biased depending on the method used.

fox<-subset(data,data$Mode!="D"&data$Mode!="R"&data$Age!="0")       # Work only on hunting and trapping and shooting to avoid biais
fox10<-subset(fox,fox$Site=="FDC10")     # Present work on FDC 10

fox10$Site<-factor(fox10$Site)
fox10$Mode<-factor(fox10$Mode)
fox10$Age<-factor(fox10$Age)

####################################################
### MODEL SELECTION ON OVERALL AGE-AT-HARVEST DATA
####################################################

# Age structure of the data

table(fox10$Age,fox10$Year)
 barplot(table(fox10$Year),xlab="Year",ylab="Frequency")

FOX10<-subset(fox10,fox10$Year!="2005")                                         # Remove 2005 year with few data
FOX10$Year<-factor(FOX10$Year)
FOX10$Mode<-factor(FOX10$Mode)

# Differences of age structure between sampling method
m<-table(FOX10$Age,FOX10$Mode)
table(FOX10$Age,FOX10$Mode)
barplot(table(FOX10$Age,FOX10$Mode),beside=T,xlab="Sampling mode",ylab="Number of fox sampled")
legend(x = 1, y = 250, legend = c("Yearlings", "Adults 2", "3","4","5","6+"),lty = c(1, 1), lwd = c(2, 2), col=c("grey10","grey20","grey30","grey40","grey50","grey60","grey70"), cex = 1)
title("Age structure of the FDC10 data set")
barplot(t(m)/colSums(m),beside=T,xlab="Sampling mode",ylab="Prop of fox sampled",legend.text=c("Hunting","Trapping","Shooting"))                                            
chisq.test(m)    # stable age structure between site in all FDC35
# Sex structure of the data
barplot(table(FOX10$Sex,FOX10$Mode),beside=T,col=c("grey","black","grey50"),xlab="Sampling mode",ylab="Number of fox sampled")     
legend(x =1, y = 250, legend = c("Female", "Male"),lty = c(1, 1), lwd = c(2, 2), col=c("grey","black","grey50"), cex = 1)
title("Sex structure of the FDC10 data set")

# Differences of age structure between site
m<-table(FOX10$Age,FOX10$GIC)
barplot(t(m)/colSums(m),beside=T,xlab="GIC",ylab="Prop of fox sampled",legend.text=c("BAR","SAR"))                                            
m1<-rbind(m[1:4,],colSums(m[5:10,]))
chisq.test(m[,c(1,4)])    # stable age structure between site in all FDC35

# Check the assumption of stable age structure

m<-table(FOX10$Age,FOX10$Year)
barplot(t(m)/colSums(m),beside=T,xlab="Age structure",ylab="Prop of fox sampled",legend.text=c("2006","2007","2008","2009","2010","2011"))                                            
m1<-rbind(m[1:4,],colSums(m[5:10,]))
chisq.test(m1)          # unstable age structure due to 2001 and 2002 where only Domagne data
chisq.test(m[,3:6])    # stable age structure between 2003 and 2007 in all FDC35

### lambda estimates from independant distance sampling data : N entre 2006 et 2011 [ rq: CV = ecart-type / moyenne= racine(variance)/mean

DENS<-read.table("Total.txt",h=T,dec=".")
names (DENS)                                                                 

lamB<- DENS$LAMBDA[DENS$GIC=="BAR"]
VarLamB<-DENS$VARLAM[DENS$GIC=="BAR"]

lamS<- DENS$LAMBDA[DENS$GIC=="SAR"]
VarLamS<- DENS$VARLAM[DENS$GIC=="SAR"]

lambda<-colMeans(rbind(lamB,lamS))
VarLam<-colMeans(rbind(VarLamB,VarLamS))

### survival estimate between 2006 and 20011 :

Tot<-table(FOX10$Age,FOX10$Year)
Xit<-Tot
nyear<- ncol(Xit)
nage<-nrow(Xit)
lambda<-lambda
VarLam<-VarLam


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


nbSit<-50                                                                        # nb of survival parameters for AIC 

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

nbSit1<-30                                                                    # nb of survival parameters for AIC 

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

nbSit2<-20                                                                       # nb of survival parameters for AIC 

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

nbSit3<-15                                                                       # nb of survival parameters for AIC 

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

nbSit4<-10                                                                      # nb of survival parameters for AIC 

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

nbSit5<-5                                                                       # nb of survival parameters for AIC 

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

# plot results

matplot(t(Sit4[1:2,]),type="b",pch=16:17,ylim=c(0,1),col=c("black","black"),lwd=2,cex=2,ylab="Survival probalility with SD",xlab="Year")
  arrows(c(1,2,3,4,5),Sit4[1,]-SD_Sit4[1,],c(1,2,3,4,5), Sit4[1,]+SD_Sit4[1,], angle=90, col="black",code=3,length=0.1)
  arrows(c(1,2,3,4,5),Sit4[2,]-SD_Sit4[2,],c(1,2,3,4,5), Sit4[2,]+SD_Sit4[2,], angle=90, col="black",code=3,length=0.1)
title("Survival Probability in FDC10")
  legend(x = 3.5, y = 1, legend = c("Yearlings", "Adults"), pch=16:17,lty = c(1, 1), lwd = c(2, 2), col=c("black","black"), bty = "n", cex = 1)

points(c(2006,2007,2008,2009,2010),t(Sit5[1,]),type="b",pch=16,ylim=c(0,1),col=c("black"),lwd=2,cex=2,ylab="Survival probalility with SD",xlab="Year")
  arrows(c(2006,2007,2008,2009,2010),Sit5[1,]-SD_Sit5[1,],c(2006,2007,2008,2009,2010), Sit5[1,]+SD_Sit5[1,], angle=90, col="black",code=3,length=0.1)
     legend(x = 2006, y = 1, legend = c("All data", "Trapping", "Shooting"), pch=16:18,lty = c(1, 1), lwd = c(2, 2), col=c("black","grey50","grey50"), bty = "n", cex = 1)
############################################
#### Sampling method influence on estimation
############################################

### survival estimate between 2002 and 2007

Hunt<-subset(FOX10,FOX10$Mode=="T")          
Ch<-table(Hunt$Age,Hunt$Year)
Trap<-subset(FOX10,FOX10$Mode=="P")            
Tr<-table(Trap$Age,Trap$Year)

Xit<-Tr                 ## change it to choice the sampling method !
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


# plot results

points(c(2006,2007,2008,2009,2010),t(Chit5[1,]),type="b",pch=17,ylim=c(0,1),col=c("grey50"),lwd=2,cex=2,ylab="Survival probalility with SD",xlab="Year")
  arrows(c(2006,2007,2008,2009,2010),Chit5[1,]-SD_Chit5[1,],c(2006,2007,2008,2009,2010), Chit5[1,]+SD_Chit5[1,], angle=90, col="grey50",code=3,length=0.1)
 
#############################
#### Site by Site estimation
#############################

names(FOX10)

### Construct life table of standing age structure

bar<-as.matrix(table(FOX10$Age[FOX10$GIC=="B"],FOX10$Year[FOX10$GIC=="B"]))
sar<-as.matrix(table(FOX10$Age[FOX10$GIC=="S"],FOX10$Year[FOX10$GIC=="S"]))


###########################################
###  BAROIS GIC   entre 2006 et 2012
###########################################

bar
lamB
VarLamB

### other data
DensB<- DENS$DENSITY[DENS$GIC=="BAR"]
CVDB<-DENS$DCV[DENS$GIC=="BAR"]/100


LSB<-DENS[17:22,36:39]                           # Litter size estimate
rownames(LSB)<-c("2006","2007","2008","2009","2010","2011")

PBB<-DENS[17:22,30:32]                     # Probability of breeding
rownames(PBB)<-c("2006","2007","2008","2009","2010","2011")


### ESTIMATE X0F, number of newborns per year based on the number of female caugth
#################

barF<-as.matrix(table(FOX10$Age[FOX10$GIC=="B"&FOX10$Sex=="F"],FOX10$Year[FOX10$GIC=="B"&FOX10$Sex=="F"]))      #nb of female caugth

X0B<-rep(NA,6)

for (t in 1:length(X0B)) {
  X0B1<- barF[1,t]*PBB[t,1]*LSB[t,1] 
  X0B2<- barF[2,t]*PBB[t,2]*LSB[t,2] 
  X0B35<- sum(barF[3:5,t])*PBB[t,2]*LSB[t,3] 
  X0B610<-sum(barF[6:10,t])*PBB[t,3]*LSB[t,4] 
  
  X0B[t]<- X0B1+X0B2+X0B35+X0B610
      }   # t
      
X0B<-round(X0B)

### survival estimate between 2002 and 2007

XitB<-rbind(X0B,bar )                                                        # add X0 to life table !
nyear<- ncol(XitB)
nage<-nrow(XitB)

X0B/colSums(XitB)
# table of Cit

CitB<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitB[i,t]<- XitB[i,t]/sum(XitB[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t


### Model1:  Assumes survival rates are equal for ages 2-10

# table of Sit

SitB1<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:2) {         		                                               # loop
      SitB1[i,t]<- CitB[i+1,t+1]*lamB[t]/CitB[i,t]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitB1[i,t]<- sum(CitB[4:nage,t+1])*lamB[t]/sum(CitB[3:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitB1<- replace(SitB1,which(is.nan(SitB1)==TRUE),0)                                 # replace impossible value by 0
SitB1<- replace(SitB1,which(is.infinite(SitB1)==TRUE),0)

nbSitB1<-15                                                                     # nb of survival parameters for AIC 

Var_SitB1<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB1[i,t]<- (SitB1[i,t]^2)*((1-CitB[i,t])/XitB[i,t]+(1-CitB[i+1,t+1])/XitB[i+1,t+1])+ VarLamB[t]*(CitB[i+1,t+1]/CitB[i,t])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB1[i,t]<- (SitB1[i,t]^2)*((1-sum(CitB[3:(nage-1),t]))/sum(XitB[3:(nage-1),t])+(1-sum(CitB[4:nage,t+1]))/sum(XitB[4:nage,t+1]))+ VarLamB[t]*(sum(CitB[4:nage,t+1])/sum(CitB[3:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitB1<- replace(Var_SitB1,which(is.nan(Var_SitB1)==TRUE),0)                     # replace impossible value by 0
Var_SitB1<- replace(Var_SitB1,which(is.infinite(Var_SitB1)==TRUE),0)

SD_SitB1 <- sqrt(Var_SitB1)                                                        # matrix of standard deviation of survival estimate 


### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitB<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitB[1,t]<- CitB[2,t+1]*lamB[t]/CitB[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitB[i,t]<- sum(CitB[3:nage,t+1])*lamB[t]/sum(CitB[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitB<- replace(SitB,which(is.nan(SitB)==TRUE),0)                                 # replace impossible value by 0
SitB<- replace(SitB,which(is.infinite(SitB)==TRUE),0)

nbSitB<-10  
                                                                   # nb of survival parameters for AIC 
Var_SitB<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitB[1,t]<- (SitB[1,t]^2)*((1-CitB[1,t])/XitB[1,t]+(1-CitB[2,t+1])/XitB[2,t+1])+ VarLamB[t]*(CitB[2,t+1]/CitB[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB[i,t]<- (SitB[i,t]^2)*((1-sum(CitB[2:(nage-1),t]))/sum(XitB[2:(nage-1),t])+(1-sum(CitB[3:nage,t+1]))/sum(XitB[3:nage,t+1]))+ VarLamB[t]*(sum(CitB[3:nage,t+1])/sum(CitB[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitB<- replace(Var_SitB,which(is.nan(Var_SitB)==TRUE),0)                     # replace impossible value by 0
Var_SitB<- replace(Var_SitB,which(is.infinite(Var_SitB)==TRUE),0)

SD_SitB <- sqrt(Var_SitB)                                                        # matrix of standard deviation of survival estimate 

  
##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitB6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitB6[i]<- sum(XitB[i,])/sum(XitB)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitB7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitB7[i]<- CitB6[i+1]*mean(lamB[1:(nyear-1)])/CitB6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitB7[i]<- sum(CitB6[3:nage])*mean(lamB[1:(nyear-1)])/sum(CitB6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitB7<- replace(SitB7,which(is.nan(SitB7)==TRUE),0)                                 # replace impossible value by 0
SitB7<- replace(SitB7,which(is.infinite(SitB7)==TRUE),0)

nbSitB7<-2

Var_SitB7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB7[i]<- (SitB7[i]^2)*((1-CitB6[i])/XitB[i]+(1-CitB6[i+1])/XitB[i+1])+ mean(VarLamB[1:(nyear-1)])*(CitB6[i+1]/CitB6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB7[i]<- (SitB7[i]^2)*((1-sum(CitB6[2:(nage-1)]))/sum(XitB[2:(nage-1)])+(1-sum(CitB6[3:nage]))/sum(XitB[3:nage]))+ mean(VarLamB[1:(nyear-1)])*(sum(CitB6[3:nage])/sum(CitB6[2:(nage-1)]))^2                               
      }   # i

Var_SitB7<- replace(Var_SitB7,which(is.nan(Var_SitB7)==TRUE),0)                     # replace impossible value by 0
Var_SitB7<- replace(Var_SitB7,which(is.infinite(Var_SitB7)==TRUE),0)

SD_SitB7 <- sqrt(Var_SitB7)                                                        # matrix of standard deviation of survival estimate 


# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitB8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitB8[i]<- CitB6[i+1]*mean(lamB[1:(nyear-1)])/CitB6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitB8[i]<- sum(CitB6[4:nage])*mean(lamB[1:(nyear-1)])/sum(CitB6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitB8<- replace(SitB8,which(is.nan(SitB8)==TRUE),0)                                 # replace impossible value by 0
SitB8<- replace(SitB8,which(is.infinite(SitB8)==TRUE),0)

nbSitB8<-3

Var_SitB8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB8[i]<- (SitB8[i]^2)*((1-CitB6[i])/XitB[i]+(1-CitB6[i+1])/XitB[i+1])+ mean(VarLamB[1:(nyear-1)])*(CitB6[i+1]/CitB6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB8[i]<- (SitB8[i]^2)*((1-sum(CitB6[3:(nage-1)]))/sum(XitB[3:(nage-1)])+(1-sum(CitB6[4:nage]))/sum(XitB[4:nage]))+ mean(VarLamB[1:(nyear-1)])*(sum(CitB6[4:nage])/sum(CitB6[3:(nage-1)]))^2                               
      }   # i

Var_SitB8<- replace(Var_SitB8,which(is.nan(Var_SitB8)==TRUE),0)                     # replace impossible value by 0
Var_SitB8<- replace(Var_SitB8,which(is.infinite(Var_SitB8)==TRUE),0)

SD_SitB8 <- sqrt(Var_SitB8)                                                        # matrix of standard deviation of survival estimate 


##### plot results

BarSit<- rbind(SitB[1:2,],SD_SitB[1:2,])
rownames(BarSit)<-c("JuvenileS0F","Adults", "SD_S0","SD_SA")

BS<-cbind(BarSit,c(SitB7[1:2],SD_SitB7[1:2]))
colnames(BS)<-c("2006","2007","2008","2009","2010","Average")

### ESTIMATE X0T, number of newborns per year based on the total number of removals and assuming a 1:1 sex ratio
#################

barF<-bar*0.5     #nb of female caugth = 0.5* total number of removals


X0B<-rep(NA,6)

for (t in 1:length(X0B)) {
  X0B1<- barF[1,t]*PBB[t,1]*LSB[t,1] 
  X0B2<- barF[2,t]*PBB[t,2]*LSB[t,2] 
  X0B35<- sum(barF[3:5,t])*PBB[t,2]*LSB[t,3] 
  X0B610<-sum(barF[6:10,t])*PBB[t,3]*LSB[t,4] 
  
  X0B[t]<- X0B1+X0B2+X0B35+X0B610
      }   # t
      
X0B<-round(X0B)

### survival estimate between 2002 and 2007

XitB<-rbind(X0B,bar)                                                        # add X0 to life table !
nyear<- ncol(XitB)
nage<-nrow(XitB)

X0B/colSums(XitB)

CitB<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitB[i,t]<- XitB[i,t]/sum(XitB[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### Model1:  Assumes survival rates are equal for ages 2-10

# table of Sit

SitB1<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:2) {         		                                               # loop
      SitB1[i,t]<- CitB[i+1,t+1]*lamB[t]/CitB[i,t]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitB1[i,t]<- sum(CitB[4:nage,t+1])*lamB[t]/sum(CitB[3:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitB1<- replace(SitB1,which(is.nan(SitB1)==TRUE),0)                                 # replace impossible value by 0
SitB1<- replace(SitB1,which(is.infinite(SitB1)==TRUE),0)


Var_SitB1<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB1[i,t]<- (SitB1[i,t]^2)*((1-CitB[i,t])/XitB[i,t]+(1-CitB[i+1,t+1])/XitB[i+1,t+1])+ VarLamB[t]*(CitB[i+1,t+1]/CitB[i,t])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB1[i,t]<- (SitB1[i,t]^2)*((1-sum(CitB[3:(nage-1),t]))/sum(XitB[3:(nage-1),t])+(1-sum(CitB[4:nage,t+1]))/sum(XitB[4:nage,t+1]))+ VarLamB[t]*(sum(CitB[4:nage,t+1])/sum(CitB[3:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitB1<- replace(Var_SitB1,which(is.nan(Var_SitB1)==TRUE),0)                     # replace impossible value by 0
Var_SitB1<- replace(Var_SitB1,which(is.infinite(Var_SitB1)==TRUE),0)

SD_SitB1 <- sqrt(Var_SitB1)                                                        # matrix of standard deviation of survival estimate 

### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitB<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitB[1,t]<- CitB[2,t+1]*lamB[t]/CitB[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitB[i,t]<- sum(CitB[3:nage,t+1])*lamB[t]/sum(CitB[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitB<- replace(SitB,which(is.nan(SitB)==TRUE),0)                                 # replace impossible value by 0
SitB<- replace(SitB,which(is.infinite(SitB)==TRUE),0)


Var_SitB<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitB[1,t]<- (SitB[1,t]^2)*((1-CitB[1,t])/XitB[1,t]+(1-CitB[2,t+1])/XitB[2,t+1])+ VarLamB[t]*(CitB[2,t+1]/CitB[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB[i,t]<- (SitB[i,t]^2)*((1-sum(CitB[2:(nage-1),t]))/sum(XitB[2:(nage-1),t])+(1-sum(CitB[3:nage,t+1]))/sum(XitB[3:nage,t+1]))+ VarLamB[t]*(sum(CitB[3:nage,t+1])/sum(CitB[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitB<- replace(Var_SitB,which(is.nan(Var_SitB)==TRUE),0)                     # replace impossible value by 0
Var_SitB<- replace(Var_SitB,which(is.infinite(Var_SitB)==TRUE),0)

SD_SitB <- sqrt(Var_SitB)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitB6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitB6[i]<- sum(XitB[i,])/sum(XitB)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitB7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitB7[i]<- CitB6[i+1]*mean(lamB[1:(nyear-1)])/CitB6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitB7[i]<- sum(CitB6[3:nage])*mean(lamB[1:(nyear-1)])/sum(CitB6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitB7<- replace(SitB7,which(is.nan(SitB7)==TRUE),0)                                 # replace impossible value by 0
SitB7<- replace(SitB7,which(is.infinite(SitB7)==TRUE),0)

Var_SitB7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB7[i]<- (SitB7[i]^2)*((1-CitB6[i])/XitB[i]+(1-CitB6[i+1])/XitB[i+1])+ mean(VarLamB[1:(nyear-1)])*(CitB6[i+1]/CitB6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB7[i]<- (SitB7[i]^2)*((1-sum(CitB6[2:(nage-1)]))/sum(XitB[2:(nage-1)])+(1-sum(CitB6[3:nage]))/sum(XitB[3:nage]))+ mean(VarLamB[1:(nyear-1)])*(sum(CitB6[3:nage])/sum(CitB6[2:(nage-1)]))^2                               
      }   # i

Var_SitB7<- replace(Var_SitB7,which(is.nan(Var_SitB7)==TRUE),0)                     # replace impossible value by 0
Var_SitB7<- replace(Var_SitB7,which(is.infinite(Var_SitB7)==TRUE),0)

SD_SitB7 <- sqrt(Var_SitB7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitB8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitB8[i]<- CitB6[i+1]*mean(lamB[1:(nyear-1)])/CitB6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitB8[i]<- sum(CitB6[4:nage])*mean(lamB[1:(nyear-1)])/sum(CitB6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitB8<- replace(SitB8,which(is.nan(SitB8)==TRUE),0)                                 # replace impossible value by 0
SitB8<- replace(SitB8,which(is.infinite(SitB8)==TRUE),0)

Var_SitB8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitB8[i]<- (SitB8[i]^2)*((1-CitB6[i])/XitB[i]+(1-CitB6[i+1])/XitB[i+1])+ mean(VarLamB[1:(nyear-1)])*(CitB6[i+1]/CitB6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitB8[i]<- (SitB8[i]^2)*((1-sum(CitB6[3:(nage-1)]))/sum(XitB[3:(nage-1)])+(1-sum(CitB6[4:nage]))/sum(XitB[4:nage]))+ mean(VarLamB[1:(nyear-1)])*(sum(CitB6[4:nage])/sum(CitB6[3:(nage-1)]))^2                               
      }   # i

Var_SitB8<- replace(Var_SitB8,which(is.nan(Var_SitB8)==TRUE),0)                     # replace impossible value by 0
Var_SitB8<- replace(Var_SitB8,which(is.infinite(Var_SitB8)==TRUE),0)

SD_SitB8 <- sqrt(Var_SitB8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

BarSitT<- rbind(SitB[1:2,],SD_SitB1[1:2,])
rownames(BarSitT)<-c("JuvenileS0T","Adults","SD_S0","SD_SA")

BST<-cbind(BarSitT,c(SitB7[1:2],SD_SitB7[1:2]))
colnames(BST)<-c("2006","2007","2008","2009","2010","Average")

plot(c(2006,2007,2008,2009,2010),BST[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2006,2010.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2006.1,2007.1,2008.1,2009.1,2010.1),BST[2,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2006,2007,2008,2009,2010),BST[1,1:5]-BST[3,1:5],c(2006,2007,2008,2009,2010), BST[1,1:5]+BST[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2006.1,2007.1,2008.1,2009.1,2010.1),BST[2,1:5]-BST[4,1:5],c(2006.1,2007.1,2008.1,2009.1,2010.1), BST[2,1:5]+BST[4,1:5], angle=90, col="black",code=3,length=0)

  
  legend(x =2006, y = 1.05, legend = c("Juveniles : 0.25 +/- 0.04", "Adults : 0.51 +/- 0.09"), pch=16:18,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2010.8, y=1,"Barois",font=2)

   # three age class constant model
   
 TCMB<- cbind(SitB8[1:3],SD_SitB8[1:3])       
TCMB
colnames(TCMB)<-c("Phi","SD_Phi")
rownames(TCMB)<-c("SJ","SY","SA")

plot(c(1,2,3),TCMB[,1],pch=16,cex=1.5,ylim=c(0,1),xaxt="n",xlim=c(0.5,3.5),ylab="Survival probalility with SD",xlab="")
    arrows(c(1,2,3),TCMB[,1]-TCMB[,2],c(1,2,3),TCMB[,1]+TCMB[,2], angle=90, col="black",code=3,length=0)
axis(1, at = 1, lab = expression("Juvenile"))
axis(1, at = 2, lab = expression("Yearling")) 
axis(1, at = 3, lab = expression("Adult")) 
text (x=2, y=1,"Barois",font=2)  
    
  # Plot difference between juvenile survival

plot(c(2006,2007,2008,2009,2010),BS[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2006,2010.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2006.1,2007.1,2008.1,2009.1,2010.1),BST[1,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2006,2007,2008,2009,2010),BS[1,1:5]-BS[3,1:5],c(2006,2007,2008,2009,2010), BS[1,1:5]+BS[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2006.1,2007.1,2008.1,2009.1,2010.1),BST[1,1:5]-BST[3,1:5],c(2006.1,2007.1,2008.1,2009.1,2010.1), BST[1,1:5]+BST[3,1:5], angle=90, col="black",code=3,length=0)

  
  legend(x =2006, y = 1.05, legend = c("S0F: 0.24 +/- 0.04","S0T : 0.25 +/- 0.04"), pch=16:17,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2009.8, y=1,"Barrois",font=2)

###########################################
###  Sarce GIC   entre 2006 et 20011
###########################################

sar
lamS
VarLamS


### other data
DensS<- DENS$DENSITY[DENS$GIC=="SAR"]
CVDS<-DENS$DCV[DENS$GIC=="SAR"]/100


LSS<-DENS[23:28,36:39]                       # Litter size estimate
rownames(LSS)<-c("2006","2007","2008","2009","2010","2011")

PBS<-DENS[23:28,30:32]                       # Probability of breeding
rownames(PBS)<-c("2006","2007","2008","2009","2010","2011")

### ESTIMATE X0F, number of newborns per year based on the number of female caugth
#################
sarF<-as.matrix(table(FOX10$Age[FOX10$GIC=="S"&FOX10$Sex=="F"],FOX10$Year[FOX10$GIC=="S"&FOX10$Sex=="F"]))      #nb of female caugth


X0S<-rep(NA,6)

for (t in 1:length(X0S)) {
  X0S1<- sarF[1,t]*PBS[t,1]*LSS[t,1] 
  X0S2<- sarF[2,t]*PBS[t,2]*LSS[t,2] 
  X0S35<- sum(sarF[3:5,t])*PBS[t,2]*LSS[t,3] 
  X0S610<-sum(sarF[6:10,t])*PBS[t,3]*LSS[t,4] 
  
  X0S[t]<- X0S1+X0S2+X0S35+X0S610
      }   # t
      
X0S<-round(X0S)

### survival estimate between 2002 and 2007

XitS<-rbind(X0S,sar)                                                        # add X0 to life table !
nyear<- ncol(XitS)
nage<-nrow(XitS)

X0S/colSums(XitS)


CitS<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitS[i,t]<- XitS[i,t]/sum(XitS[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### Model1:  Assumes survival rates are equal for ages 2-10

# table of Sit

SitS1<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:2) {         		                                               # loop
      SitS1[i,t]<- CitS[i+1,t+1]*lamS[t]/CitS[i,t]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitS1[i,t]<- sum(CitS[4:nage,t+1])*lamS[t]/sum(CitS[3:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitS1<- replace(SitS1,which(is.nan(SitS1)==TRUE),0)                                 # replace impossible value by 0
SitS1<- replace(SitS1,which(is.infinite(SitS1)==TRUE),0)


Var_SitS1<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS1[i,t]<- (SitS1[i,t]^2)*((1-CitS[i,t])/XitS[i,t]+(1-CitS[i+1,t+1])/XitS[i+1,t+1])+ VarLamS[t]*(CitS[i+1,t+1]/CitS[i,t])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS1[i,t]<- (SitS1[i,t]^2)*((1-sum(CitS[3:(nage-1),t]))/sum(XitS[3:(nage-1),t])+(1-sum(CitS[4:nage,t+1]))/sum(XitS[4:nage,t+1]))+ VarLamS[t]*(sum(CitS[4:nage,t+1])/sum(CitS[3:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitS1<- replace(Var_SitS1,which(is.nan(Var_SitS1)==TRUE),0)                     # replace impossible value by 0
Var_SitS1<- replace(Var_SitS1,which(is.infinite(Var_SitS1)==TRUE),0)

SD_SitS1 <- sqrt(Var_SitS1)                                                        # matrix of standard deviation of survival estimate 

### ModelD:  Assumes survival rates are equal for ages 1-8

# table of Sit

SitS<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitS[1,t]<- CitS[2,t+1]*lamS[t]/CitS[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitS[i,t]<- sum(CitS[3:nage,t+1])*lamS[t]/sum(CitS[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitS<- replace(SitS,which(is.nan(SitS)==TRUE),0)                                 # replace impossible value by 0
SitS<- replace(SitS,which(is.infinite(SitS)==TRUE),0)


Var_SitS<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitS[1,t]<- (SitS[1,t]^2)*((1-CitS[1,t])/XitS[1,t]+(1-CitS[2,t+1])/XitS[2,t+1])+ VarLamS[t]*(CitS[2,t+1]/CitS[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS[i,t]<- (SitS[i,t]^2)*((1-sum(CitS[2:(nage-1),t]))/sum(XitS[2:(nage-1),t])+(1-sum(CitS[3:nage,t+1]))/sum(XitS[3:nage,t+1]))+ VarLamS[t]*(sum(CitS[3:nage,t+1])/sum(CitS[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitS<- replace(Var_SitS,which(is.nan(Var_SitS)==TRUE),0)                     # replace impossible value by 0
Var_SitS<- replace(Var_SitS,which(is.infinite(Var_SitS)==TRUE),0)

SD_SitS <- sqrt(Var_SitS)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitS6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitS6[i]<- sum(XitS[i,])/sum(XitS)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitS7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitS7[i]<- CitS6[i+1]*mean(lamS[1:(nyear-1)])/CitS6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitS7[i]<- sum(CitS6[3:nage])*mean(lamS[1:(nyear-1)])/sum(CitS6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitS7<- replace(SitS7,which(is.nan(SitS7)==TRUE),0)                                 # replace impossible value by 0
SitS7<- replace(SitS7,which(is.infinite(SitS7)==TRUE),0)

Var_SitS7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS7[i]<- (SitS7[i]^2)*((1-CitS6[i])/XitS[i]+(1-CitS6[i+1])/XitS[i+1])+ mean(VarLamS[1:(nyear-1)])*(CitS6[i+1]/CitS6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS7[i]<- (SitS7[i]^2)*((1-sum(CitS6[2:(nage-1)]))/sum(XitS[2:(nage-1)])+(1-sum(CitS6[3:nage]))/sum(XitS[3:nage]))+ mean(VarLamS[1:(nyear-1)])*(sum(CitS6[3:nage])/sum(CitS6[2:(nage-1)]))^2                               
      }   # i

Var_SitS7<- replace(Var_SitS7,which(is.nan(Var_SitS7)==TRUE),0)                     # replace impossible value by 0
Var_SitS7<- replace(Var_SitS7,which(is.infinite(Var_SitS7)==TRUE),0)

SD_SitS7 <- sqrt(Var_SitS7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitS8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitS8[i]<- CitS6[i+1]*mean(lamS[1:(nyear-1)])/CitS6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitS8[i]<- sum(CitS6[4:nage])*mean(lamS[1:(nyear-1)])/sum(CitS6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitS8<- replace(SitS8,which(is.nan(SitS8)==TRUE),0)                                 # replace impossible value by 0
SitS8<- replace(SitS8,which(is.infinite(SitS8)==TRUE),0)

Var_SitS8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS8[i]<- (SitS8[i]^2)*((1-CitS6[i])/XitS[i]+(1-CitS6[i+1])/XitS[i+1])+ mean(VarLamS[1:(nyear-1)])*(CitS6[i+1]/CitS6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS8[i]<- (SitS8[i]^2)*((1-sum(CitS6[3:(nage-1)]))/sum(XitS[3:(nage-1)])+(1-sum(CitS6[4:nage]))/sum(XitS[4:nage]))+ mean(VarLamS[1:(nyear-1)])*(sum(CitS6[4:nage])/sum(CitS6[3:(nage-1)]))^2                               
      }   # i

Var_SitS8<- replace(Var_SitS8,which(is.nan(Var_SitS8)==TRUE),0)                     # replace impossible value by 0
Var_SitS8<- replace(Var_SitS8,which(is.infinite(Var_SitS8)==TRUE),0)

SD_SitS8 <- sqrt(Var_SitS8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

SarSit<- rbind(SitS[1:2,],SD_SitS[1:2,])
rownames(SarSit)<-c("JuvenileS0T","Adults", "SD_S0","SD_S1")

SS<-cbind(SarSit,c(SitS7[1:2],SD_SitS7[1:2]))
colnames(SS)<-c("2006","2007","2008","2009","2010","Average")


### ESTIMATE X0T, number of newborns per year based on the total number of removals and assuming a 1:1 sex ratio
#################

sarF<-sar*0.5     #nb of female caugth = 0.5* total number of removals


X0S<-rep(NA,6)

for (t in 1:length(X0S)) {
  X0S1<- sarF[1,t]*PBS[t,1]*LSS[t,1] 
  X0S2<- sarF[2,t]*PBS[t,2]*LSS[t,2] 
  X0S35<- sum(sarF[3:5,t])*PBS[t,2]*LSS[t,3] 
  X0S610<-sum(sarF[6:10,t])*PBS[t,3]*LSS[t,4] 
  
  X0S[t]<- X0S1+X0S2+X0S35+X0S610
      }   # t
      
X0S<-round(X0S)

### survival estimate between 2002 and 2007

XitS<-rbind(X0S,sar)                                                        # add X0 to life table !
nyear<- ncol(XitS)
nage<-nrow(XitS)

X0S/colSums(XitS)
# table of Cit

CitS<- matrix(0, nage, nyear)

        for (t in 1:nyear) {
      for (i in 1:nage) {         			                                         # loop
      CitS[i,t]<- XitS[i,t]/sum(XitS[,t])                                           # Eq (2) from Udevitz 2012
      }   # i
}		# t

### Model1:  Assumes survival rates are equal for ages 2-10

# table of Sit

SitS1<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      for (i in 1:2) {         		                                               # loop
      SitS1[i,t]<- CitS[i+1,t+1]*lamS[t]/CitS[i,t]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitS1[i,t]<- sum(CitS[4:nage,t+1])*lamS[t]/sum(CitS[3:(nage-1),t])        # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitS1<- replace(SitS1,which(is.nan(SitS1)==TRUE),0)                                 # replace impossible value by 0
SitS1<- replace(SitS1,which(is.infinite(SitS1)==TRUE),0)


Var_SitS1<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS1[i,t]<- (SitS1[i,t]^2)*((1-CitS[i,t])/XitS[i,t]+(1-CitS[i+1,t+1])/XitS[i+1,t+1])+ VarLamS[t]*(CitS[i+1,t+1]/CitS[i,t])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS1[i,t]<- (SitS1[i,t]^2)*((1-sum(CitS[3:(nage-1),t]))/sum(XitS[3:(nage-1),t])+(1-sum(CitS[4:nage,t+1]))/sum(XitS[4:nage,t+1]))+ VarLamS[t]*(sum(CitS[4:nage,t+1])/sum(CitS[3:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitS1<- replace(Var_SitS1,which(is.nan(Var_SitS1)==TRUE),0)                     # replace impossible value by 0
Var_SitS1<- replace(Var_SitS1,which(is.infinite(Var_SitS1)==TRUE),0)

SD_SitS1 <- sqrt(Var_SitS1)                                                        # matrix of standard deviation of survival estimate 

### ModelD:  Assumes survival rates are equal for ages 1-8

# taSle of Sit

SitS<- matrix(0, nage, nyear-1)

for (t in 1:(nyear-1)) {
      SitS[1,t]<- CitS[2,t+1]*lamS[t]/CitS[1,t]                                # Eq (6) from Udevitz 2012
      
      for (i in 2:(nage-1)) {         		                                       # loop
      SitS[i,t]<- sum(CitS[3:nage,t+1])*lamS[t]/sum(CitS[2:(nage-1),t])                 # Eq (6) from Udevitz 2012
      }   # i

}		# t

SitS<- replace(SitS,which(is.nan(SitS)==TRUE),0)                                 # replace impossible value by 0
SitS<- replace(SitS,which(is.infinite(SitS)==TRUE),0)


Var_SitS<-  matrix(0, nage, nyear-1)                                             # matrix of variance of survival estimate   

for (t in 1:(nyear-1)) {                                                         # loop
      Var_SitS[1,t]<- (SitS[1,t]^2)*((1-CitS[1,t])/XitS[1,t]+(1-CitS[2,t+1])/XitS[2,t+1])+ VarLamS[t]*(CitS[2,t+1]/CitS[1,t])^2                               

      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS[i,t]<- (SitS[i,t]^2)*((1-sum(CitS[2:(nage-1),t]))/sum(XitS[2:(nage-1),t])+(1-sum(CitS[3:nage,t+1]))/sum(XitS[3:nage,t+1]))+ VarLamS[t]*(sum(CitS[3:nage,t+1])/sum(CitS[2:(nage-1),t]))^2                               
      }                                                  
}		# t

Var_SitS<- replace(Var_SitS,which(is.nan(Var_SitS)==TRUE),0)                     # replace impossible value by 0
Var_SitS<- replace(Var_SitS,which(is.infinite(Var_SitS)==TRUE),0)

SD_SitS <- sqrt(Var_SitS)                                                        # matrix of standard deviation of survival estimate 

##### TIME INVARIANT MODELS ASSUMING STABLE AGE STRUCTURE: 
# table of Cit

CitS6<- rep(0, nage)

for (i in 1:nage) {         			                                               # loop
      CitS6[i]<- sum(XitS[i,])/sum(XitS)                                            # Eq (12) from Udevitz 2012
      }   # i

# Model 7 Time invariant : assumptions of equal survival between 1-8
# table of Sit

SitS7<- matrix(0, nage)  

      for (i in 1:1) {         		                                               # loop
      SitS7[i]<- CitS6[i+1]*mean(lamS[1:(nyear-1)])/CitS6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 2:(nage-1)) {         		                                       # loop
      SitS7[i]<- sum(CitS6[3:nage])*mean(lamS[1:(nyear-1)])/sum(CitS6[2:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitS7<- replace(SitS7,which(is.nan(SitS7)==TRUE),0)                                 # replace impossible value by 0
SitS7<- replace(SitS7,which(is.infinite(SitS7)==TRUE),0)

Var_SitS7<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:1) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS7[i]<- (SitS7[i]^2)*((1-CitS6[i])/XitS[i]+(1-CitS6[i+1])/XitS[i+1])+ mean(VarLamS[1:(nyear-1)])*(CitS6[i+1]/CitS6[i])^2                               
      }   # i
      for (i in 2:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS7[i]<- (SitS7[i]^2)*((1-sum(CitS6[2:(nage-1)]))/sum(XitS[2:(nage-1)])+(1-sum(CitS6[3:nage]))/sum(XitS[3:nage]))+ mean(VarLamS[1:(nyear-1)])*(sum(CitS6[3:nage])/sum(CitS6[2:(nage-1)]))^2                               
      }   # i

Var_SitS7<- replace(Var_SitS7,which(is.nan(Var_SitS7)==TRUE),0)                     # replace impossible value by 0
Var_SitS7<- replace(Var_SitS7,which(is.infinite(Var_SitS7)==TRUE),0)

SD_SitS7 <- sqrt(Var_SitS7)                                                        # matrix of standard deviation of survival estimate 

# Model 8 Time invariant : assumptions of equal survival between 2-8
# table of Sit

SitS8<- matrix(0, nage)  

      for (i in 1:2) {         		                                               # loop
      SitS8[i]<- CitS6[i+1]*mean(lamS[1:(nyear-1)])/CitS6[i]                              # Eq (6) from Udevitz 2012
      }   # i
      for (i in 3:(nage-1)) {         		                                       # loop
      SitS8[i]<- sum(CitS6[4:nage])*mean(lamS[1:(nyear-1)])/sum(CitS6[3:(nage-1)] )       # Eq (6) from Udevitz 2012
      }   # i
      
SitS8<- replace(SitS8,which(is.nan(SitS8)==TRUE),0)                                 # replace impossible value by 0
SitS8<- replace(SitS8,which(is.infinite(SitS8)==TRUE),0)

Var_SitS8<-  matrix(0, nage)                                                      # matrix of variance of survival estimate   

      for (i in 1:2) {         		                                               # Eq (7) from Udevitz 2012
      Var_SitS8[i]<- (SitS8[i]^2)*((1-CitS6[i])/XitS[i]+(1-CitS6[i+1])/XitS[i+1])+ mean(VarLamS[1:(nyear-1)])*(CitS6[i+1]/CitS6[i])^2                               
      }   # i
      for (i in 3:(nage-1)) {         		                                       # Eq (7) from Udevitz 2012
      Var_SitS8[i]<- (SitS8[i]^2)*((1-sum(CitS6[3:(nage-1)]))/sum(XitS[3:(nage-1)])+(1-sum(CitS6[4:nage]))/sum(XitS[4:nage]))+ mean(VarLamS[1:(nyear-1)])*(sum(CitS6[4:nage])/sum(CitS6[3:(nage-1)]))^2                               
      }   # i

Var_SitS8<- replace(Var_SitS8,which(is.nan(Var_SitS8)==TRUE),0)                     # replace impossible value by 0
Var_SitS8<- replace(Var_SitS8,which(is.infinite(Var_SitS8)==TRUE),0)

SD_SitS8 <- sqrt(Var_SitS8)                                                        # matrix of standard deviation of survival estimate 

##### plot results

SarSitT<- rbind(SitS[1:2,],SD_SitS[1:2,])
rownames(SarSitT)<-c("JuvenileS0T","Adults", "SD_S0","SD_SA")

SST<-cbind(SarSitT,c(SitS7[1:2],SD_SitS7[1:2]))
colnames(SST)<-c("2006","2007","2008","2009","2010","Average")

plot(c(2006,2007,2008,2009,2010),SST[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2006,2010.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2006.1,2007.1,2008.1,2009.1,2010.1),SST[2,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2006,2007,2008,2009,2010),SST[1,1:5]-SST[3,1:5],c(2006,2007,2008,2009,2010), SST[1,1:5]+SST[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2006.1,2007.1,2008.1,2009.1,2010.1),SST[2,1:5]-SST[4,1:5],c(2006.1,2007.1,2008.1,2009.1,2010.1), SST[2,1:5]+SST[4,1:5], angle=90, col="black",code=3,length=0)
  
  legend(x =2006, y = 1.05, legend = c("Juveniles : 0.23 +/- 0.06","Adults : 0.41 +/- 0.11"), pch=16:18,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2010, y=1,"Sarce",font=2)
  
  
   # Plot difference between juvenile survival

plot(c(2006,2007,2008,2009,2010),SS[1,1:5],type="b",pch=16,ylim=c(0,1),xlim=c(2006,2010.2),lwd=1,cex=1.5,ylab="Survival probalility with SD",xlab="Year")
points(c(2006.1,2007.1,2008.1,2009.1,2010.1),SST[1,1:5],type="b",pch=17,ylim=c(0,1),lwd=1,cex=1.5)
  arrows(c(2006,2007,2008,2009,2010),SS[1,1:5]-SS[3,1:5],c(2006,2007,2008,2009,2010), SS[1,1:5]+SS[3,1:5], angle=90, col="black",code=3,length=0)
  arrows(c(2006.1,2007.1,2008.1,2009.1,2010.1),SST[1,1:5]-SST[3,1:5],c(2006.1,2007.1,2008.1,2009.1,2010.1), SST[1,1:5]+SST[3,1:5], angle=90, col="black",code=3,length=0)

  
  legend(x =2006, y = 1.05, legend = c("S0F: 0.24+/- 0.06","S0T : 0.23 +/- 0.06"), pch=16:17,lty = c(1, 1), lwd = c(1, 1), col=c("black","black"), bty = "n", cex = 1)
  text (x=2009.8, y=1,"Sarce",font=2)

      # three age class constant model
   
TCMS<- cbind(SitS8[1:3],SD_SitS8[1:3])       
TCMS
colnames(TCMS)<-c("Phi","SD_Phi")
rownames(TCMS)<-c("SJ","SY","SA")

plot(c(1,2,3),TCMS[,1],pch=16,cex=1.5,ylim=c(0,1),xaxt="n",xlim=c(0.5,3.5),ylab="Survival probalility with SD",xlab="")
    arrows(c(1,2,3),TCMS[,1]-TCMS[,2],c(1,2,3),TCMS[,1]+TCMS[,2], angle=90, col="black",code=3,length=0)
axis(1, at = 1, lab = expression("Juvenile"))
axis(1, at = 2, lab = expression("Yearling")) 
axis(1, at = 3, lab = expression("Adult")) 
text (x=2, y=1,"Sarce",font=2)  
 
 
DENB <-cbind(BS[1,],t(BST))
DENSS <-cbind(SS[1,],t(SST))


DEN <-rbind(DENB,DENS)
colnames(DEN)<-c("SurYear","SJF","SJT","SA","SD_SJ","SD_SA")          

                                                   