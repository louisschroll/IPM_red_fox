#' Calculate posterior summaries of parameter correlations with density
#'
#' @param MCMC.samples an mcmc list containing posterior samples from a model run.
#' @param Tmax integer. Number of years to consider in the analysis
#'
#' @return a vector of file names. The files are saved to the root directory. 
#' @export
#'
#' @examples

calculate_p.hoc_param.corr <- function(MCMC.samples, 
                                       Tmax){
  
  ## Convert posterior samples to matrix format
  out.mat <- as.matrix(MCMC.samples)
  
  # Determine year range
  yearIdxs <- (1:Tmax)
  
  # Set up vector for storing correlation coefficients
  corr_lambda <- corr_mH <- corr_mO <- corr_immR <-corr_Psi <- corr_rho <- rep(NA, nrow(out.mat))
  corr_mH.mO <- corr_mH.immR <- corr_mH.immR.delay <- corr_mH.mO.delay <-corr_mO.immR <-  corr_mH.rep <- corr_mO.rep <-  rep(NA, nrow(out.mat))
  
  for(i in 1:nrow(out.mat)){
    
    # Extract parameter time series
    age_class <- 2 #we just choose age class 2 for now
    
    densT <- unname(out.mat[i,  paste0("N.tot[", yearIdxs,"]")])
    
    lambdaT <- c(densT[2:length(densT)]/densT[1:(length(densT)-1)], NA)
    
    mHT <- unname(out.mat[i, paste0("mH[",age_class,", " , yearIdxs, "]")])
    mOT <- unname(out.mat[i, paste0("mO[",age_class,", " , yearIdxs, "]")])
    PsiT <- unname(out.mat[i, paste0("Psi[",age_class,", " , yearIdxs, "]")])
    rhoT <- unname(out.mat[i, paste0("rho[",age_class,", " , yearIdxs, "]")])
    R.totT <- unname(out.mat[i, paste0("R.tot[" , yearIdxs, "]")])
    immRT <- unname(out.mat[i, paste0("immR[" , yearIdxs, "]")])
    
    
    # Calculate correlation coefficients
    
    #-- simple density dependence
    corr_lambda[i] <- cor.test(lambdaT[2:(Tmax-1)], densT[2:(Tmax-1)])$estimate        #here we use indexing to not use the first and last year for which there is no growth rate  
    corr_mH[i] <- cor.test(mHT[2:Tmax], densT[2:Tmax])$estimate                        #the first year pop size is set to 0 so remove the first year with indexing
    corr_mO[i] <- cor.test(mOT[2:Tmax], densT[2:Tmax])$estimate                        #the first year pop size is set to 0 so remove the first year with indexing
    corr_immR[i] <- cor.test(immRT[2:Tmax], densT[2:Tmax])$estimate                    #the first year pop size is set to 0 so remove the first year with indexing
    corr_Psi[i] <- cor.test(PsiT[2:Tmax], (densT[2:Tmax]- R.totT[2:Tmax]))$estimate   #N-R = adult population size, because reproduction should not depend on offspring that are not there yet
    corr_rho[i] <- cor.test(rhoT[2:Tmax], (densT[2:Tmax]- R.totT[2:Tmax]))$estimate   #N-R = adult population size, because reproduction should not depend on offspring that are not there yet)
    
    #-- compensatory mechanisms
    corr_mH.mO[i] <-cor.test(mHT, mOT)$estimate           #this one I would expect to be compensatory
    corr_mH.immR[i] <-cor.test(mHT, immRT)$estimate       #this one I would expect to be compensatory, but I have more faith in the delayed version because immigration occurs afterwards in the next autumn
    
    #But, i think compensation by harvest more likely to occur in the next year, at least for immigration if we say that immigration early in the year and much of harvest late in the year
    #so I added these 2:
    corr_mH.immR.delay[i] <-cor.test(mHT[1:(Tmax-1)], immRT[2:Tmax])$estimate         #I have most faith in this one to be compensatory
    corr_mH.mO.delay[i]   <-cor.test(mHT[1:(Tmax-1)], mOT[2:Tmax])$estimate           #I have equal faith in this one as the non delayed version to be compensatory
    
    corr_mO.immR[i] <-cor.test(mOT, immRT)$estimate                                   #this one I would not expect to be compensatory but additive, because immigration and Natural mortality likely driven by same environment
    corr_mH.rep[i]  <-cor.test(mHT[1:(Tmax-1)], (PsiT[2:Tmax]*rhoT[2:Tmax]))$estimate #note indexing, checks if reproductive output is higher/lower following winter with more mortality
    corr_mO.rep[i]  <-cor.test(mOT[1:(Tmax-1)], (PsiT[2:Tmax]*rhoT[2:Tmax]))$estimate #same as above for nat mort. but this one I would expect to be additive rather than compensatory because dependent on same environment
    
  }
  
  par(mfrow = c(3, 2))
  
  hist(corr_lambda, main = paste("Lambda~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mH, main = paste("Harvest~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mO, main = paste("Nat.mort~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_immR, main = paste("Immigration~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_Psi, main = paste("Breeding~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_rho, main = paste("Litter~N"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  
  
  par(mfrow = c(3, 2))
  
  hist(corr_mH.mO, main = paste("Harvest~Nat.mort"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mH.immR, main = paste("Harvest~Imm"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mH.immR.delay, main = paste("Harvest~Imm[t+1]"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mH.mO.delay, main = paste("Harvest~Nat.mort[t+1]"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mO.immR, main = paste("Nat.mort~Imm"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  hist(corr_mH.rep, main = paste("Harvest~Reproduction"), xlab=NULL, ylab=NULL)
  abline(v=0, col="blue")
  
  # Summarize posteriors in data frame
  corr_data <- rbind(quantile(corr_lambda, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mH, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mO, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_immR, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_Psi, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_rho, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     
                     quantile(corr_mH.mO, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mH.immR, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mH.immR.delay, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mH.mO.delay, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mO.immR, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mH.rep, probs = c(0.5, 0.025, 0.25, 0.75, 0.975)),
                     quantile(corr_mO.rep, probs = c(0.5, 0.025, 0.25, 0.75, 0.975))) %>%
    tibble::as_tibble() %>%
    dplyr::mutate(Parameter = c("Lambda~N", "Harvest~N", "Nat.mort~N", "Immigration~N", "Breeding~N", "Litter~N", 
                                "Harvest~Nat.mort", "Harvest~Imm", "Harvest~Imm[t+1]", "Harvest~Nat.mort[t+1]", "Nat.mort~Imm", "Harvest~Reproduction", "Nat.mort~Reproduction"),
                  N_years = length(yearIdxs),
                  .before = `50%`) %>%
    dplyr::mutate(Evidence = dplyr::case_when(sign(`2.5%`) == sign(`97.5%`) ~ "**",
                                              sign(`25%`) == sign(`75%`) ~ "*",
                                              TRUE ~ "-"))
  
## Save results as RDS and csv
  saveRDS(corr_data, file = "p.hoc_param.corr.rds")
  write.csv(corr_data, "p.hoc_param.corr.csv", row.names = FALSE)


## Return filepaths
filepaths <- c("p.hoc_param.corr.rds", "p.hoc_param.corr.csv")

return(filepaths)
}

