#' Plot decomposition of variance in mO into covariates (rodent, reindeer and their interaction) and random effect
#'
#' @param MCMC.samples an mcmc list containing posterior samples from a model run.
#' @param Tmax integer. Number of years to consider in the analysis
#' @param minYear integer. First year to consider in the analysis
#'
#' @return character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".

plotVariance_comp_mO <- function(MCMC.samples, Tmax, minYear){
  
  sam.mat <- as.matrix(MCMC.samples)

  #For natural mortality there are two covariates and three effects, i.e. the formula is: 
  #mO[a, t] <- exp(log(Mu.mO[a]) + betaRd.mO*Reindeer[t] + betaR.mO*RodentAbundance[t+1] + betaRxRd.mO*Reindeer[t]*RodentAbundance[t+1] + epsilon.mO[t])
  #One can choose any age here because the covariates and random effect are modelled independent of age and age is therefore cancelled out
  #we just choose age class 2 for now
  a <-2
  testData.mO <- data.frame()
  nSamples <- nrow(sam.mat)

  #Within a t loop to match with the right rodent and reindeer abundance
  
  for(t in 1:Tmax){
    
    #Posterior of parameter rodent effect on natural mortality
    betaR.mO <- sam.mat[ , "betaR.mO"]
    RodentAbundance <- sam.mat[ , paste0("RodentAbundance[", t+1, "]")] # t+1 because see formula above or in modelcode
    RodentEffect <- betaR.mO * RodentAbundance                          #see formula above or modelcode
    
    #Posterior of parameter reindeer effect on natural mortality
    betaRd.mO <- sam.mat[ , "betaRd.mO"]
    Reindeer <- sam.mat[ , paste0("Reindeer[", t, "]")]                 #t because see formula above or in modelcode 
    ReindeerEffect <- betaRd.mO * Reindeer                              #see formula above or modelcode
    
    #Posterior of parameter rodent x reindeer interaction effect on natural mortality
    betaRxRd.mO <- sam.mat[ , "betaRxRd.mO"]
    InteractionEffect <- betaRxRd.mO * Reindeer * RodentAbundance   #see formula above or modelcode 
    
    #Random effect
    epsilon.mO <- sam.mat[, paste0("epsilon.mO[", t,"]")]
    
    testData_temp <- data.frame(
      Year = t,
      Component = rep(c("Rodent covariate", "Reindeer covariate", "Rodent x Reindeer interaction" , "Random effect"), each = nSamples),
      Value = c(RodentEffect, ReindeerEffect, InteractionEffect, epsilon.mO)
    )
    testData.mO <- rbind(testData.mO, testData_temp)
  }
  
  #This saves your results in a dataframe from which you can then make the comparison. You can, for example, use density plots (or ridgeplots, ggridges is great). 


  testData.mO$Year <- testData.mO$Year+minYear-1 #paste0(testData.mO$Year+minYear-1,"-", testData.mO$Year+minYear)
  
  plot.cols <- paletteer::paletteer_c("grDevices::Temps", length(unique(testData.mO$Component)))
  
  p.mO.decomp <- ggplot(testData.mO) + 
    geom_density(aes(x = Value, color = Component, fill = Component), alpha = 0.5) + 
    scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
    theme_bw() + theme(panel.grid = element_blank()) + 
    facet_wrap(~Year, scales = "free_y", ncol = 3)+
    #ylim(0, 3)+
    xlim(-3, 3)+
    labs(title= "Decomposition of natural mortality (mO) estimates into covariate and random effect")
  
  #Instead of a graphical calculation, you could of course also do a numerical calculation. For that, you'd want to summarize your data frame across samples, and then compare the median [95% CI]. 
  
  pdf("Plots/RedfoxIPM_Variance_decomp_mO.pdf", width = 8, height = 8)
  print(p.mO.decomp)
  dev.off()


  ## Return list of plots
  plotList <- c("Plots/RedfoxIPM_Variance_decomp_mO.pdf")

  return(plotList)
} 