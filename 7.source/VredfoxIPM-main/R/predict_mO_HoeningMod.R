#' Predict parameters for prior distribution of adult natural mortality using 
#' the Hoenig model
#'
#' Tom Porteus et al. developed an approach to generate informative priors for
#' natural mortality based on maximum (observed) age. It involves a predictive
#' model called the "Hoenig model", the parameters of which were estimated for
#' a range of families using a phylogenetic meta-analysis. 
#' We received the posterior samples from Tom Porteus' published model run upon
#' request.
#' 
#' Reference:
#' Porteus, T. A., Reynolds, J. C., & McAllister, M. K. (2018). 
#' Establishing Bayesian priors for natural mortality rate in carnivore populations. 
#' The Journal of Wildlife Management, 82(8), 1645-1657.
#' 
#' @param datafile character string. Path/file name for file containing posterior
#' samples from Tom Porteus' model run. 
#' @param mu.t.max numeric. Offset parameter for maximum age (= 22.61062). 
#' Provided by Tom Porteus. 
#' @param maxAge integer. Maximum recorded age of harvested animals. 
#' @param nsim integer. Number of simulation replicates for each posterior sample. 
#' @param plot logical. If TRUE, plots predicted prior distributions for natural
#' mortality using the exact posterior samples for Hoenig model parameters vs. 
#' de novo simulation from extracted log-mean and log-sd. 
#'
#' @return a list containing the log mean and log sd of a prior distribution
#' for adult natural mortality (canid familiy).
#' @export
#'
#' @examples

predict_mO_HoenigMod <- function(datafile, mu.t.max, maxAge, nsim, plot){
  
  
  ## Load posterior samples from Hoenig model
  postSamples <- read.table(datafile, header = T)
  
  ## Extract posterior samples for relevant parameters
  param.a <- postSamples[,'a_1']
  param.b <- postSamples[,'b_1']
  sigma.z <- postSamples[,'sigma.z_1']

  ## Simulate prior from maximum recorded age
  # Set the number of posterior samples
  nsamples <- length(param.a)

  # Prepare a data frame to store results
  data <- data.frame()
  
  # Loop over posterior samples and simulate mO
  for(i in 1:nsamples){
    
    log_mu.mO <- param.a[i] + param.b[i] * (log(maxAge) - log(mu.t.max))
    
    mO <- rlnorm(nsim, meanlog = log_mu.mO, sdlog = sigma.z[i])
    data <- rbind(data, data.frame(sampleNo = rep(i, nsim), 
                                   simNo = 1:nsim, 
                                   mO_est = mO))
  }
  
  
  ## Extract log-mean and log-sd for prior distribution
  log_mean <- mean(log(data$mO_est)) 
  log_sd <- sd(log(data$mO_est)) 
  
  ## Optional: Visualize simulated distributions
  if(plot){
    mO_sim <- rlnorm(nrow(data), meanlog = log_mean, sdlog = log_sd)
    
    mO_Comp <- data.frame(
      mO_est = c(data$mO_est, mO_sim),
      Source = rep(c('Posterior estimation', 'Simulation'), each = nrow(data))
    )
    
    ggplot(mO_Comp, aes(x = mO_est)) + 
      geom_density(aes(color = Source, fill = Source), alpha = 0.25) + 
      xlab('Natural mortality prior') + 
      xlim(0, 1) + 
      scale_color_manual(values = c("grey40", "#00C7A9")) +
      scale_fill_manual(values = c("grey40", "#00C7A9")) + 
      theme_bw() + theme(panel.grid = element_blank(), legend.position = 'top')
  }

  
  ## Save and return results
  mO.prior <- list(mnat.logmean = log_mean, mnat.logsd = log_sd)
  
  saveRDS(mO.prior, file = "mO_prior_Parameters.rds")
  
  return(mO.prior)
  
}
