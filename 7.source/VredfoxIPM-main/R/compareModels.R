#' Compare a number of different models
#'
#' @param Amax integer. Number of age classes to consider in analyses.
#' @param Tmax integer. The number of years to consider in analyses.
#' @param minYear integer. First year to consider in analyses.
#' @param maxYear integer. Last year to display in the time series plots. If not 
#' provided, defaults to minYear+Tmax-1.
#' @param logN logical. If TRUE, plots population size at the log scale. If 
#' FALSE (default), plots population size on the natural scale. 
#' @param post.filepaths character vectors containing paths to .rds files 
#' containing posterior samples from models to compare. Can be provided instead
#' of post.list. 
#' @param post.list list containing posterior samples from models to compare
#' in mcmc.list format. Can be procided instead of post.filepaths. 
#' @param model.names character vector with user-defined names for models to 
#' compare. 
#' @param censusCollapse logical vector. Determines for each model whether June 
#' or October should be used as the population census in comparisons. This was
#' necessary for comparing earlier model versions without summer harvest (and hence
#' no October census) to models with summer harvest. In the majority of cases,
#' this will not be relevant any longer and can be ignored (if not provided, it
#' will automatically default to plot the June census population size).
#' @param plotFolder character string containing the path to the folder in which
#' to store comparison plots.  
#' @param returnSumData logical. If TRUE, returns a data frame containing 
#' posterior samples for all compared model as an object into the R workspace. 
#' If FALSE (default), no data is returned.
#'
#' @return a dataframe of posterior samples from all compared models, provided
#' thata returnSumData is set to TRUE. 
#' @export
#'
#' @examples
#' 
compareModels <- function(Amax, Tmax, minYear, maxYear, logN = FALSE, 
                          post.filepaths, post.list, 
                          model.names, censusCollapse, 
                          plotFolder, returnSumData = FALSE){
  
  ## Check models are specified correctly
  if((missing(post.filepaths) & missing(post.list)) |
     (missing(post.filepaths) & missing(post.list))){
    stop("Models have to be specified either via file paths (post.filepaths) or 
         using object names (post.objects).")
  }
  
  ## Make plotting directory if it does not exist already
  if(!dir.exists(plotFolder)){
    dir.create(plotFolder)
  }
  
  ## Count number of models
  nModels <- length(model.names)
  
  ## Make censusCollapse object if not provided
  if(missing(censusCollapse)){
    censusCollapse <- rep(TRUE, nModels)
  }
  
  ## Set maxYear if not provided
  if(missing(maxYear)){
    maxYear <- minYear + Tmax - 1
  }
  
  ## Reformat posterior samples
  post.data <- data.frame()
  for(i in 1:nModels){
    
    # Extract samples for relevant model
    if(!missing(post.list)){
      samples <- as.matrix(post.list[[i]])
    }else{
      samples <- as.matrix(readRDS(post.filepaths[i]))
    }
    
    # Collapse censuses if necessary
    # if(censusCollapse[i]){
    #   
    #   samples[, "N.tot[1]"] <- rowSums(samples[, paste0("octN[", 1:Amax, ", 1]")])
    #   
    #   for(t in 2:(Tmax+1)){
    #     samples[, paste0("N.tot[", t, "]")] <- samples[, paste0("N.tot[", t, "]")] + samples[, paste0("Imm[", t, "]")]
    #   }
    # }
    
    if(!censusCollapse[i]){
      
      samples[, "N.tot[1]"] <- 0
      
      for(t in 2:(Tmax+1)){
        samples[, paste0("N.tot[", t, "]")] <- samples[, paste0("N.tot[", t, "]")] - samples[, paste0("Imm[", t, "]")]
      }
    }
    
    # Change format and add to list
    model.data <- reshape2::melt(samples)
    colnames(model.data) <- c("Sample", "Parameter", "Value")
    model.data$Model <- model.names[i]
    post.data <- rbind(post.data, model.data)
  }
  
  ## Summarize posterior samples into median + 95% CI
  sum.data <- post.data %>%
    
    dplyr::group_by(Parameter, Model) %>%
    
    dplyr::summarise(median = median(Value, na.rm = TRUE),
                     lCI = quantile(Value, probs = 0.025, na.rm = TRUE),
                     uCI = quantile(Value, probs = 0.975, na.rm = TRUE),
                     .groups = "keep") %>%
    
    dplyr::mutate(Parameter = as.character(Parameter)) %>%
    
    dplyr::ungroup()
  
  ## Extract and add age and year information
  idx.data <- data.frame(cbind(unique(sum.data$Parameter), stringr::str_extract_all(unique(sum.data$Parameter), pattern = "\\d+", simplify = TRUE))) %>%
    dplyr::rename("Parameter" = "X1",
                  "Idx1" = "X2", 
                  "Idx2" = "X3") %>%
    
    dplyr::mutate(Idx1 = as.numeric(ifelse(Idx1 %in% c("", 0), NA, Idx1)),
                  Idx2 = as.numeric(ifelse(Idx2 %in% c("", 0), NA, Idx2)),
                  YearIdx = dplyr::case_when(!is.na(Idx2) ~ Idx2,
                                             !is.na(Idx1) & !grepl("Mu", Parameter) ~ Idx1),
                  AgeIdx = dplyr::case_when(!is.na(Idx2) ~ Idx1, 
                                            !is.na(Idx1) & grepl("Mu", Parameter) ~ Idx1),
                  Year = YearIdx + minYear - 1,
                  Age = ifelse(AgeIdx == Amax, paste0(Amax-1, "+"), AgeIdx-1),
                  ParamName = stringr::word(Parameter, 1, sep = "\\["))

    
  sum.data <- sum.data %>%
    dplyr::left_join(idx.data, by = "Parameter")
    
  ## Set parameter groups for plotting posterior density overlaps
  plot.params <- list(
    VRmeans = c(paste0("Mu.mH[", 1:Amax, "]"),
                paste0("Mu.mHs[", 1:Amax, "]"),
                paste0("Mu.mO[", 1:Amax, "]"), 
                paste0("Mu.Psi[", 2:Amax, "]"), 
                paste0("Mu.rho[", 2:Amax, "]"),  
                "Mu.S0", "Mu.immR", "Mu.Imm"),
    
    VReffects = c("sigma.mH", "sigma.mHs", "sigma.mO", "sigma.Psi", "sigma.rho", 
                  "sigma.immR", "logsigma.Imm",
                  "betaHE.mH", 
                  "betaR.Psi", paste0("betaR.Psi[", 2:3, "]"),
                  "betaR.rho", paste0("betaR.rho[", 2:3, "]"),
                  "betaRd.mO", "betaR.mO", "betaRxRd.mO"),
    
    Imm = paste0("Imm[", 1:Tmax, "]"),
    
    Ntot = paste0("N.tot[", 1:Tmax, "]"),
    
    Btot = paste0("B.tot[", 1:Tmax, "]"),
    
    Rtot = paste0("R.tot[", 1:Tmax, "]")
  )
  
  ## Set parameters plotting time series of posterior summaries
  plotTS.paramsAge <- list(
    ParamNames = c("mO", "S", "mH", "mHs", "Psi", "rho", "immR"),
    ParamLabels = c("Natural mortality", "Survival", 
                    "Harvest mortality", "Summer harvest mortality", "Pregnancy rate", "# fetuses/female")
  )

  plotTS.params <- list(
    ParamNames = c("N.tot", "B.tot", "R.tot", "Imm", "immR"),
    ParamLabels = c("Female population size", "# breeding females", "# female recruits", 
                    "# female immigrants", "Immigration rate")
  )
  
  ## Optional: convert population size estimates to log scale
  if(logN){
  
    ## Entire posterior samples
    popN.params <- c(plot.params$Imm, plot.params$Ntot, plot.params$Btot, plot.params$Rtot)
    
    for(i in 1:length(popN.params)){
      post.data$Value[which(post.data$Parameter == popN.params[i])] <- log(post.data$Value[which(post.data$Parameter == popN.params[i])])
    }
    post.data$Value[which(post.data$Value == -Inf)] <- 0
    
    ## Summaries
    popN.rows <- which(sum.data$ParamName %in% c("Imm", "N", "N.tot", "B", "B.tot", "R", "R.tot"))
    sum.data$median[popN.rows] <- ifelse(sum.data$median[popN.rows] == 0, 0, log(sum.data$median[popN.rows]))
    sum.data$lCI[popN.rows] <- ifelse(sum.data$lCI[popN.rows] == 0, 0, log(sum.data$lCI[popN.rows]))
    sum.data$uCI[popN.rows] <- ifelse(sum.data$uCI[popN.rows] == 0, 0, log(sum.data$uCI[popN.rows]))
  }
  
  
  ## Set plotting colors
  plot.cols <- paletteer::paletteer_c("grDevices::Temps", length(model.names))

  ## Plot posterior overlaps
  pdf(paste0(plotFolder, "/PosteriorDensities.pdf"), width = 9, height = 6)
  for(x in 1:length(plot.params)){
    
    print(
      ggplot(subset(post.data, Parameter %in% plot.params[[x]]), aes(x = Value, color = Model, fill = Model)) + 
        geom_density(alpha = 1/nModels) + 
        facet_wrap(~Parameter, scales = "free") + 
        #scale_fill_viridis_d() + scale_color_viridis_d() + 
        scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
        theme_bw() + theme(panel.grid = element_blank())
    )
  }
  dev.off()
  
  ## Plot posterior summary time series for age-specific parameters
  pdf(paste0(plotFolder, "/PosteriorSummaries_TimeSeriesAge.pdf"), width = 8, height = 8)
  for(x in 1:length(plotTS.paramsAge$ParamNames)){
    
    print(
      ggplot(subset(sum.data, ParamName == plotTS.paramsAge$ParamNames[x] & Year <= maxYear), aes(group = Model)) + 
        geom_line(aes(x = Year, y = median, color = Model)) + 
        geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 1/nModels) + 
        scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        facet_wrap(~ Age, ncol = 1, scales = "free_y") + 
        ggtitle(plotTS.paramsAge$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
    
  }
  dev.off()
  
  ## Plot posterior summary time series for age-specific parameters
  pdf(paste0(plotFolder, "/PosteriorSummaries_TimeSeries.pdf"), width = 8, height = 4)
  for(x in 1:length(plotTS.params$ParamNames)){
    
    print(
      ggplot(subset(sum.data, ParamName == plotTS.params$ParamNames[x] & Year > minYear & Year <= maxYear), aes(group = Model)) + 
        geom_line(aes(x = Year, y = median, color = Model)) + 
        geom_ribbon(aes(x = Year, ymin = lCI, ymax = uCI, fill = Model), alpha = 1/nModels) + 
        scale_fill_manual(values = plot.cols) + scale_color_manual(values = plot.cols) + 
        scale_x_continuous(breaks = c(minYear:maxYear), labels = c(minYear:maxYear)) + 
        ggtitle(plotTS.params$ParamLabels[x]) +  
        theme_bw() + theme(panel.grid.minor = element_blank(), 
                           panel.grid.major.y = element_blank(), 
                           axis.text.x = element_text(angle = 45, vjust = 0.5))
    )
    
  }
  dev.off()
  
  ## Optional: return summary data
  if(returnSumData){
    return(sum.data)
  }
}
