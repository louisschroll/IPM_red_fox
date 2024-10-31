#' Plots transient sensitivites and elasticities
#'
#' @param sensitivities a list of lists containing posterior samples for transient 
#' sensitivities and elasticities for all vital rate parameters as well as
#' population structure (n) and population sizes per age class (N). 
#' @param Amax integer. Number of age classes. 
#'
#' @return a character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotSensitivities <- function(sensitivities, Amax){
  
  ## Create plot folder if not present
  if(!dir.exists("Plots")){dir.create("Plots")}
  
  for(i in 1:2){
    
    ## Select relevant part of data
    params <- sensitivities[[i]]$samples
    
    ## Drop pre-fix for generalising
    names(params) <- stringr::str_split_fixed(names(params), pattern = "_", n = 2)[,2]
    
    ## Extract number of samples
    nosamples <- dim(params[[1]])[1]
    
    #--------------------------#
    # Assemble summarised data #
    #--------------------------#
    
    ## Assemble data (summarized survival/mortality)
    sum.data <- data.frame(
      type = rep(c("Annual survival", "Harvest mortality", "Natural mortality",
                   "Pregnancy rate", "Fetus number",
                   "Denning survival", "Denning mortality",
                   "Immigration rate", "Population structure"), each = nosamples),
      estimate = c(rowSums(params$S) + rowSums(params$Ss), 
                   rowSums(params$mH) + rowSums(params$mHs), rowSums(params$mO),
                   rowSums(params$Psi), rowSums(params$rho),
                   params$S0, params$m0,
                   params$immR, rowSums(params$n))
      
    )
    
    ## Re-order factor levels
    sum.data$type <- factor(sum.data$type, levels = c("Annual survival", "Harvest mortality", "Natural mortality",
                                                      "Pregnancy rate", "Fetus number",
                                                      "Denning survival", "Denning mortality",
                                                      "Immigration rate", "Population structure"))
    
    
    #----------------------------#
    # Assemble age-specific data #
    #----------------------------#
    
    ## Bind all data into a data frame
    age.data <- data.frame(rlist::list.cbind(params))
    
    ## Change column names
    colnames(age.data) <- c(
      paste0("S_", 1:Amax), paste0("mH_", 1:Amax), paste0("mO_", 1:Amax),
      paste0("Psi_", 1:Amax), paste0("rho_", 1:Amax),
      "S0", "m0", 
      paste0("Ss_", 1:Amax), paste0("mHs_", 1:Amax), 
      "immR", paste0("n_", 1:Amax), paste0("N_", 1:Amax)
    )
    
    ## Convert to longitudinal format
    age.data <- reshape2::melt(age.data)
    
    ## Prepare and add information on season
    collapseVars_S <- c(paste0("S_", 1:Amax), paste0("Ss_", 1:Amax))
    collapseVars_mH <- c(paste0("mH_", 1:Amax), paste0("mHs_", 1:Amax))
    summerVars <- c(paste0("Ss_", 1:Amax), paste0("mHs_", 1:Amax))
    winterVars <- c(paste0("S_", 1:Amax), paste0("mH_", 1:Amax))
    
    seasonInfo <- data.frame(variable = c(paste0("S_", 1:Amax), paste0("Ss_", 1:Amax), paste0("mH_", 1:Amax), paste0("mHs_", 1:Amax)),
                             variable2 = c(rep(paste0("S_", 1:Amax), 2), rep(paste0("mH_", 1:Amax), 2)),
                             Season = rep(rep(c("Oct-Jun", "Jul-Sep"), each = Amax), 2))
               
    age.data <- age.data %>%
      dplyr::left_join(seasonInfo, by = "variable")
    
    ## Remove unnecessary data
    age.data <- age.data %>%
      dplyr::filter(!(variable %in% c("Psi_1", "rho_1", paste0("N_", 1:Amax)))) %>%
      dplyr::rename(Parameter = variable,
                    Parameter2 = variable2,
                    Estimate = value)
    
    
    #---------------------------------#
    # Plot sensitivities/elasticities #
    #---------------------------------#
    
    ## Plot colors
    temp.colors <- paletteer::paletteer_c("grDevices::Temps", length(unique(sum.data$type))-3)
    plot.colors <- c("#047993FF", "#005F94FF", temp.colors[1:3], rep(temp.colors[4], 2), temp.colors[5:6])

    ## Summed estimates for all parameters
    addline_format <- function(x,...){
      gsub('\\s','\n',x)
    }
    
    p.sum <- ggplot(sum.data, aes(x = type, y = estimate, group = type)) + 
      geom_violin(aes(fill = type, color = type), alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_fill_manual(values = plot.colors) + 
      scale_color_manual(values = plot.colors) + 
      scale_x_discrete(labels = addline_format(c("Annual survival", "Harvest mortality", "Natural mortality",
                                                 "Pregnancy rate", "Fetus number",
                                                 "Denning survival", "Denning mortality",
                                                 "Immigration rate", "Population structure"))) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))

    
    ## Survival panel
    p.S <- ggplot(subset(age.data, Parameter2 %in% paste0("S_", 1:Amax))) + 
      geom_violin(aes(x = Parameter2, y = Estimate, fill = Season), color = plot.colors[1], alpha = 0.5, scale = 'width', draw_quantiles = 0.5, position = "dodge") + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(S[1], S[2], S[3], S[4], S[5])) + 
      scale_fill_manual(values = c("white", plot.colors[1])) +
      theme_bw() + 
      theme(legend.position = c(0.8, 0.8), 
            legend.key.size = unit(0.5, 'cm'),
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'), 
            legend.title = element_blank(), 
            legend.text = element_text(size = 8),
            panel.grid = element_blank(), 
            axis.text.x = element_text(size = 12), 
            axis.title = element_text(size = 12))
    
    ## Harvest mortality panel
    p.mH <- ggplot(subset(age.data, Parameter2 %in% paste0("mH_", 1:Amax))) + 
      geom_violin(aes(x = Parameter2, y = Estimate, fill = Season), color = plot.colors[2], alpha = 0.5, scale = 'width', draw_quantiles = 0.5, position = "dodge") + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H)) + 
      scale_fill_manual(values = c("white", plot.colors[2])) +
      theme_bw() + 
      theme(legend.position = c(0.8, 0.2), 
            legend.key.size = unit(0.5, 'cm'),
            legend.key.height = unit(0.5, 'cm'),
            legend.key.width = unit(0.5, 'cm'), 
            legend.title = element_blank(), 
            legend.text = element_text(size = 8),
            panel.grid = element_blank(), 
            axis.text.x = element_text(size = 12), 
            axis.title = element_text(size = 12))

    ## Natural mortality panel
    p.mO <- ggplot(subset(age.data, Parameter %in% paste0("mO_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[3], color = plot.colors[3], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    
    ## Pregnancy rate panel
    p.Psi <- ggplot(subset(age.data, Parameter %in% paste0("Psi_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[4], color = plot.colors[4], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(Psi[2], Psi[3], Psi[4], Psi[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    
    ## Fetus number panel
    p.rho <- ggplot(subset(age.data, Parameter %in% paste0("rho_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[5], color = plot.colors[5], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(rho[2], rho[3], rho[4], rho[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))
    
    
    ## Population structure panel
    p.n <- ggplot(subset(age.data, Parameter %in% paste0("n_", 1:Amax)), aes(x = Parameter, y = Estimate, group = Parameter)) + 
      geom_violin(fill = plot.colors[length(plot.colors)], color = plot.colors[length(plot.colors)], alpha = 0.5, scale = 'width', draw_quantiles = 0.5) + 
      ylab(ifelse(i == 1, "Sensitivity", "Elasticity")) + 
      xlab('') + 
      scale_x_discrete(labels = expression(n[1], n[2], n[3], n[4], n[5])) + 
      theme_bw() + 
      theme(legend.position = 'none', panel.grid = element_blank(), axis.text.x = element_text(size = 12), axis.title = element_text(size = 12))

    ## Combine panels and save to pdf
    pdf(paste0("Plots/RedFoxIPM_", ifelse(i == 1, "Sensitivities", "Elasticities"), "_sum.pdf"), width = 10, height = 6)
    print(
      p.sum
    )
    dev.off()
    
    pdf(paste0("Plots/RedFoxIPM_", ifelse(i == 1, "Sensitivities", "Elasticities"), "_age.pdf"), width = 7, height = 8)
    print(
      (p.mH + labs(tag = 'a)') | p.mO + labs(tag = 'b)')) / (p.S + labs(tag = 'c)')| p.Psi + labs(tag = 'd)')) / (p.rho  + labs(tag = 'e)')| p.n + labs(tag = 'f)'))
    )
    dev.off()
    
  }  
  
  ## Return list of plots
  plotList <- c(paste0("Plots/RedFoxIPM_", c("Sensitivities", "Elasticities"), "_sum.pdf"),
                paste0("Plots/RedFoxIPM_", c("Sensitivities", "Elasticities"), "_age.pdf"))
  
  return(plotList)
  
}
