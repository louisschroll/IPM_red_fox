#' Plots results from a random design transient LTRE
#'
#' @param LTRE_results a list of lists containing results of a fixeddesign
#' transient LTRE.
#' @param Amax integer. Number of age classes. 
#' @param Tmax integer. Number of years in the analysis.
#' @param minYear integer. First year in the analysis.
#' @param HazardRates logical. If TRUE (default), plots results of an LTRE with 
#' mortality hazard rates, else (FALSE) of an LTRE with survival probabilities. 
#' @param PopStructure logical. If TRUE (default), plots results of an LTRE with 
#' population proportions (n), else (FALSE) of an LTRE with age-specific population numbers (N).
#'
#' @return a character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotLTRE_fixedDesign <- function(LTRE_results, Amax, Tmax, minYear, HazardRates, PopStructure){
  
  ## Create plot folder if not present
  if(!dir.exists("Plots")){dir.create("Plots")}
  
  #-------------#
  # Format data #
  #-------------#
  
  ## Select relevant part of data
  contData <- LTRE_results$contData_summary
  
  
  ## Add years
  contData <- contData %>%
    dplyr::mutate(y1 = t1 + minYear - 1,
                  y2 = t2 + minYear - 1) %>%
    dplyr::mutate(yplot = y1 + 1)
  
  
  ## Split off and format summed data
  contData_sum <- contData %>%
    dplyr::filter(Variable %in% c("Stot_sum", "mHtot_sum", "mO_sum",
                                  "Psi_sum", "rho_sum",
                                  "S0", "m0", "immR", 
                                  "n_sum", "N_sum")) %>%
    dplyr::mutate(type = dplyr::case_when(Variable == "Stot_sum" ~ "Annual survival",
                                          Variable == "mHtot_sum" ~ "Harvest mortality",
                                          Variable == "mO_sum" ~ "Natural mortality",
                                          Variable == "Psi_sum" ~ "Pregnancy rate",
                                          Variable == "rho_sum" ~ "Fetus number",
                                          Variable == "S0" ~ "Denning survival",
                                          Variable == "m0" ~ "Denning mortality",
                                          Variable == "immR" ~ "Immigration rate",
                                          Variable == "n_sum" ~ "Population structure",
                                          Variable == "N_sum" ~ "Population size/structure"))
  
  
  ## Make ordered list of present parameter types
  if(HazardRates){
    typeList <- c("Harvest mortality", "Natural mortality",
                  "Pregnancy rate", "Fetus number",
                  "Denning mortality", "Immigration rate")
  }else{
    typeList <- c("Annual survival", 
                  "Pregnancy rate", "Fetus number",
                  "Denning survival", "Immigration rate")
  }
  
  if(PopStructure){
    typeList <- c(typeList, "Population structure")
  }else{
    typeList <- c(typelist, "Population size/structure")
  }
  
  
  ## Re-order factor levels
  contData_sum$type <- factor(contData_sum$type, levels = typeList)
  
  contData$Variable <- factor(contData$Variable, levels = c(paste0("S_", 1:Amax), paste0("mH_", 1:Amax), paste0("mO_", 1:Amax),
                                                            paste0("Psi_", 1:Amax), paste0("rho_", 1:Amax), 
                                                            "S0", "m0", 
                                                            paste0("Ss_", 1:Amax), paste0("mHs_", 1:Amax), 
                                                            "immR", paste0("n_", 1:Amax), paste0("N_", 1:Amax),
                                                            "S_sum", "Stot_sum", "mH_sum", "mHs_sum", "mHtot_sum", "mO_sum", "Psi_sum", "rho_sum", "n_sum", "N_sum"))
  
  contData_sum$Variable <- factor(contData_sum$Variable, levels = c("Stot_sum", "mHtot_sum", "mO_sum", 
                                                                    "Psi_sum", "rho_sum", "S0", "m0",
                                                                    "immR", "n_sum", "N_sum"))
  
  ## Prepare and add information on season to complete data
  collapseVars_S <- c(paste0("S_", 1:Amax), paste0("Ss_", 1:Amax))
  collapseVars_mH <- c(paste0("mH_", 1:Amax), paste0("mHs_", 1:Amax))
  summerVars <- c(paste0("Ss_", 1:Amax), paste0("mHs_", 1:Amax))
  winterVars <- c(paste0("S_", 1:Amax), paste0("mH_", 1:Amax))
  
  seasonInfo <- data.frame(Variable = c(paste0("S_", 1:Amax), paste0("Ss_", 1:Amax), paste0("mH_", 1:Amax), paste0("mHs_", 1:Amax)),
                           Variable2 = c(rep(paste0("S_", 1:Amax), 2), rep(paste0("mH_", 1:Amax), 2)),
                           Season = rep(rep(c("Oct-Jun", "Jul-Sep"), each = Amax), 2))
  
  contData <- contData %>%
    dplyr::left_join(seasonInfo, by = "Variable")
  
  #---------------------------------------------#
  # Plot contributions - Stacked bar/line plots #
  #---------------------------------------------#
  
  ## Plot colors
  temp.colors <- paletteer::paletteer_c("grDevices::Temps", 6)
  
  if(HazardRates){
    plot.colors <- c("#005F94FF", temp.colors)
  }else{
    plot.colors <- c("#047993FF", temp.colors[2:6])
  }
  
  col.offset <- ifelse(HazardRates, 1, 0)
  
  
  ## Barplots for summed contributions
  p.bar <- ggplot(contData_sum, aes(x = yplot, y = median, group = type)) + 
    geom_bar(aes(fill = type, color = type), stat = 'identity', position = 'stack') + 
    ylab("Contribution") + 
    xlab("Year") +
    scale_fill_manual(values = plot.colors) + 
    scale_color_manual(values = plot.colors)  + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    geom_hline(yintercept = 0, linetype = 'dotted') + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) 
  
  pdf(paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_Bars.pdf"), height = 5, width = 8.3)
  print(p.bar)
  dev.off()
  
  
  ## Scaled barplot for absolute summed contributions
  p.bar.sc <- ggplot(contData_sum, aes(x = yplot, y = abs(median), group = type)) + 
    geom_bar(aes(fill = type, color = type), stat = 'identity', position = 'fill') + 
    ylab("Absolute contribution") + 
    xlab("Year") +
    scale_fill_manual(values = plot.colors) + 
    scale_color_manual(values = plot.colors)  + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) 
  
  pdf(paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_ScaledBars.pdf"), height = 5, width = 8.3)
  print(p.bar.sc)
  dev.off()
  
  
  ## Stacked lineplot for absolute summed contributions
  p.stacklines.sc <- ggplot(contData_sum, aes(x = yplot, y = abs(median), group = type)) + 
    geom_area(aes(fill = type, color = type), stat = 'identity', position = 'fill') + 
    ylab("Absolute contribution") + 
    xlab("Year") +
    scale_fill_manual(values = plot.colors) + 
    scale_color_manual(values = plot.colors)  + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) 
  
  pdf(paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_StackedLines.pdf"), height = 5, width = 8.3)
  print(p.stacklines.sc)
  dev.off()
  
  
  ## Stacked lineplots for age-specific parameters
  
  if(!HazardRates){

    # Survival probability
    summer.color <- "#785F94"
    S.params <- c(paste0("S_", 1:Amax), paste0("Ss_", 1:Amax))
    S.labels <- expression(S[1], S[2], S[3], S[4], S[5],
                           Ss[1], Ss[2], Ss[3], Ss[4], Ss[5])
    S.colors <- c(plot.colors[1],
                  alpha(plot.colors[1], 0.8),
                  alpha(plot.colors[1], 0.6),
                  alpha(plot.colors[1], 0.4),
                  alpha(plot.colors[1], 0.2),
                  summer.color,
                  alpha(summer.color, 0.8),
                  alpha(summer.color, 0.6),
                  alpha(summer.color, 0.4),
                  alpha(summer.color, 0.2))
    
    p.stacklines.S <- ggplot(subset(contData, Variable %in% S.params), aes(x = yplot, y = abs(median), group = Variable)) + 
      geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
      ylab('Absolute contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = S.colors, labels = S.labels) + 
      scale_color_manual(values = S.colors, labels = S.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  }else{
    
    # Harvest mortality
    summer.color <- "#785F94"
    mH.params <- c(paste0("mH_", 1:Amax), paste0("mHs_", 1:Amax))
    mH.labels <- expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H,
                            m[1]^Hs, m[2]^Hs, m[3]^Hs, m[4]^Hs, m[5]^Hs)
    mH.colors <- c(plot.colors[1],
                   alpha(plot.colors[1], 0.8),
                   alpha(plot.colors[1], 0.6),
                   alpha(plot.colors[1], 0.4),
                   alpha(plot.colors[1], 0.2),
                   summer.color,
                   alpha(summer.color, 0.8),
                   alpha(summer.color, 0.6),
                   alpha(summer.color, 0.4),
                   alpha(summer.color, 0.2))
    
    p.stacklines.mH <- contData %>%
      dplyr::filter(Variable %in% mH.params) %>%
      dplyr::mutate(Variable = factor(Variable, levels = mH.params)) %>%
      ggplot(aes(x = yplot, y = abs(median))) + 
      geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
      ylab('Absolute contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = mH.colors, labels = mH.labels) + 
      scale_color_manual(values = mH.colors, labels = mH.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    # Natural mortality
    mO.params <- paste0("mO_", 1:Amax)
    mO.labels <- expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)
    mO.colors <- c(plot.colors[2], alpha(plot.colors[2], 0.8), alpha(plot.colors[2], 0.6), alpha(plot.colors[2], 0.4), alpha(plot.colors[2], 0.2))
    
    p.stacklines.mO <- ggplot(subset(contData, Variable %in% mO.params), aes(x = yplot, y = abs(median), group = Variable)) + 
      geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
      ylab('Absolute contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = mO.colors, labels = mO.labels) + 
      scale_color_manual(values = mO.colors, labels = mO.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  }

  # Pregnancy rate
  Psi.params <- paste0("Psi_", 2:Amax)
  Psi.labels <- expression(Psi[2], Psi[3], Psi[4], Psi[5])
  Psi.colors <- c(alpha(plot.colors[2+col.offset], 0.8), alpha(plot.colors[2+col.offset], 0.6), alpha(plot.colors[2+col.offset], 0.4), alpha(plot.colors[2+col.offset], 0.2))
  
  p.stacklines.Psi <- ggplot(subset(contData, Variable %in% Psi.params), aes(x = yplot, y = abs(median), group = Variable)) + 
    geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
    ylab('Absolute contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = Psi.colors, labels = Psi.labels) + 
    scale_color_manual(values = Psi.colors, labels = Psi.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  # Fetus number
  rho.params <- paste0("rho_", 2:Amax)
  rho.labels <- expression(rho[2], rho[3], rho[4], rho[5])
  rho.colors <- c(alpha(plot.colors[3+col.offset], 0.8), alpha(plot.colors[3+col.offset], 0.6), alpha(plot.colors[3+col.offset], 0.4), alpha(plot.colors[3+col.offset], 0.2))
  
  p.stacklines.rho <- ggplot(subset(contData, Variable %in% rho.params), aes(x = yplot, y = abs(median), group = Variable)) + 
    geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
    ylab('Absolute contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = rho.colors, labels = rho.labels) + 
    scale_color_manual(values = rho.colors, labels = rho.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  # Population structure
  if(PopStructure){
    n.params <- paste0("n_", 1:Amax)
    n.labels <- expression(n[1], n[2], n[3], n[4], n[5])
  }else{
    n.params <- paste0("N_", 1:Amax)
    n.labels <- expression(N[1], N[2], N[3], N[4], N[5])
  }
  n.colors <- c(plot.colors[length(plot.colors)], alpha(plot.colors[length(plot.colors)], 0.8), alpha(plot.colors[length(plot.colors)], 0.6), alpha(plot.colors[length(plot.colors)], 0.4), alpha(plot.colors[length(plot.colors)], 0.2))
  
  p.stacklines.n <- ggplot(subset(contData, Variable %in% n.params), aes(x = yplot, y = abs(median), group = Variable)) + 
    geom_area(aes(fill = Variable), color = NA, stat = 'identity', position = 'fill') + 
    ylab('Absolute contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = n.colors, labels = n.labels) + 
    scale_color_manual(values = n.colors, labels = n.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  # Plot separate panels together
  if(HazardRates){
    pdf("Plots/RedFoxIPM_fixedLTRE_MHR_StackedLines_Groups.pdf", width = 7, height = 10)
    print(p.stacklines.mH / p.stacklines.mO / p.stacklines.Psi / p.stacklines.rho / p.stacklines.n + plot_layout(heights = c(1.75, rep(1, 4))))
    dev.off()
  }else{
    pdf("Plots/RedFoxIPM_fixedLTRE_SP_StackedLines_Groups.pdf", width = 7, height = 10)
    print(p.stacklines.S / p.stacklines.Psi / p.stacklines.rho / p.stacklines.n / plot_spacer()  + plot_layout(heights = c(1.75, rep(1, 4))))
    dev.off()
  }

  
  #------------------------------------------#
  # Plot contributions - Line & ribbon plots #
  #------------------------------------------#
  
  
  ## Line & ribbon plot for absolute summed contributions
  p.ribbon <- ggplot(contData_sum, aes(x = yplot, y = median, group = type)) + 
    geom_line(aes(color = type)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = type), alpha = 0.3) +
    ylab("Contribution") + 
    xlab("Year") +
    scale_fill_manual(values = plot.colors) + 
    scale_color_manual(values = plot.colors)  + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5)) 
  
  pdf(paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_Ribbon.pdf"), height = 5, width = 8.3)
  print(p.ribbon)
  dev.off()
  
  
  ## Line & ribbon plots for age-specific parameters
  
  col.offset <- ifelse(HazardRates, 1, 0)
  
  if(!HazardRates){
    
    # Survival probability
    S.params <- paste0("S_", 1:Amax)
    S.labels <- expression(S[1], S[2], S[3], S[4], S[5])
    S.colors <- plot.colors[1:Amax]
    
    p.ribbon.S <- ggplot(subset(contData, Variable %in% S.params), aes(x = yplot, y = median, group = Variable)) + 
      geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
      ylab('Contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = S.colors, labels = S.labels) + 
      scale_color_manual(values = S.colors, labels = S.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

    # Summer survival probability
    Ss.params <- paste0("Ss_", 1:Amax)
    Ss.labels <- expression(Ss[1], Ss[2], Ss[3], Ss[4], Ss[5])
    Ss.colors <- plot.colors[1:Amax]
    
    p.ribbon.Ss <- ggplot(subset(contData, Variable %in% Ss.params), aes(x = yplot, y = median, group = Variable)) + 
      geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
      ylab('Contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = Ss.colors, labels = Ss.labels) + 
      scale_color_manual(values = Ss.colors, labels = Ss.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
    
  }else{
    
    # Harvest mortality
    mH.params <- paste0("mH_", 1:Amax)
    mH.labels <- expression(m[1]^H, m[2]^H, m[3]^H, m[4]^H, m[5]^H)
    mH.colors <- plot.colors[1:Amax]
    
    p.ribbon.mH <- ggplot(subset(contData, Variable %in% mH.params), aes(x = yplot, y = median, group = Variable)) + 
      geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
      ylab('Contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = mH.colors, labels = mH.labels) + 
      scale_color_manual(values = mH.colors, labels = mH.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

    # Summer harvest mortality
    mHs.params <- paste0("mHs_", 1:Amax)
    mHs.labels <- expression(m[1]^Hs, m[2]^Hs, m[3]^Hs, m[4]^Hs, m[5]^Hs)
    mHs.colors <- plot.colors[1:Amax]
    
    p.ribbon.mHs <- ggplot(subset(contData, Variable %in% mHs.params), aes(x = yplot, y = median, group = Variable)) + 
      geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
      ylab('Contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = mHs.colors, labels = mHs.labels) + 
      scale_color_manual(values = mHs.colors, labels = mHs.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
    
    # Natural mortality
    mO.params <- paste0("mO_", 1:Amax)
    mO.labels <- expression(m[1]^O, m[2]^O, m[3]^O, m[4]^O, m[5]^O)
    mO.colors <- plot.colors[1:Amax]
    
    p.ribbon.mO <- ggplot(subset(contData, Variable %in% mO.params), aes(x = yplot, y = median, group = Variable)) + 
      geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
      ylab('Contribution') + 
      xlab("Year") + 
      scale_fill_manual(values = mO.colors, labels = mO.labels) + 
      scale_color_manual(values = mO.colors, labels = mO.labels) + 
      scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
      theme_bw() + 
      theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  }
  
  # Pregnancy rate
  Psi.params <- paste0("Psi_", 2:Amax)
  Psi.labels <- expression(Psi[2], Psi[3], Psi[4], Psi[5])
  Psi.colors <- plot.colors[2:Amax]
  
  p.ribbon.Psi <- ggplot(subset(contData, Variable %in% Psi.params), aes(x = yplot, y = median, group = Variable)) + 
    geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
    ylab('Contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = Psi.colors, labels = Psi.labels) + 
    scale_color_manual(values = Psi.colors, labels = Psi.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  # Fetus number
  rho.params <- paste0("rho_", 2:Amax)
  rho.labels <- expression(rho[2], rho[3], rho[4], rho[5])
  rho.colors <- plot.colors[2:Amax]
  
  p.ribbon.rho <- ggplot(subset(contData, Variable %in% rho.params), aes(x = yplot, y = median, group = Variable)) + 
    geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
    ylab('Contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = rho.colors, labels = rho.labels) + 
    scale_color_manual(values = rho.colors, labels = rho.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  # Population structure
  if(PopStructure){
    n.params <- paste0("n_", 1:Amax)
    n.labels <- expression(n[1], n[2], n[3], n[4], n[5])
  }else{
    n.params <- paste0("N_", 1:Amax)
    n.labels <- expression(N[1], N[2], N[3], N[4], N[5])
  }
  n.colors <- plot.colors[1:Amax]
  
  p.ribbon.n <- ggplot(subset(contData, Variable %in% n.params), aes(x = yplot, y = median, group = Variable)) + 
    geom_line(aes(color = Variable)) + geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = Variable), alpha = 0.3) +
    ylab('Contribution') + 
    xlab("Year") + 
    scale_fill_manual(values = n.colors, labels = n.labels) + 
    scale_color_manual(values = n.colors, labels = n.labels) + 
    scale_x_continuous(breaks = c(minYear:max(contData_sum$yplot)), labels = c(minYear:max(contData_sum$yplot))) + 
    theme_bw() + 
    theme(legend.title = element_blank(), panel.grid = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))

  
  # Plot separate panels together
  if(HazardRates){
    pdf("Plots/RedFoxIPM_fixedLTRE_MHR_Ribbon_Groups.pdf", width = 7, height = 12)
    print(p.ribbon.mH / p.ribbon.mHs / p.ribbon.mO / p.ribbon.Psi / p.ribbon.rho / p.ribbon.n)
    dev.off()
  }else{
    pdf("Plots/RedFoxIPM_fixedLTRE_SP_Ribbon_Groups.pdf", width = 7, height = 12)
    print(p.ribbon.S / p.ribbon.Ss / p.ribbon.Psi / p.ribbon.rho / p.ribbon.n / plot_spacer())
    dev.off()
  }
  
  ## Return list of plots
  plotList <- c(paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_Bars.pdf"),
                paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_ScaledBars.pdf"),
                paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_StackedLines.pdf"),
                paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_StackedLines_Groups.pdf"),
                paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_Ribbon.pdf"),
                paste0("Plots/RedFoxIPM_fixedLTRE_", ifelse(HazardRates, "MHR", "SP"), "_Ribbon_Groups.pdf"))
  
  return(plotList)
  
}