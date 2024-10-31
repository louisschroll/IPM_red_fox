#' Plot vital rate and population-level estimated from IPM
#'
#' @param MCMC.samples an mcmc.list containing the output of the fitted IPM.  
#' @param nim.data a list containing the input data for the model. 
#' @param Amax integer. Number of age classes. 
#' @param Tmax integer. Number of years in the analysis.
#' @param minYear integer. First year in the analysis. 
#' @param logN logical. If TRUE, plots population-level quantities (total and
#' breeding population sizes, numbers of local recruits and immigrants) on the
#' log- instead of natural scale. This makes sense with simulated scenarios. 
#'
#' @returna character vector of plot names. The plots themselves are saved
#' as pdf's in the subfolder "Plots".
#' @export
#'
#' @examples

plotIPM_basicOutputs <- function(MCMC.samples, nim.data, Amax, Tmax, minYear, logN = FALSE){
  
  #-------------------#
  # Data reformatting #
  #-------------------#
  
  ## Reformat posterior samples as longitudinal dataframe
  results <- as.matrix(MCMC.samples) %>%
    as.data.frame() %>%
    dplyr::mutate(`Mu.S[1]` = exp(-(`Mu.mH[1]` + `Mu.mO[1]`)),
                  `Mu.S[2]` = exp(-(`Mu.mH[2]` + `Mu.mO[2]`)),
                  `Mu.S[3]` = exp(-(`Mu.mH[3]` + `Mu.mO[3]`)),
                  `Mu.S[4]` = exp(-(`Mu.mH[4]` + `Mu.mO[4]`)),
                  `Mu.S[5]` = exp(-(`Mu.mH[5]` + `Mu.mO[5]`)),
                  `Mu.Ss[1]` = exp(-(`Mu.mHs[1]`)),
                  `Mu.Ss[2]` = exp(-(`Mu.mHs[2]`)),
                  `Mu.Ss[3]` = exp(-(`Mu.mHs[3]`)),
                  `Mu.Ss[4]` = exp(-(`Mu.mHs[4]`)),
                  `Mu.Ss[5]` = exp(-(`Mu.mHs[5]`)),
                  SampleID = 1:nrow(as.matrix(MCMC.samples))) %>%
    tidyr::pivot_longer(cols = -SampleID) %>%
    dplyr::rename(Parameter = name,
                  Value = value)
  
  ## Extract and add age and year information
  idx.data <- data.frame(cbind(unique(results$Parameter), stringr::str_extract_all(unique(results$Parameter), pattern = "\\d+", simplify = TRUE))) %>%
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
  
  
  results <- results %>%
    dplyr::left_join(idx.data, by = "Parameter")
  

  ## Make dataframe of posterior summaries
  results.sum <- results %>%
    dplyr::group_by(Parameter, Year, Age, ParamName) %>%
    dplyr::summarise(median = median(Value, na.rm = TRUE),
                     lCI = quantile(Value, probs = 0.025, na.rm = TRUE),
                     uCI = quantile(Value, probs = 0.975, na.rm = TRUE),
                     .groups = "keep") %>%
    dplyr::ungroup()
  
  ## Make dataframe containing population proportions
  Ntot.data <- results %>%
    dplyr::filter(ParamName == "N.tot") %>%
    dplyr::select(SampleID, Value, Year) %>%
    dplyr::rename(Ntot = Value)
  
  Btot.data <- results %>%
    dplyr::filter(ParamName == "B.tot") %>%
    dplyr::select(SampleID, Value, Year) %>%
    dplyr::rename(Btot = Value)
  
  results.pp <- results %>%
    dplyr::filter(ParamName %in% c("B", "N")) %>%
    dplyr::left_join(Ntot.data, by = c("Year", "SampleID")) %>%
    dplyr::left_join(Btot.data, by = c("Year", "SampleID")) %>%
    dplyr::mutate(prop = dplyr::case_when(ParamName == "N" ~ Value / Ntot,
                                          ParamName == "B" ~ Value / Btot)) %>%
    dplyr::group_by(Parameter, Year, Age, ParamName) %>%
    dplyr::summarise(median = median(prop, na.rm = TRUE),
                     lCI = quantile(prop, probs = 0.025, na.rm = TRUE),
                     uCI = quantile(prop, probs = 0.975, na.rm = TRUE),
                     .groups = "keep") %>%
    dplyr::ungroup()
  
  ## Optional: convert population-level quantities to log scale
  if(logN){
    suppressWarnings(
      results.sum <- results.sum %>%
        dplyr::mutate(median = dplyr::case_when(!(ParamName %in% c("N.tot", "B.tot", "R.tot", "Imm")) ~ median,
                                                 median > 0 ~ log(median),
                                                 TRUE ~ 0),
                      lCI = dplyr::case_when(!(ParamName %in% c("N.tot", "B.tot", "R.tot", "Imm")) ~ lCI,
                                              lCI > 0 ~ log(lCI),
                                              TRUE ~ 0),
                      uCI = dplyr::case_when(!(ParamName %in% c("N.tot", "B.tot", "R.tot", "Imm")) ~ uCI,
                                              uCI > 0 ~ log(uCI),
                                              TRUE ~ 0))
    )
  }
  
  ## Define plot colors
  plot.colors.param <- c("#047993FF", "#005F94FF", paletteer::paletteer_c("grDevices::Temps", 6))
  plot.colors.age <- paletteer::paletteer_dynamic("cartography::harmo.pal", Amax)
  
  
  #-------------------------------------------#
  # Visualizations of posterior distributions #
  #-------------------------------------------#
  
  ## Age-specific survival and mortality

  # Survival panel
  S.labs <- c("Annual survival (Oct-Jun)", "Summer survival (Jul-Sep)")
  names(S.labs) <- c("Mu.S", "Mu.Ss")
  
  p.S_age <- results %>% 
    dplyr::filter(ParamName %in% c("Mu.S", "Mu.Ss")) %>%
    ggplot(aes(x = Value, group = Age)) + 
    geom_density(aes(y = after_stat(scaled), color = Age, fill = Age), alpha = 0.5) + 
    scale_color_manual(values = plot.colors.age) + 
    scale_fill_manual(values = plot.colors.age) + 
    ylab("Density") +
    facet_wrap(~ ParamName, ncol = 1, scales = "free_y", labeller = labeller(ParamName = S.labs)) + 
    theme_bw() + theme(panel.grid = element_blank(), axis.title.x = element_blank(), legend.position = "none")
  
  # Mortality panel
  m.labs <- c("Harvest mortality (Oct-Jun)", "Natural mortality (Oct-Jun)", "Harvest mortality (Jul-Sep)")
  names(m.labs) <- c("Mu.mH", "Mu.mO", "Mu.mHs")
  
  p.m_age <- results %>% 
    dplyr::filter(ParamName %in% c("Mu.mHs", "Mu.mH", "Mu.mO") & Value < 2) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("Mu.mO", "Mu.mH", "Mu.mHs"))) %>%
    ggplot(aes(x = Value, group = Age)) + 
    geom_density(aes(y = after_stat(scaled), color = Age, fill = Age), alpha = 0.5) + 
    scale_color_manual(values = plot.colors.age) + 
    scale_fill_manual(values = plot.colors.age) + 
    facet_wrap(~ ParamName, ncol = 1, scales = "free_y", labeller = labeller(ParamName = m.labs)) + 
    theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank())
  
  # Combined panel
  pdf("Plots/ParamAvg_Survival&Mortality.pdf", width = 8, height = 5)
  print(p.S_age | p.m_age) + plot_layout(ncol = 2, widths = c(1, 2)) 
  dev.off()
  
  ## Other vital rates
  o.labs <- c("Pregnancy rate", "Fetus number", "Denning survival", "Immigration rate")
  names(o.labs) <- c("Mu.Psi", "Mu.rho", "Mu.S0", "Mu.immR")
  
  p.o_age <- results %>% 
    dplyr::filter(ParamName %in% c("Mu.Psi", "Mu.rho", "Mu.S0", "Mu.immR")) %>%
    dplyr::filter(Age != 0 | is.na(Age)) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("Mu.Psi", "Mu.rho", "Mu.S0", "Mu.immR"))) %>%
    ggplot(aes(x = Value, group = Age)) + 
    geom_density(aes(color = Age, fill = Age), alpha = 0.5) + 
    scale_color_manual(values = plot.colors.age[2:5]) + 
    scale_fill_manual(values = plot.colors.age[2:5]) + 
    facet_wrap(~ ParamName, ncol = 2, scales = "free", labeller = labeller(ParamName = o.labs)) + 
    theme_bw() + theme(panel.grid = element_blank(), axis.title = element_blank())
  
  pdf("Plots/ParamAvg_Reproduction&Immigration.pdf", width = 6, height = 4)
  print(p.o_age) 
  dev.off()
  
  
  #---------------------------------------#
  # Visualizations of posterior summaries #
  #---------------------------------------#
  
  ## Total population size (+ harvest count)
  
  # Population size panel
  p.N_time <- results.sum %>%
    dplyr::filter(ParamName == "N.tot" & Year < minYear+Tmax & Year > minYear) %>%
    ggplot(aes(x = Year)) + 
    geom_line(aes(y = median), color = "#4D004B") + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI), fill = "#4D004B", alpha = 0.5) + 
    scale_x_continuous(limits = c(minYear, minYear+Tmax), breaks = c(minYear:(minYear+Tmax)), labels = c(minYear:(minYear+Tmax))) + 
    ylab("Population size (females)") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  # Harvest count panel
  p.H_time <- data.frame(Hcount = c(colSums(nim.data$C_w)/nim.data$pData_w, colSums(nim.data$C_s)/nim.data$pData_s), 
                         Year = c(as.numeric(colnames(nim.data$C_w)), as.numeric(colnames(nim.data$C_s))),
                         Period = c(rep("Oct-May", ncol(nim.data$C_w)), rep("Jul-Sep", ncol(nim.data$C_s)))) %>%
    dplyr::mutate(Year_offset = Year + 0.5) %>% 
    ggplot(aes(x = Year_offset, y = Hcount, group = Period)) +
    geom_bar(aes(fill = Period), stat = "identity", position = "dodge") + 
    scale_x_continuous(breaks = c(minYear:(minYear+Tmax)), labels = c(minYear:(minYear+Tmax))) + 
    scale_fill_manual(values = c("#BED68AFF", "#005F94FF")) + 
    ylab("# shot (females)") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                       axis.text.x = element_blank(), axis.title.x = element_blank(),
                       legend.position = c(0.1, 0.7))
  
  
  # Combined panel
  pdf("Plots/ParamTime_PopSize.pdf", width = 9, height = 7)
  print(p.H_time / p.N_time + plot_layout(heights = c(0.4, 1))) 
  dev.off()
  
  
  ## Population structure (entire vs. breeding only)
  
  # Total population panel
  p.n_time <- results.pp %>%
    dplyr::filter(ParamName == "N" & Year < minYear+Tmax) %>%
    ggplot(aes(x = Year, y = median, group = Age)) + 
    geom_area(aes(fill = Age)) + 
    ggtitle("a) Population composition") +
    scale_fill_manual(values = plot.colors.age) + 
    scale_x_continuous(breaks = c(minYear:(minYear+Tmax-1)), labels = c(minYear:(minYear+Tmax-1))) + 
    ylab("Proportion/age class (females)") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                       axis.text.x = element_text(angle = 45, vjust = 0.5),
                       axis.title.x = element_blank())
  
  # Breeding population panel
  p.b_time <- results.pp %>%
    dplyr::filter(ParamName == "B" & Year < minYear+Tmax & Age != "0") %>%
    ggplot(aes(x = Year, y = median, group = Age)) + 
    geom_area(aes(fill = Age)) + 
    ggtitle("b) Breeding population composition") + 
    scale_fill_manual(values = plot.colors.age[2:5]) + 
    scale_x_continuous(breaks = c(minYear:(minYear+Tmax-1)), labels = c(minYear:(minYear+Tmax-1))) + 
    ylab("Proportion/age class (females)") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  # Combined panel
  pdf("Plots/ParamTime_PopStructure.pdf", width = 7, height = 6)
  print(p.n_time / p.b_time) 
  dev.off()
  
  ## Local recruits vs. immigrants
  p.RvsImm_time <- results.sum %>%
    dplyr::filter(ParamName %in% c("R.tot", "Imm") & Year < minYear+Tmax & Year > minYear) %>%
    ggplot(aes(x = Year, group = ParamName)) + 
    geom_line(aes(y = median, color = ParamName)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = ParamName), alpha = 0.5) + 
    scale_color_manual(name = "Origin", labels = c("Immigrants", "Local Recruits"), values = c("#E79069FF", "#5DBE85FF")) +
    scale_fill_manual(name = "Origin", labels = c("Immigrants", "Local Recruits"), values = c("#E79069FF", "#5DBE85FF")) +
    scale_x_continuous(breaks = c(minYear:(minYear+Tmax-1)), labels = c(minYear:(minYear+Tmax-1))) + 
    ylab("Number of females") + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                       axis.text.x = element_text(angle = 45, vjust = 0.5))
  
  pdf("Plots/ParamTime_RecruitsVSImmigrants.pdf", width = 8, height = 4)
  print(p.RvsImm_time) 
  dev.off()
  
  
  ## Vital rates over time
  results.sum.VR <- results.sum %>%
    dplyr::filter(ParamName %in% c("S", "mH", "mO", "Psi", "rho", "mHs", "immR") & Year < minYear+Tmax & Age %in% c(NA, 1)) %>%
    dplyr::mutate(ParamName = factor(ParamName, levels = c("S", "mH", "mO", "Psi", "rho", "mHs", "immR"))) %>%
    dplyr::filter(!(ParamName == "mHs" & Year == minYear))
  
  VR.labs <- c("Annual survival", "Winter harvest mortality", "Natural mortality", "Pregnancy rate (1 year old)", "Fetus number (1 year old)", "Summer harvest mortality", "Immigration rate")
  names(VR.labs) <- c("S", "mH", "mO", "Psi", "rho", "mHs", "immR")
  VR.cols <- c(plot.colors.param[1:5], "#785F94", plot.colors.param[7])
  
  p.VRs_time <- results.sum.VR %>%
    ggplot(aes(x = Year, group = ParamName)) + 
    geom_line(aes(y = median, color = ParamName)) + 
    geom_ribbon(aes(ymin = lCI, ymax = uCI, fill = ParamName), alpha = 0.5) + 
    scale_color_manual(values = VR.cols) +
    scale_fill_manual(values = VR.cols) +
    scale_x_continuous(breaks = c(minYear:(minYear+Tmax-1)), labels = c(minYear:(minYear+Tmax-1))) + 
    ylab("Estimate") + 
    facet_wrap(~ ParamName, ncol = 1, scales = "free_y", labeller = labeller(ParamName = VR.labs)) + 
    theme_bw() + theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(), 
                       axis.text.x = element_text(angle = 45, vjust = 0.5), legend.position = "none")
  
  
  pdf("Plots/ParamTime_VitalRates.pdf", width = 7, height = 10)
  print(p.VRs_time) 
  dev.off()
  
  ## Return list of plots
  plotList <- c("Plots/ParamAvg_Survival&Mortality.pdf",
                "Plots/ParamAvg_Reproduction&Immigration.pdf",
                "Plots/ParamTime_PopSize.pdf",
                "Plots/ParamTime_PopStructure.pdf",
                "Plots/ParamTime_RecruitsVSImmigrants.pdf",
                "Plots/ParamTime_VitalRates.pdf")
  
  return(plotList)
}
