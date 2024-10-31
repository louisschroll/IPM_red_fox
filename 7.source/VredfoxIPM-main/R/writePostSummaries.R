#' Summarise posteriors and write to file
#'
#' @param MCMC.samples an mcmc.list containing the output of the fitted IPM. 
#' @param Amax integer. Number of age classes. 
#' @param minYear integer. First year in the analysis. 
#' 
#' @return character vector of file names. The files themselves are saved
#' in the working directory. 
#' @export 
#'
#' @examples

writePostSummaries <- function(MCMC.samples, Amax, minYear){
  
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
  
  ## Save to .rds and .csv
  saveRDS(results.sum, file = "RedfoxIPM_PostSummaries.rds")
  write.csv(results.sum, file = "RedfoxIPM_PostSummaries.csv", row.names = FALSE)
  
  ## Return filenames
  fileNames <- c("RedfoxIPM_PostSummaries.rds", "RedfoxIPM_PostSummaries.csv")
  return(fileNames)
  
}