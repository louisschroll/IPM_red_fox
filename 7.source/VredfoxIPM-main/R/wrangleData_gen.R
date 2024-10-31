
wrangleData_gen <- function(datapath, minYear, onlyFemales = FALSE, poolYrs.genData, threshold = 0.2){
  
  ## Read in data (output from GeneClass analyses)
  gen.data <- read.delim(datapath)
  
  ## Optional: keep only females
  if(onlyFemales){
    gen.data <- subset(gen.data, v_sex == "female")
  }
  
  ## Determine variable (= analytical approach) to use
  gen.data$pvar <- gen.data$pvarV2
  
  ## Remove NAs and transform year variables
  gen.data <- gen.data %>% 
    dplyr::filter(!is.na(pvar)) %>%
    dplyr::mutate(yearIdx_shot = as.numeric(stringr::str_sub(t_hunting_year, 1, 4)) - minYear + 1,
                  yearIdx_birth = v_year_birth - minYear + 1,
                  yearIdx_birth_pre = ifelse(v_year_birth < minYear, v_year_birth - min(v_year_birth) + 1, NA))
  
  ## Make derived variables (immigrant status given a pre-defined threshold, re-scaled probabilities, likelihood proportion)
  gen.data <- gen.data %>%
    dplyr::mutate(pImm = 1 - pvar,
                  pImm_rescaled = 1 - (pvar/max(pvar)),
                  immStatus = ifelse(pvar < threshold, 1, 0),
                  pImm_LL = llotherV2 / (llvarV2 + llotherV2))
  
  ## Extract and list relevant data
  
  if(poolYrs.genData){
    
    # Time-independent analysis (all data pooled)
    gen.data.out <- list(
      pImm = gen.data$pImm,
      pImm_rescaled = gen.data$pImm_rescaled,
      pImm_LL = gen.data$pImm_LL,
      genObs_Imm = sum(gen.data$immStatus),
      genObs_Res = nrow(gen.data) - sum(gen.data$immStatus),
      Xgen = nrow(gen.data)
    )
    
  }else{
    
    # Time-dependent analysis (data separated into "in study period" vs. "before study period")
    inStudy <- which(gen.data$yearIdx_birth > 0)
    
    gen.data.sum <- gen.data %>%
      dplyr::mutate(inStudy = ifelse(yearIdx_birth > 0, 1, 0)) %>%
      dplyr::group_by(inStudy, yearIdx_birth, yearIdx_birth_pre) %>%
      dplyr::summarise(genObs_Imm = sum(immStatus),
                       genObs_All = n(), .groups = "keep") %>%
      dplyr::mutate(genObs_Res = genObs_All - genObs_Imm)
    
    gen.data.out <- list(
      
      pImm_in = gen.data$pImm[inStudy],
      pImm_rescaled_in = gen.data$pImm_rescaled[inStudy],
      pImm_LL_in = gen.data$pImm_LL[inStudy],
      pImm_yrsB_in = gen.data$yearIdx_birth[inStudy],
      pImm_yrsH_in = gen.data$yearIdx_shot[inStudy],
      Xgen_in = length(gen.data$pvar[inStudy]),
      
      genObs_Imm_in = gen.data.sum$genObs_Imm[which(gen.data.sum$inStudy == 1)],
      genObs_Res_in = gen.data.sum$genObs_Res[which(gen.data.sum$inStudy == 1)],
      genObsT_in = nrow(subset(gen.data.sum, inStudy == 1)),
      
      pImm_pre = gen.data$pImm[-inStudy],
      pImm_rescaled_pre = gen.data$pImm_rescaled[-inStudy],
      pImm_LL_pre = gen.data$pImm_LL[-inStudy],
      pImm_yrsB_pre = gen.data$yearIdx_birth_pre[-inStudy],
      pImm_yrsH_pre = gen.data$yearIdx_shot[-inStudy],
      Xgen_pre = length(gen.data$pvar[-inStudy]),
      
      genObs_Imm_pre = gen.data.sum$genObs_Imm[which(gen.data.sum$inStudy == 0)],
      genObs_Res_pre = gen.data.sum$genObs_Res[which(gen.data.sum$inStudy == 0)],
      genObsT_pre = nrow(subset(gen.data.sum, inStudy == 0))
    )
  }
  
  ## Return data
  return(gen.data.out)
}
