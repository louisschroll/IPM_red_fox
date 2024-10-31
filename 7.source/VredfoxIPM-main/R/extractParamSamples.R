#' Extract posterior samples for vital rates and population-level quantitites
#'
#' @param MCMC.samples MCMClist object. The run integrated population model as
#' returned by nimbleMCMC().
#' @param Amax integer. Number of age classes. 
#' @param Tmax integer. Number of years in analysis. 
#'
#' @return a list of lists containing posterior samples for all vital rates and
#' population-level quantities for t > 2. The sublist "t" contains time-specific 
#' parameters while the sublist "t_mean" contains time-average parameters. The 
#' latter are needed for evaluating transient sensitivities. 
#' 
#' @export
#'
#' @examples

extractParamSamples <- function(MCMC.samples, Amax, Tmax){
  
  ###############
  #### SETUP ####
  ###############
  
  ## Convert MCMC samples to matrix
  out.mat <- as.matrix(MCMC.samples)
  
  
  ## Extract number of samples
  nosamples <- dim(out.mat)[1]
  
  
  ## Prepare arrays to rearrange samples - Vital rates & population sizes
  
  # Time-varying vital rates
  S <- mH <- mO <- array(NA, dim = c(nosamples, Amax, Tmax-1))
  Psi <- rho <- array(NA, dim = c(nosamples, Amax, Tmax))
  S0 <- m0 <- immR <- matrix(NA, nrow = nosamples, ncol = Tmax)
  mHs <- array(NA, c(nosamples, Amax, Tmax))
  
  # Time-varying population sizes and growth rates
  N <- n <- array(NA, dim = c(nosamples, Amax, Tmax))
  N_tot <- matrix(NA, nrow = nosamples, ncol = Tmax)
  lambda <- matrix(NA, nrow = nosamples, ncol = Tmax-1)
  
  B <- L <- R <- array(NA, dim = c(nosamples, Amax, Tmax))
  B_tot <- L_tot <- R_tot <- matrix(NA, nrow = nosamples, ncol = Tmax)
  
  # Time-varying immigrant numbers
  Imm <- matrix(NA, nrow = nosamples, ncol = Tmax)
  
  
  ## Fill samples into vectors and matrices
  for(i in 1:nosamples){
    
    for(t in 1:Tmax){
      
      # Time-varying vital rates   
      for(a in 1:Amax){
        Psi[i, a, t] <- out.mat[i, paste0("Psi[", a, ", ", t, "]")]
        rho[i, a, t] <- out.mat[i, paste0("rho[", a, ", ", t, "]")]
        mHs[i, a, t] <- out.mat[i, paste0("mHs[", a, ", ", t, "]")]
        
        if(t < Tmax){
          mH[i, a, t] <- out.mat[i, paste0("mH[", a, ", ", t, "]")]
          mO[i, a, t] <- out.mat[i, paste0("mO[", a, ", ", t, "]")]
          S[i, a, t] <- out.mat[i, paste0("S[", a, ", ", t, "]")]
        }
      }
      
      #S0[i, t] <- out.mat[i, paste0("Mu.S0[", t, "]")]
      S0[i, t] <- out.mat[i, "Mu.S0"]
      m0[i, t] <- -log(S0[i, t])
      
      immR[i, t] <- out.mat[i, paste0("immR[", t, "]")]
      
      # Time-varying population sizes
      for(a in 1:Amax){
        N[i, a, t] <- out.mat[i, paste0("N[", a, ", ", t, "]")]
        n[i, a, t] <- out.mat[i, paste0("N[", a, ", ", t, "]")]/out.mat[i, paste0("N.tot[", t, "]")]
        
        B[i, a, t] <- out.mat[i, paste0("B[", a, ", ", t, "]")]
        L[i, a, t] <- out.mat[i, paste0("L[", a, ", ", t, "]")]
        R[i, a, t] <- out.mat[i, paste0("R[", a, ", ", t, "]")]
      }
      
      N_tot[i, t] <- out.mat[i, paste0("N.tot[", t, "]")]
      
      B_tot[i, t] <- sum(B[i, , t])
      L_tot[i, t] <- sum(L[i, , t])
      R_tot[i, t] <- sum(R[i, , t])
      
      
      # Time-varying immigrant numbers
      Imm[i, t] <- out.mat[i, paste0("Imm[", t, "]")]
      
      # Population growth rate
      if(t < Tmax){
        lambda[i, t] <- out.mat[i, paste0("N.tot[", t+1, "]")]/out.mat[i, paste0("N.tot[", t, "]")]
      }
    }
  }
  
  
  ## Calculate time-average population sizes
  message("The first year is dropped from summaries as no June population size was estimated.")
  
  N_mean <- apply(N[, , 2:Tmax], c(1, 2), mean)
  N_tot_mean <- rowMeans(N_tot)
  n_mean <- N_mean / N_tot_mean
  
  B_mean <- apply(B[, , 2:Tmax], c(1, 2), mean)
  B_tot_mean <- rowMeans(B_tot[, 2:Tmax])
  L_mean <- apply(L[, , 2:Tmax], c(1, 2), mean)
  L_tot_mean <- rowMeans(L_tot[, 2:Tmax])
  R_mean <- apply(R[, , 2:Tmax], c(1, 2), mean)
  R_tot_mean <- rowMeans(R_tot[, 2:Tmax])
  
  
  ## Calculate average immigrant numbers
  Imm_mean <- rowMeans(Imm[, 2:(Tmax-1)])
  
  
  ## Calculate time-average vital rates
  S_mean <- apply(S[, , 2:(Tmax-1)], c(1, 2), mean)
  mH_mean <- apply(mH[, , 2:(Tmax-1)], c(1, 2), mean)
  mO_mean <- apply(mO[, , 2:(Tmax-1)], c(1, 2), mean)
  
  Psi_mean <- apply(Psi[, , 3:Tmax], c(1, 2), mean)
  rho_mean <- apply(rho[, , 3:Tmax], c(1, 2), mean)
  
  mHs_mean <- apply(mHs[, , 2:(Tmax-1)], c(1, 2), mean)
  
  S0_mean <- rowMeans(S0[, 3:Tmax], na.rm = T)
  m0_mean <- rowMeans(m0[, 3:Tmax], na.rm = T)
  
  immR_mean <- rowMeans(immR[, 2:(Tmax-1)], na.rm = T)
  
  ## Make time-average population growth rate
  lambda_mean <- rowMeans(lambda[, 2:(Tmax-1)], na.rm = T)
  
  
  ## Collect parameter samples in a list
  paramSamples <- list(
    
    t = list(S = S,
             mH = mH,
             mO = mO,
             Psi = Psi,
             rho = rho,
             S0 = S0,
             m0 = m0,
             Ss = exp(-mHs),
             mHs = mHs,
             immR = immR,
             n = n,
             N = N, 
             N_tot = N,
             lambda = lambda,
             B = B,
             B_tot = B,
             L = L,
             L_tot = L,
             R = R,
             R_tot = R),
    
    t_mean = list(S = S_mean,
                  mH = mH_mean,
                  mO = mO_mean,
                  Psi = Psi_mean,
                  rho = rho_mean,
                  S0 = S0_mean,
                  m0 = m0_mean,
                  Ss = exp(-mHs_mean),
                  mHs = mHs_mean,
                  immR = immR_mean,
                  n = n_mean,
                  N = N_mean, 
                  N_tot = N_tot_mean,
                  lambda = lambda_mean,
                  B = B_mean,
                  B_tot = B_tot_mean,
                  L = L_mean,
                  L_tot = L_tot_mean,
                  R = R_mean,
                  R_tot = R_tot_mean)
  )
  
  return(paramSamples)
}
