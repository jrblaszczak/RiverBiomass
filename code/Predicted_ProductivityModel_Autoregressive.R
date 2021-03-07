## Growth Model 1 - Data simulation 

PM_AR <- function(phi, alpha, beta, sig_p, df) {
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  obs_err <- GPP_sd
  
  ## Vectors for model output
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- GPP[1]
  l_pred_GPP <- numeric(Ndays)
  l_pred_GPP[1] <- log(GPP[1])
  
  ## Process model
  for (j in 2:Ndays) {
    l_pred_GPP[j] = MCMCglmm::rtnorm(1, mean=phi*l_pred_GPP[j-1] + alpha*light[j] + beta*tQ[j], sd = sig_p, lower=-5, upper=5)
  }
  
  for (i in 2:Ndays){
  pred_GPP[i] <- rtnorm(1, mean = exp(l_pred_GPP[i]), sd = obs_err[i], lower=0)
  }
  
  return(pred_GPP)
}

