## Growth Model 1 - Data simulation 

PM1 <- function(phi, alpha, beta, sig_p, df) {
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  #proc_err <- rlnorm(Ndays, meanlog = 0, sdlog = sig_p) # no longer here, because moved it down to the process model, however process model not working when a log normal distribution
  obs_err <- GPP_sd
  
  ## Vectors for model output
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  l_pred_GPP <- numeric(Ndays)
  l_pred_GPP[1] <- log(GPP[1])
  
  ## Process model
  for (j in 2:Ndays) {
    l_pred_GPP[j] = rnorm(1, mean=phi*l_pred_GPP[j-1] + alpha*light[j] + beta*tQ[j], sd = sig_p)
  }
  
  for (i in 2:Ndays){
  pred_GPP[i] <- rtnorm(1, mean = exp(l_pred_GPP[i]), sd = obs_err[i], lower=0)
  }
  
  return(pred_GPP)
}


