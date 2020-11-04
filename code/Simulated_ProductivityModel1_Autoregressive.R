## Growth Model 1 - Data simulation

PM1 <- function(phi, alpha, beta, sig_p, df) {
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ
  Q95 <- df$Q95
  
  ## Error
  proc_err <- rnorm(Ndays, mean = 0, sd = sig_p)
  obs_err <- GPP_sd
  
  ## Vectors for model output
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  l_pred_GPP <- numeric(Ndays)
  l_pred_GPP[1] <- log(GPP[1])
  
  ## Process model
  for (j in 2:Ndays) {
    l_pred_GPP[j] = phi*l_pred_GPP[j-1] + alpha*light[j] + beta*tQ[j] + proc_err
  }
  pred_GPP <- exp(l_pred_GPP) + obs_err
  return(pred_GPP)
}


