## Growth Model 1 - Data simulation

PM1 <- function(phi, alpha, beta_r, sig_o, sig_p, df) {
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  light <- df$light_rel
  tQ <- df$tQ
  Q95 <- df$Q95
  
  ## Error
  proc_err <- rnorm(Ndays, mean = 0, sd = sig_p)
  obs_err<-numeric(length(Ndays))
  for(i in 1:Ndays){
    obs_err[i] = rnorm(Ndays, mean=0, sd = sig_o)
  }
  
  ## Vectors for model output
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  l_pred_GPP <- numeric(Ndays)
  l_pred_GPP[1] <- log(GPP[1])
  
  ## Process model
  for (j in 2:Ndays) {
    l_pred_GPP[j] = (phi*l_pred_GPP[j-1] + alpha*light[j] - beta_r)*Q95[j] + beta_r + proc_err
  }
  pred_GPP <- exp(l_pred_GPP) + obs_err
  return(pred_GPP)
}


