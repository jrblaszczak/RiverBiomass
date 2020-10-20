## Growth Model 2 - Data simulation

PM2 <- function(r, K, beta_r, sig_o, sig_p, df) {
  
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
  
  ## Vectors for model output of B and pred_GPP
  B<-numeric(Ndays)
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = (B[j-1]*exp(r*B[(j-1)]*(1-(B[(j-1)]/K))) - beta_r)*Q95[j] + beta_r + proc_err
  }
  
  pred_GPP <- light*exp(B) + obs_err
  return(pred_GPP)
}


