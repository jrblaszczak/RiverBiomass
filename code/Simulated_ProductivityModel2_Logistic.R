## Growth Model 2 - Data simulation

PM2 <- function(rmax, K, s, c, sig_p, df) {
  
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
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = (B[j-1]*exp(rmax*B[(j-1)]*(1-(B[(j-1)]/K))))*P[j] + rnorm(1,0,proc_err)
  }
  
  pred_GPP <- light*exp(B) + rnorm(1,0,obs_err)
  return(pred_GPP)
}


