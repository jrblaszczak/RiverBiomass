## Growth Model 5 - Data simulation

PM5 <- function(beta_0, beta_1, s, c, sig_p, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  obs_err <- GPP_sd
  
  ## Vectors for model output of P, B, pred_GPP, r
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])

  ## Process Model
  for (j in 2:Ndays){
    B[j] = norm(1, mean = (beta_0 + beta_1*B[(j-1)])*P[j], sd = sig_p)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- norm(1, mean = light[i]*exp(B[i]), sd = obs_err[i])
  }
  
  return(pred_GPP)
  
}

