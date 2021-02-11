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
  B[1] <- 0
  pred_GPP<-numeric(Ndays)
  #pred_GPP[1] <- exp(B[1])*light[j]

  ## Process Model
  for (j in 2:Ndays){
    B[j] = rtnorm(1, mean = (beta_0 + beta_1*B[(j-1)])*P[j], sd = sig_p, upper=3.5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = obs_err[i], lower = 0)
  }
  
  return(pred_GPP)
}



PM5_BL <- function(beta_0, beta_1, s, c, a, sig_p, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  depth <- df$depth
  turb <- df$turb
  
  ## Error
  obs_err <- GPP_sd
  
  ## Vectors for model output of P, B, pred_GPP, r
  P <- numeric(Ndays)
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  
  ## Benthic light
  ben_light <- numeric(Ndays)
  for(i in 1:length(light)){
    ben_light[i]=light[i]*exp(-1*a*turb[i]*depth[i])
  }
  
  B<-numeric(Ndays)
  B[1] <- 0
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- exp(B[1])*ben_light[j]
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = rtnorm(1, mean = (beta_0 + beta_1*B[(j-1)])*P[j], sd = sig_p, upper=3.5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = exp(B[i]), sd = obs_err[i], lower = 0)
  }
  
  return(pred_GPP)
}
