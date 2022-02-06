## Growth Model 3 - Data simulation

PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  #GPP_sd <- df$GPP_sd
  light <- df$light_rel_PAR
  tQ <- df$tQ # discharge standardized to max value

  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}



## Biomass
PM_Ricker_B <- function(r, lambda, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel_PAR
  tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(B)
}



PM_Ricker_lv <- function(r, lambda, s, c, sig_p, sig_o, df, light_version) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- light_version
  tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

