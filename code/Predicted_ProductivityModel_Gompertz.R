## Growth Model Gompertz - Data simulation

PM_Gompertz <- function(beta_0, beta_1, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value

  
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
    B[j] = MCMCglmm::rtnorm(1, mean = (beta_0 + beta_1*B[(j-1)])*P[j], sd = sig_p, upper=5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0)
  }
  
  return(pred_GPP)
  
}


PM_Gompertz_B <- function(beta_0, beta_1, s, c, sig_p, sig_o, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value

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
    B[j] = MCMCglmm::rtnorm(1, mean = (beta_0 + beta_1*B[(j-1)])*P[j], sd = sig_p, upper=5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0)
  }
  
  return(B)
  
}



