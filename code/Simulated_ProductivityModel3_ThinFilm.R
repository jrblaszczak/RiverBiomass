## Growth Model 3 - Thin Film (Density Independent Growth)

PM3 <- function(alpha, gamma, s, c, sig_p, df) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ
  Q95 <- df$Q95
  
  ## Error
  proc_err <- rnorm(Ndays, mean = 0, sd = sig_p)
  obs_err<-GPP_sd
  
  ## Vectors for model output
  P <- numeric(Ndays)
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  pred_GPP<-numeric(Ndays)
  
  b <- alpha*light
  d <- numeric(Ndays)
  B[1] <- GPP[1]/b[1]
  
  ant_b <- numeric(Ndays)
  for (j in 1:Ndays){
    if (j < 21){
      ant_b[j] = mean(b[1:j])
    } else { ant_b[j] = mean(b[(j-20):j])}
  }
  
  ## Process Model
  for (j in 2:Ndays) {
    
    B[j] = (B[(j-1)] + B[(j-1)]*(b[j] - ant_b[j]*(gamma+(1-gamma)*B[(j-1)])))*P[j] + rnorm(1,0,proc_err)
  }
  
  pred_GPP <- rtnorm(1,b*B,obs_err)
  
  return(pred_GPP)
  
}
