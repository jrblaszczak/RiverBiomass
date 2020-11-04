## Growth Model 3 - Thin Film (Density Independent Growth)

PM3 <- function(alpha, gamma, beta_r, sig_o, sig_p, df) {
  
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
    
    B[j] = (B[(j-1)] + B[(j-1)]*(b[j] - ant_b[j]*(gamma+(1-gamma)*B[(j-1)])) - beta_r)*Q95[j] + beta_r + proc_err
    
  }
    pred_GPP <- b*B + obs_err
  return(pred_GPP)
  
}
