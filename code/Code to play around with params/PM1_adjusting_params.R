## Growth Model 1 - Data simulation 

#PM1 <- function(phi, alpha, beta, sig_p, df) {
library(MCMCglmm) ## this is needed to use the function rtnorm (there are other packages)

phi <- 0.95
alpha <- 0.02
beta <- -0.01
sig_p <- 0.3

  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  #proc_err <- rlnorm(Ndays, meanlog = 0, sdlog = sig_p)
  obs_err <- GPP_sd
  
  ## Vectors for model output
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- log(GPP[1])
  l_pred_GPP <- numeric(Ndays)
  l_pred_GPP[1] <- log(GPP[1])
  
  ## Process model
  for (j in 2:Ndays) {
    l_pred_GPP[j] = rnorm(1, mean=phi*l_pred_GPP[j-1] + alpha*light[j] + beta*tQ[j], sd = sig_p) ## When change to rlnorm, it goes crazy
  }
  plot(l_pred_GPP)
  l_pred_GPP
  
  
  0.95*l_pred_GPP[2]
  0.02*light[2]
  (-0.01)*tQ[2]
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = exp(l_pred_GPP[i]), sd = obs_err[i], lower=0)
  }
  plot(pred_GPP)
  lines(GPP)
  
#  return(pred_GPP)
#}


