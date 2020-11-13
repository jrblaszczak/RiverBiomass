## Growth Model 1 - Data simulation 

#PM1 <- function(phi, alpha, beta, sig_p, df) {

phi <- 0.95
alpha <- 2
beta <- -0.5
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
    l_pred_GPP[j] = rnorm(1, mean=0.9*l_pred_GPP[j-1] + 2*light[j] + (-2)*tQ[j], sd = sig_p) ## on another scale?
  }
  plot(l_pred_GPP)
  l_pred_GPP
  
  
  0.9*l_pred_GPP[1]
  2*light[2]
  (-2)*tQ[2]
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = exp(l_pred_GPP[i]), sd = obs_err[i], lower=0)
  }
  plot(pred_GPP)
  lines(GPP)
  
#  return(pred_GPP)
#}


