PM5(pars5$beta_0[i],pars5$beta_1[i],pars5$beta_2[i],pars5$s[i],pars5$c[i],pars5$sig_p[i],df)

#PM5 <- function(beta_0, beta_1, beta_2, s, c, sig_p, df) {
i <- 2
beta_0 <- pars5$beta_0[i]
beta_1 <- pars5$beta_1[i]
beta_2 <- pars5$beta_2[i]
s <- pars5$s[i]
c <- pars5$c[i]
sig_p <- pars5$sig_p[i]



  
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
  pred_GPP[1] <- exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = rtnorm(1, mean = (beta_0 + beta_1*exp(B[(j-1)])+beta_2*light[j])*P[j], sd = sig_p, upper=3.5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = exp(B[i]), sd = obs_err[i], lower = 0)
  }
  
  plot(pred_GPP)
#}