#PM2 <- function(r, K, s, c, sig_p, df) {

r <- 0.2
K <- 6
s <- 10
c <- 0.5
sig_p <- 0.3


  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  proc_err <- rnorm(Ndays, mean = 0, sd = sig_p)
  obs_err <- GPP_sd
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  plot(P)
  
  B<-numeric(Ndays)
  B[1] <- 0
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = rnorm(1, mean = (B[j-1]*exp(0.3*B[(j-1)]*(1-(B[(j-1)]/1))))*P[j], sd = sig_p)
  }
  B
  plot(B)
  plot(exp(B))
  

  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = light[i]*exp(B[i]), sd = obs_err[i]*.2, lower = 0)
  }
  plot(pred_GPP)
  lines(GPP)
 
plot(exp(pred_GPP))
   
#return(exp(pred_GPP))
#}