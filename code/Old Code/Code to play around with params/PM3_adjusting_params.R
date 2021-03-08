#PM3 <- function(r, lambda, s, c, sig_p, df) {
 
x <- df_Ricker_all[[2]] ## from Productivity model simulations - 10 rivers

r=x$r
lambda=x$lambda
s=x$s
c=x$c
sig_p=x$sig_p
df=x







 
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ # discharge standardized to max value
  
  ## Error
  #proc_err <- rnorm(Ndays, mean = 0, sd = sig_p)
  obs_err <- GPP_sd
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- 0
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper=3.5, lower=-3.5)
  }
  
  plot(B)
  
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = light[i]*exp(B[i]), sd = obs_err[i], lower = 0.01)
  }
  
  plot(pred_GPP)
#}
