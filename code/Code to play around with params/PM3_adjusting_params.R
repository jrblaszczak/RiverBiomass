#PM3 <- function(alpha, gamma, s, c, sig_p, df) {

#Download df from 

  alpha <- 2
  gamma <- 0.3
  s <- 10
  c <- 0.5
  sig_p <- 0.3
  
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- df$light_rel
  tQ <- df$tQ
  
  ## Error
  proc_err <- rlnorm(Ndays, meanlog = 0, sdlog = sig_p)
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
  for (j in 2:Ndays){
    B[j] = rnorm(1, mean = (B[(j-1)] + B[(j-1)]*(b[j] - ant_b[j]*(gamma+(1-gamma)*B[(j-1)])))*P[j], sd = sig_p)
  }
  B
  plot(B)
  plot(exp(B))
  
  
  for (i in 2:Ndays){
    pred_GPP[i] <- rtnorm(1, mean = b[i]*B[i], sd = obs_err[i], lower = 0)
  }
  plot(pred_GPP)
  lines(GPP)
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
#  return(pred_GPP)