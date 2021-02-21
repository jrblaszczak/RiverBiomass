## Growth Model 3 - Data simulation

#PM3 <- function(r, lambda, s, c, sig_p, df) {

pars3<-extract(stan_model_output_Ricker[[1]], c("r","lambda","s","c","B","P","pred_GPP","sig_p"))

df <- dat$nwis_01645762
r <- pars3$r[50]
lambda <- pars3$lambda[50]
s <- pars3$s[50]
c <- pars3$c[50]
sig_p <- pars3$sig_p[50]


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
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] = rnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = obs_err[i], lower=0.001)
  }
 
plot(exp(B))  
plot(pred_GPP)
