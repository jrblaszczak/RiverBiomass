#PM2 <- function(r, K, s, c, sig_p, df) {

i <- 2
#r=pars2$r[i]
r=0.3
K=pars2$K[i]
s=pars2$s[i]
c=pars2$c[i]
sig_p=pars2$sig_p[i]


## Data
Ndays<-length(df$GPP)
GPP <- df$GPP
GPP_sd <- df$GPP_sd
light <- df$light_rel
tQ <- df$tQ # discharge standardized to max value
depth <- df$depth
turb <- df$turb

## Error
obs_err <- GPP_sd

## Vectors for model output of P, B, pred_GPP
P <- numeric(Ndays)
P[1] <-1
for(i in 2:length(tQ)){
  P[i] = exp(-exp(s*(tQ[i] - c)))
}
plot(P)

## Benthic light
#ben_light <- numeric(Ndays)
#for(i in 1:length(light)){
#  ben_light[i]=light[i]*exp(-1*a*turb[i]*depth[i])
#}

B<-numeric(Ndays)
B[1] <- 0
pred_GPP<-numeric(Ndays)
pred_GPP[1] <- light[1]*exp(B[1])

## Process Model
for (j in 2:Ndays){
  B[j] = rtnorm(1, mean = (B[j-1]*exp(r*B[(j-1)]*(1-(B[(j-1)]/K))))*P[j], sd = sig_p, upper=3.5)
}
B


for (i in 2:Ndays){
  pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = obs_err[i], lower = 0)
}

plot(pred_GPP)
   
#return(exp(pred_GPP))
#}