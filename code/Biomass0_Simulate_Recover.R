##==============================================================================
## Script to simulate data and test parameter recovery
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools","shinystan",
         "MCMCglmm"), require, character.only=T)

theme_set(theme(legend.position = "none",
                  panel.background = element_rect(color = "black", fill=NA, size=1),
                  axis.text = element_text(size=12),
                  axis.title = element_text(size=12)))

#############################
## Source & visualize data
#############################
source("DataSource_6rivers_StreamLight.R")

## Choose example data
ex <- df$nwis_01608500
ex2 <- df$nwis_01649190
  
## Visualize GPP
plot_grid(
  ggplot(ex, aes(date, GPP))+
    geom_line(color="chartreuse4")+
    geom_point(color="chartreuse4", size=2)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="South Branch Potomac, WV"),
  
  ggplot(ex2, aes(date, GPP))+
    geom_line(color="chartreuse4")+
    geom_point(color="chartreuse4", size=2)+
    labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title="Paint Branch, MD"),
  ncol=2)

## Visualize light_rel differences
plot_grid(
  ggplot(ex, aes(date, light_rel_PPFD))+
    geom_line()+
    geom_line(aes(date, light_rel_PAR), color="blue")+
    labs(title="South Branch Potomac, WV"),
  ggplot(ex2, aes(date, light_rel_PPFD))+
    geom_line()+
    geom_line(aes(date, light_rel_PAR), color="blue")+
    labs(title="Paint Branch, MD"),
  ncol=2
)


## Simulate data using deterministic function
PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df, light_version) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- light_version
  tQ <- df$tQ # discharge standardized to max value
  
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
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

## Predict nwis_01608500 using previously fit parameters (easy example)
ex_pred.GPP.PAR <- PM_Ricker(r = 0.3, lambda = -0.03,
                         s = 124, c = 0.28,
                         sig_p = 0.21, sig_o = 0.9,
                         df = ex, light_version = ex$light_rel_PAR)
ex_pred.GPP.PPFD <- PM_Ricker(r = 0.3, lambda = -0.03,
                             s = 124, c = 0.28,
                             sig_p = 0.21, sig_o = 0.9,
                             df = ex, light_version = ex$light_rel_PPFD)
## Plot comparison
plot(ex$GPP, pch=19)
lines(ex_pred.GPP.PPFD, col="blue")
lines(ex_pred.GPP.PAR, col="red") #does better job but parameters are chosen based on fit


## Predict nwis_01649190 using previously fit parameters (challenging example)
ex_pred.GPP.PAR <- PM_Ricker(r = 0.3, lambda = -0.03,
                             s = 124, c = 0.28,
                             sig_p = 0.21, sig_o = 0.9,
                             df = ex, light_version = ex$light_rel_PAR)
ex_pred.GPP.PPFD <- PM_Ricker(r = 0.3, lambda = -0.03,
                              s = 124, c = 0.28,
                              sig_p = 0.21, sig_o = 0.9,
                              df = ex, light_version = ex$light_rel_PPFD)
## Plot comparison
plot(ex_pred.GPP.PAR, type="l")
lines(ex_pred.GPP.PPFD, col="blue")















