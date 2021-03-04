# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse", "reshape2","PerformanceAnalytics",
         "rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("Oregon_ClackamasData_Source.R")

# source simulation models
# input variables: GPP, GPP_sd, light, tQ
source("Simulated_ProductivityModel1_Autoregressive.R") # estimated parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # estimated parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_ThinFilm.R") # estimated parameters: alpha, gamma, s, c, sig_p

# subset the data
df <- dat$Clackamas_OR
df <- df[which(df$year == "2011"),]
df <- df[1:170,]
plot(df$GPP)

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP, GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data <- stan_data_compile(df)

###############################
## Run Stan to get estimates ##
###############################
PM1_DataOutput <- stan("Stan_ProductivityModel1_Autoregressive.stan", data=stan_data, chains=4, iter=1000)
PM2_DataOutput <- stan("Stan_ProductivityModel2_Logistic.stan", data=stan_data, chains=4, iter=1000)
PM3_DataOutput <- stan("Stan_ProductivityModel3_ThinFilm.stan", data=stan_data, chains=4, iter=1000)

PM1_extract <- extract(PM1_DataOutput, c("phi","alpha","beta","l_pred_GPP","sig_p"))
PM2_extract <- extract(PM2_DataOutput, c("r","K","s","c","B","P","pred_GPP","sig_p"))
PM3_extract <- extract(PM3_DataOutput, c("alpha","gamma","s","c","N","P","pred_GPP","sig_p"))

source("StanParameterExtraction_Source.R")
PM1_medpar <- phenom_extract_medians(extract(PM1_DataOutput, c("phi","alpha","beta","l_pred_GPP","sig_p")))
PM2_medpar <- mechB_extract_medians(extract(PM2_DataOutput, c("r","K","s","c","B","P","pred_GPP","sig_p")))
PM3_medpar <- mechN_extract_medians(extract(PM3_DataOutput, c("alpha","gamma","s","c","N","P","pred_GPP","sig_p")))

###################
## Simulate Data ##
###################
## Simulate data for each growth model
df$simGPP_PM1 <- PM1(phi=PM1_medpar$par$phi,
                     alpha=PM1_medpar$par$alpha,
                     beta=PM1_medpar$par$beta,
                     sig_p=PM1_medpar$par$sig_p, df=df)
df$simGPP_PM2 <- PM2(r=PM2_medpar$par$r,
                     K=PM2_medpar$par$K,
                     s=PM2_medpar$par$s,
                     c=PM2_medpar$par$c,
                     sig_p=PM2_medpar$par$sig_p, df=df)
df$simGPP_PM3 <- PM3(alpha=PM3_medpar$par$alpha,
                     gamma=PM3_medpar$par$gamma,
                     s=PM3_medpar$par$s,
                     c=PM3_medpar$par$c,
                     sig_p=PM3_medpar$par$sig_p, df=df)










