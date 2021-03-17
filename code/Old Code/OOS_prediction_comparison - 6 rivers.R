## Out of sample prediction

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2"), require, character.only=T)

# source simulation models
source("Simulated_ProductivityModel1_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

##############################
## Data Import & Processing
##############################
## Import within sample data
source("DataSource_9rivers_oos.R")

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"
PM4.col <- "#743731"

#######################################################################
## Simulate data for each model using the previous parameter estimates
## (modified 'Simulation_matrix_output - 6 rivers.R')
#######################################################################

## Import stan fits - simulate one at a time
stan_model_output_AR <- readRDS("stan_6riv_output_AR.rds")
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker.rds")
stan_model_output_Gompertz <- readRDS("stan_6riv_output_Gompertz.rds")

#####################################
## Model 1 - AR OOS Prediction
####################################
names(dat_oos); names(stan_model_output_AR)
AR_list <- Map(c, stan_model_output_AR, dat_oos)

AR_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars1 <- extract(output, c("phi","alpha","beta","sig_p"))
  simmat1<-matrix(NA,length(df$GPP),length(unlist(pars1$phi)))
  rmsemat1<-matrix(NA,length(df$GPP),1)
  
  # Simulate
  for (i in 1:length(pars1$phi)){
    simmat1[,i]<-PM1(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p[i],df)
    rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat1, rmsemat1)
  return(l)
  
}

AR_sim <- lapply(AR_list, function(x) AR_sim_fxn(x))

## Save simulation
saveRDS(AR_sim, "Sim_6riv_AR_oos.rds")
## If previously simulated
simmat1_list <- readRDS("Sim_6riv_AR_oos.rds")

# For every day extract median and CI
median_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat1 <- ldply(dat_oos, data.frame)
df_sim1 <- as.data.frame(cbind(dat1$site_name, as.character(dat1$date), dat1$GPP, median_simmat1$X..i.., lower_simmat1$X..i.., upper_simmat1$X..i..))
colnames(df_sim1) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
df_sim1[,3:6] <- apply(df_sim1[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim1 <- left_join(df_sim1, site_info[,c("site_name","short_name")])
df_sim1$short_name <- factor(df_sim1$short_name, levels=c("Silver Creek, UT",
                                                          "Medina River, TX",
                                                          "Anacostia River, MD",
                                                          "West Fork River, WV",
                                                          "St. John's River, FL",
                                                          "Clackamas River, OR"))
## Plot
df_sim1_plot <- ggplot(df_sim1, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM1.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: GPP Out-of-Sample Prediction")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM1.col, alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #scale_y_continuous(limits=c(0,30))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)
df_sim1_plot


###############################
## Model 3 Output - Ricker
###############################
names(dat_oos); names(stan_model_output_Ricker)
Ricker_list <- Map(c, stan_model_output_Ricker, dat_oos[(names(stan_model_output_Ricker))])

Ricker_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars3<-extract(output, c("r","lambda","s","c","B","P","pred_GPP","sig_p"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM3(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],df)
    biomat3[,i]<-PM3_B(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim <- lapply(Ricker_list, function(x) Ricker_sim_fxn(x))

## Save simulation
saveRDS(Ricker_sim, "./rds files/Sim_6riv_Ricker_oos.rds")








## If previously simulated
simmat3_list <- readRDS("Sim_6riv_Ricker_oos.rds")


# For every day extract median and CI
median_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat3 <- ldply(dat_oos, data.frame)
df_sim3 <- as.data.frame(cbind(dat3$site_name, as.character(dat3$date),
                               dat3$GPP, median_simmat3$X..i.., lower_simmat3$X..i.., upper_simmat3$X..i..))
colnames(df_sim3) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim3$Date <- as.POSIXct(as.character(df_sim3$Date), format="%Y-%m-%d")
df_sim3[,3:6] <- apply(df_sim3[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim3 <- left_join(df_sim3, site_info[,c("site_name","short_name")])
df_sim3$short_name <- factor(df_sim3$short_name, levels=c("Silver Creek, UT",
                                                          "Medina River, TX",
                                                          "Anacostia River, MD",
                                                          "West Fork River, WV",
                                                          "St. John's River, FL",
                                                          "Clackamas River, OR"))

## Plot
df_sim3_plot <- ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM3.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM3: Ricker Out-of-Sample Prediction")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM3.col, alpha=0.2, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)

df_sim3_plot


###############################
## Model 4 Output - Gompertz
###############################
names(dat_oos); names(stan_model_output_Gompertz)
Gompertz_list <- Map(c, stan_model_output_Gompertz, dat_oos)

Gompertz_sim_fxn <- function(x){
  #separate data
  output <- x[[1]]
  df <- x
  
  # extract
  pars5<-extract(output, c("beta_0","beta_1","s","c","B","P","pred_GPP","sig_p"))
  simmat5<-matrix(NA,length(df$GPP),length(unlist(pars5$sig_p)))
  rmsemat5<-matrix(NA,length(df$GPP),1)
  #Simulate
  for (i in 1:length(pars5$beta_0)){
    simmat5[,i]<-PM5(pars5$beta_0[i],pars5$beta_1[i],pars5$s[i],pars5$c[i],pars5$sig_p[i],df)
    rmsemat5[i]<-sqrt(sum((simmat5[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat5, rmsemat5)
  return(l)
  
}

Gompertz_sim <- lapply(Gompertz_list, function(x) Gompertz_sim_fxn(x))

## Save simulation
saveRDS(Gompertz_sim, "Sim_6riv_Gompertz_oos.rds")
## If previously simulated
simmat4_list <- readRDS("Sim_6riv_Gompertz_oos.rds")

# For every day extract median and CI
median_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat4 <- ldply(dat_oos, data.frame)
df_sim4 <- as.data.frame(cbind(dat4$site_name, as.character(dat4$date),
                               dat4$GPP, median_simmat4$X..i.., lower_simmat4$X..i.., upper_simmat4$X..i..))
colnames(df_sim4) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim4$Date <- as.POSIXct(as.character(df_sim4$Date), format="%Y-%m-%d")
df_sim4[,3:6] <- apply(df_sim4[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim4 <- left_join(df_sim4, site_info[,c("site_name","short_name")])
df_sim4$short_name <- factor(df_sim4$short_name, levels=c("Silver Creek, UT",
                                                          "Medina River, TX",
                                                          "Anacostia River, MD",
                                                          "West Fork River, WV",
                                                          "St. John's River, FL",
                                                          "Clackamas River, OR"))

## Plot
df_sim4_plot <- ggplot(df_sim4, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM4.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM4: Gompertz Out-of-Sample Prediction")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM4.col, alpha=0.2, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)

df_sim4_plot



