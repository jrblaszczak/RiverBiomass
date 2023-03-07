##==============================================================================
## Script for within-sample fit figure
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")
## Source functions to extract and summarize parameter estimates
source("StanParameterExtraction_Source.R")

## Create data frame
dat <- ldply(df, data.frame)

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"

###################################
## Plot within sample fit for S-TS
###################################
## Import stan model fit
stan_model_output_STS <- readRDS("./rds files/stan_6riv_output_AR_2022_02_22.rds")

PM1_medpar <- ldply(lapply(stan_model_output_STS,
                           function(x) STS_extract_medians(rstan::extract(x,c("l_pred_GPP")))),
                    data.frame)

## GPP fit
sts_GPP <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, PM1_medpar$pred_GPP, PM1_medpar$pred_GPP_Q.025, PM1_medpar$pred_GPP_Q.975))
colnames(sts_GPP) <- c("site_name","Date","GPP","pred_GPP","pred_GPP_lower","pred_GPP_upper")
sts_GPP$Date <- as.POSIXct(as.character(sts_GPP$Date), format="%Y-%m-%d")
sts_GPP[,3:6] <- apply(sts_GPP[,3:6],2,function(x) as.numeric(as.character(x)))
## Arrange rivers by river order
sts_GPP <- left_join(sts_GPP, site_info[,c("site_name","short_name")])
sts_GPP$short_name <- factor(sts_GPP$short_name, levels=site_order_list)

## Plot GPP pred
ggplot(sts_GPP, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, pred_GPP), color=PM_AR.col, size=0.75)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="Within-Sample S-TS Model Fit")+
  geom_ribbon(aes(ymin=pred_GPP_lower,ymax=pred_GPP_upper),
              fill=PM_AR.col, alpha=0.5, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #coord_cartesian(ylim=c(0,5))+
  facet_wrap(~short_name, scales = "free", ncol = 2)



###################################
## Plot within sample fit for LB-TS
###################################
## Import stan Ricker model fit
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")


PM3_medpar <- ldply(lapply(stan_model_output_Ricker,
                           function(x) LBTS_extract_medians(rstan::extract(x,c("B","P","pred_GPP")))),
                    data.frame)

## GPP fit
lbts_GPP <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, PM3_medpar$pred_GPP, PM3_medpar$pred_GPP_Q.025, PM3_medpar$pred_GPP_Q.975))
colnames(lbts_GPP) <- c("site_name","Date","GPP","pred_GPP","pred_GPP_lower","pred_GPP_upper")
lbts_GPP$Date <- as.POSIXct(as.character(lbts_GPP$Date), format="%Y-%m-%d")
lbts_GPP[,3:6] <- apply(lbts_GPP[,3:6],2,function(x) as.numeric(as.character(x)))
## Arrange rivers by river order
lbts_GPP <- left_join(lbts_GPP, site_info[,c("site_name","short_name")])
lbts_GPP$short_name <- factor(lbts_GPP$short_name, levels=site_order_list)

## Plot GPP pred
ggplot(lbts_GPP, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, pred_GPP), color=PM_Ricker.col, size=0.75)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="Within-Sample LB-TS Model Fit")+
  geom_ribbon(aes(ymin=pred_GPP_lower,ymax=pred_GPP_upper),
              fill=PM_Ricker.col, alpha=0.5, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #coord_cartesian(ylim=c(0,5))+
  facet_wrap(~short_name, scales = "free", ncol = 2)



## Latent Biomass fit
df_modB3 <- as.data.frame(cbind(dat$site_name, as.character(dat$date), PM3_medpar$B, PM3_medpar$B_Q.025, PM3_medpar$B_Q.975))
colnames(df_modB3) <- c("site_name","Date","B","B_lower","B_upper")
df_modB3$Date <- as.POSIXct(as.character(df_modB3$Date), format="%Y-%m-%d")
df_modB3[,3:5] <- apply(df_modB3[,3:5],2,function(x) as.numeric(as.character(x)))
## Arrange rivers by river order
df_modB3 <- left_join(df_modB3, site_info[,c("site_name","short_name")])
df_modB3$short_name <- factor(df_modB3$short_name, levels=site_order_list)

## Plot latent biomass predictions
ggplot(df_modB3, aes(Date, exp(B)))+
  geom_line(size=0.75, color="chartreuse4")+
  labs(y="Latent Biomass",title="Within-Sample LB-TS Model Estimates")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.5, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)

