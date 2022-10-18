##==============================================================================
## Script for within-sample predictions figure
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

# source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Predicted_ProductivityModel_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#1C474D"

##########################
## Model 1 Output - AR
#########################
simmat1_list <- readRDS("./rds files/Sim_6riv_AR_ws_2022_02_27.rds")

# For every day extract median and CI
median_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat1 <- ldply(lapply(simmat1_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat1 <- ldply(df, data.frame)
df_sim1 <- as.data.frame(cbind(dat1$site_name, as.character(dat1$date), dat1$GPP, median_simmat1$X..i.., lower_simmat1$X..i.., upper_simmat1$X..i..))
colnames(df_sim1) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
df_sim1[,3:6] <- apply(df_sim1[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim1 <- left_join(df_sim1, site_info[,c("site_name","short_name")])
df_sim1$short_name <- factor(df_sim1$short_name, levels=site_order_list)

## Plot
df_sim1_plot <- ggplot(df_sim1, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM_AR.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="Within-Sample S-TS Model GPP Predictions")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_AR.col, alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #coord_cartesian(ylim = c(0,30))+
  facet_wrap(~short_name, scales = "free", ncol = 2)
df_sim1_plot

## Save image


###############################
## Model 2 Output - Ricker
###############################
simmat3_list <- readRDS("./rds files/Sim_6riv_Ricker_ws_2022_02_27.rds")

# For every day extract median and CI
median_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat3 <- ldply(df, data.frame)
df_sim3 <- as.data.frame(cbind(dat3$site_name, as.character(dat3$date),
                               dat3$GPP, median_simmat3[,2], lower_simmat3[,2], upper_simmat3[,2]))
colnames(df_sim3) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim3$Date <- as.POSIXct(as.character(df_sim3$Date), format="%Y-%m-%d")
df_sim3[,3:6] <- apply(df_sim3[,3:6],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_sim3 <- left_join(df_sim3, site_info[,c("site_name","short_name")])
df_sim3$short_name <- factor(df_sim3$short_name, levels=site_order_list)

## Plot
df_sim3_plot <- ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM_Ricker.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="Within-Sample LB-TS Ricker Model GPP Predictions")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_Ricker.col, alpha=0.2, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #coord_cartesian(ylim=c(0,5))+
  facet_wrap(~short_name, scales = "free", ncol = 2)

df_sim3_plot

### Overlain
ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  
  geom_line(data=df_sim1, aes(Date, sim_GPP), color=PM_AR.col, size=1.2)+
  geom_ribbon(data=df_sim1, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_AR.col, alpha=0.2, show.legend = FALSE)+
  
  geom_line(aes(Date, sim_GPP), color=PM_Ricker.col, size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_Ricker.col, alpha=0.2, show.legend = FALSE)+

  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  #coord_cartesian(ylim=c(0,max(df_sim1$GPP)*1.5))+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
       title="Within-Sample S-TS & LB-TS Model GPP Predictions")+
  facet_wrap(~short_name, scales = "free", ncol = 2)



################################################
## Plot latent biomass predictions
#################################################

# For every day extract median and CI
median_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) median(x))), data.frame)
lower_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_biomat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[3]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated biomass
df_bio3 <- as.data.frame(cbind(dat3$site_name, as.character(dat3$date),
                               median_biomat3[,2], lower_biomat3[,2], upper_biomat3[,2]))
colnames(df_bio3) <- c("site_name","Date","B","B_lower","B_upper")
df_bio3$Date <- as.POSIXct(as.character(df_bio3$Date), format="%Y-%m-%d")
df_bio3[,3:5] <- apply(df_bio3[,3:5],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_bio3 <- left_join(df_bio3, site_info[,c("site_name","short_name")])
df_bio3$short_name <- factor(df_bio3$short_name, levels=site_order_list)

df_modB3 <- df_bio3

## Add in K estimates
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")
Ricker_pars<-ldply(lapply(stan_model_output_Ricker, function(x) extract(x, c("r","lambda"))), data.frame)
med_pars <- Ricker_pars %>%
  group_by(.id) %>%
  summarise(r.med = median(r, na.rm = TRUE), l.med = median(lambda, na.rm = TRUE))
med_pars$K <- (-1*med_pars$r.med)/med_pars$l.med
# add in info to merge with other data
colnames(med_pars)[which(colnames(med_pars) == ".id")] <- "site_name"
med_pars <- left_join(med_pars, site_info[,c("site_name","short_name")], by="site_name")
med_pars$short_name <- factor(med_pars$short_name, levels=site_order_list)

## plot latent biomass model predictions
df_modB3_plot <- ggplot(df_modB3, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="Within-Sample LB-TS Ricker Model Latent Biomass Predictions")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  coord_cartesian(ylim=c(0,30))+
  geom_hline(data = med_pars, aes(yintercept = K), linetype = "dashed", size=0.9)+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)
df_modB3_plot




###############################
## Model 3 Output - Gompertz
###############################
simmat4_list <- readRDS("./rds files/Sim_9riv_Gompertz_ws.rds")

# For every day extract median and CI
median_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat4 <- ldply(lapply(simmat4_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat4 <- ldply(df, data.frame)
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
                                                          "Clackamas River, OR",
                                                          "Little Difficult Run, VA",
                                                          "Paint Branch Creek, MD",
                                                          "Au Sable River, MI"))

## Plot
df_sim4_plot <- ggplot(df_sim4, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM_Gompertz.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM4: Gompertz")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_Gompertz.col, alpha=0.2, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  coord_cartesian(ylim=c(0,20))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)

df_sim4_plot


## Plot latent biomass model fit
PM4_medpar <- ldply(lapply(stan_model_output_Gompertz,
                           function(x) mechB_extract_medians(rstan::extract(x,c("beta_0","beta_1",
                                                                                "s","c","B",
                                                                                "P","pred_GPP","sig_p")))),
                    data.frame)

df_modB4 <- as.data.frame(cbind(dat4$site_name, as.character(dat4$date), PM4_medpar$B, PM4_medpar$B_Q.025, PM4_medpar$B_Q.975))
colnames(df_modB4) <- c("site_name","Date","B","B_lower","B_upper")
df_modB4$Date <- as.POSIXct(as.character(df_modB4$Date), format="%Y-%m-%d")
df_modB4[,3:5] <- apply(df_modB4[,3:5],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_modB4 <- left_join(df_modB4, site_info[,c("site_name","short_name")])
df_modB4$short_name <- factor(df_modB4$short_name, levels=c("Silver Creek, UT",
                                                            "Medina River, TX",
                                                            "Anacostia River, MD",
                                                            "West Fork River, WV",
                                                            "St. John's River, FL",
                                                            "Clackamas River, OR",
                                                            "Little Difficult Run, VA",
                                                            "Paint Branch Creek, MD",
                                                            "Au Sable River, MI"))

## Plot latent biomass predictions
df_modB4_plot <- ggplot(df_modB4, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM4: Gompertz")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  scale_y_continuous(limits=c(0,25))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)
df_modB4_plot



