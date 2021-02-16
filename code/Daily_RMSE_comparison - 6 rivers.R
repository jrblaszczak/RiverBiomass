## RMSE comparison figure

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")

# source simulation models
source("Simulated_ProductivityModel1_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p

# for parameter extraction
source("StanParameterExtraction_Source.R")

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"
PM4.col <- "#743731"

## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01649500"="Anacostia River, MD",
                                                                               "nwis_02234000"="St. John's River, FL",
                                                                               "nwis_03058000"="West Fork River, WV",
                                                                               "nwis_08180700"="Medina River, TX",
                                                                               "nwis_10129900"="Silver Creek, UT",
                                                                               "nwis_14211010"="Clackamas River, OR"))


##################################################
## Compare RMSE of all models
###################################################
## extract RMSE
simmat1_list <- readRDS("Sim_6riv_AR.rds")
simmat2_list <- readRDS("Sim_6riv_Logistic.rds")
simmat3_list <- readRDS("Sim_6riv_Ricker.rds")
simmat4_list <- readRDS("Sim_6riv_Gompertz.rds")

## extract only daily rmse
rmsemat1 <- ldply(lapply(simmat1_list, function(x) return(x[[2]])), data.frame)
rmsemat2 <- ldply(lapply(simmat2_list, function(x) return(x[[2]])), data.frame)
rmsemat3 <- ldply(lapply(simmat3_list, function(x) return(x[[2]])), data.frame)
rmsemat4 <- ldply(lapply(simmat4_list, function(x) return(x[[2]])), data.frame)


## Combine
rmse_comp <- as.data.frame(as.matrix(cbind(rmsemat1, rmsemat2$X..i.., rmsemat3$X..i.., rmsemat4$X..i..)))
colnames(rmse_comp) <- c("site_name","PM1: GPP RMSE","PM2: Logistic RMSE",
                         "PM3: Ricker RMSE","PM4: Gompertz RMSE")
rmse_comp$site_name <- as.factor(rmse_comp$site_name)
rmse_comp_long <- melt(rmse_comp, id.vars = "site_name")
rmse_comp_long$value <- as.numeric(as.character(rmse_comp_long$value))

rmse_comp_mean <- rmse_comp_long %>%
  group_by(site_name, variable) %>%
  summarise(rating.mean = mean(na.omit(value)))

rmse_list <- split(rmse_comp_long, rmse_comp_long$site_name)
lapply(rmse_list, function(x) summary(aov(value ~ variable, data=x)))

ggplot(rmse_comp_long, aes(value, fill=as.factor(variable)))+
  geom_boxplot(alpha=0.5)+ coord_flip()+
  facet_wrap(~site_name, ncol = 2)+
  scale_x_continuous(trans = "log")

  scale_fill_manual("",values=c("PM1: GPP RMSE" = PM1.col,"PM2: Logistic RMSE" = PM2.col,
                                "PM3: Ricker RMSE" = PM3.col,"PM4: Gompertz RMSE" = PM4.col))+
  scale_x_continuous(trans="log", limits=c(1,40), breaks = c(1,3,5,10,30), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0.01))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=15), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15),
        legend.position = "none",#c(.5, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=14))+
  labs(y="Density",x="Daily RMSE")+
  facet_wrap(~site_name, scales = "free_x", ncol = 2)




  geom_vline(xintercept = rmse_comp_mean$rating.mean[1],
             color=PM1.col, linetype = "dashed", size = 1)
  geom_vline(xintercept = rmse_comp_mean$rating.mean[2],
             color=PM2.col, linetype = "dashed", size = 1)
  geom_vline(xintercept = rmse_comp_mean$rating.mean[3],
             color=PM3.col, linetype = "dashed", size = 1)
  geom_vline(xintercept = rmse_comp_mean$rating.mean[4],
             color=PM4.col, linetype = "dashed", size = 1)





