## Observation error comparison
## JR Blaszczak

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc"), require, character.only=T)

## Source data
source("DataSource_6rivers.R")

# colors
PM_AR.col <- "#d95f02"
PM_Ricker.col <- "#7570b3"
PM_Gompertz.col <- "#3C7373"

## Import stan fits - simulate one at a time
stan_model_output_AR <- readRDS("./rds files/stan_6riv_output_AR_2021_05_17.rds")
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2021_05_17.rds")
#stan_model_output_Gompertz <- readRDS("./rds files/stan_6riv_output_Gompertz_2021_05_16.rds")

#############################################
## Compare observation error distributions
############################################
AR_obs <-ldply(lapply(stan_model_output_AR, function(x) extract(x, c("sig_o"))), data.frame)
Ricker_obs <-ldply(lapply(stan_model_output_Ricker, function(x) extract(x, c("sig_o"))), data.frame)
#Gompertz_obs <-ldply(lapply(stan_model_output_Gompertz, function(x) extract(x, c("sig_o"))), data.frame)
data_GPP_sd <- ldply(lapply(df, function(x) return(x$GPP_sd)), data.frame)
colnames(data_GPP_sd) <- c(".id","sig_o")## except this is really GPP_sd from the streamMetabolizer ppd

## compile together
AR_obs$source <- "STS"
Ricker_obs$source <- "LB - Ricker"
#Gompertz_obs$source <- "LB - Gompertz"
data_GPP_sd$source <- "StreamMetabolizer"
obs_dist <- rbind(AR_obs, Ricker_obs, 
                  #Gompertz_obs, 
                  data_GPP_sd)

## add short_name
site_info$.id <- site_info$site_name
obs_dist <- merge(obs_dist, site_info[,c(".id","short_name")], by=".id")

obs_dist$short_name <- factor(obs_dist$short_name, levels= site_order_list)
obs_dist$source <- factor(obs_dist$source, levels= c("StreamMetabolizer","STS","LB - Ricker")) #,"LB - Gompertz"))

## Create list to split
obs_dist_l <- split(obs_dist, obs_dist$source)

## plot
ggplot(obs_dist, aes(sig_o, fill=source))+
  geom_histogram(data=obs_dist_l$AR,aes(y=stat(count/sum(count))),alpha=0.5, color="black")+
  geom_histogram(data=obs_dist_l$Ricker,aes(y=stat(count/sum(count))),alpha=0.5, color="black")+
  #geom_histogram(data=obs_dist_l$Gompertz,aes(y=stat(count/sum(count))),alpha=0.5, color="black")+
  geom_histogram(data=obs_dist_l$StreamMetabolizer,aes(y=stat(count/sum(count))),alpha=0.7, color="black")+
  theme_bw()+
  facet_wrap(~short_name,ncol=2)+
  scale_fill_manual("Source", values = c("STS" = PM_AR.col,
                                "LB - Ricker" = PM_Ricker.col,
                                #"LB - Gompertz" = PM_Gompertz.col,
                                "StreamMetabolizer" = "#BABDBF"))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=14), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=12))+
  labs(x=expression(sigma['obs']),y="Proportion of Total")








