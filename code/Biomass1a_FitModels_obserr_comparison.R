## Fitting models to data
## JR Blaszczak

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "parallel","tidyverse","rstan","devtools",
         "bayesplot","shinystan","Metrics","MCMCglmm","tictoc"), require, character.only=T)

## Source data
source("DataSource_9rivers.R")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_l <- lapply(df, function(x) stan_data_compile(x))

#########################################
## Run Stan to get parameter estimates - initial tests
#########################################

## Initial tests
#AR
test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
             data=stan_data_l$nwis_08180700,
             chains=3,iter=5000, control=list(max_treedepth=12))
launch_shinystan(test_ar)

#Ricker
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}
test_ricker <- stan("Stan_ProductivityModel3_Ricker_fixedinit_obserr.stan",
             data=stan_data_l$nwis_08180700,
             init = init_Ricker,
             chains=3,iter=2000, control=list(max_treedepth=12))
launch_shinystan(test_ricker)


#########################################
## Run Stan to get parameter estimates - all sites
#########################################

## PM 1 - Phenomenological
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x,chains=3,iter=5000, control=list(max_treedepth=12))) #5000 seconds
PM_AR_elapsedtime <- lapply(PM_outputlist_AR, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_AR, "./rds files/stan_9riv_output_AR_2021_04_18_obserr.rds")
saveRDS(PM_AR_elapsedtime, "./rds files/stan_9riv_AR_time_2021_04_18_obserr.rds")

## PM 3 - Ricker
init_Ricker <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel3_Ricker_fixedinit_obserr.stan",
                                                data=x,chains=3,iter=5000,init = init_Ricker,
                                                control=list(max_treedepth=12)))
PM_Ricker_elapsedtime <- lapply(PM_outputlist_Ricker, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker, "./rds files/stan_9riv_output_Ricker_2021_04_18_obserr.rds")
saveRDS(PM_Ricker_elapsedtime, "./rds files/stan_9riv_Ricker_time_2021_04_18_obserr.rds")

#############################################
## Compare observation error distributions
############################################

AR_obs <-ldply(lapply(PM_outputlist_AR, function(x) extract(x, c("sig_o"))), data.frame)
Ricker_obs <-ldply(lapply(PM_outputlist_Ricker, function(x) extract(x, c("sig_o"))), data.frame)
data_GPP_sd <- ldply(lapply(df, function(x) return(x$GPP_sd)), data.frame)
colnames(data_GPP_sd) <- c(".id","sig_o")## except this is really GPP_sd from the streamMetabolizer ppd

## compile together
AR_obs$source <- "AR"
Ricker_obs$source <- "Ricker"
data_GPP_sd$source <- "StreamMetabolizer"
obs_dist <- rbind(AR_obs, Ricker_obs, data_GPP_sd)

## Subset to six sites used in final comparison
six_sites <- c("nwis_01649500","nwis_02234000","nwis_03058000",
               "nwis_08180700","nwis_10129900","nwis_14211010")

obs_dist <- obs_dist[which(obs_dist$.id %in% six_sites),]

## add short_name
site_info$.id <- site_info$site_name
obs_dist <- merge(obs_dist, site_info[,c(".id","short_name")], by=".id")

obs_dist$short_name <- factor(obs_dist$short_name, levels=c("Silver Creek, UT",
                                              "Medina River, TX",
                                              "Anacostia River, MD",
                                              "West Fork River, WV",
                                              "St. John's River, FL",
                                              "Clackamas River, OR"))


## Create list to split
obs_dist_l <- split(obs_dist, obs_dist$source)

## plot
ggplot(obs_dist, aes(sig_o, fill=source))+
  geom_histogram(data=obs_dist_l$AR,aes(y=stat(count/sum(count))),alpha=0.7)+
  geom_histogram(data=obs_dist_l$Ricker,aes(y=stat(count/sum(count))),alpha=0.7)+
  geom_histogram(data=obs_dist_l$StreamMetabolizer,aes(y=stat(count/sum(count))),alpha=0.7)+
  theme_bw()+
  facet_wrap(~short_name,ncol=2)+
  scale_y_continuous(limits=c(0,0.12))+
  theme(axis.text = element_text(size=12),
        axis.title = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))








