## Comparison of Ricker with and without fixed initial biomass

## See simulation_matrix_output.R for code to create figures
# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","rstan","shinystan","MCMCglmm"), require, character.only=T)

#################
## Source data
#################
source("DataSource_9rivers.R")
df <- dat

####################
## Stan data prep 
####################
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ, B_int=log(x$GPP[1]/x$light_rel[1]))
  return(data)
}

stan_data_l <- lapply(dat, function(x) stan_data_compile(x))

########################
## Run Stan
########################

## Limit to just 2 sites for comparison
stan_data_l <- stan_data_l[8:9]

## Ricker - free initial biomass
output_Ricker <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel3_Ricker.stan",
                                                             data=x,chains=3,iter=5000,
                                                             control=list(max_treedepth=12)))

## Ricker - fixed initial biomass
output_Ricker_fixedinit <- lapply(stan_data_l, function(x) stan("Stan_ProductivityModel3_Ricker_fixedinit.stan",
                                                      data=x,chains=3,iter=5000,
                                                      control=list(max_treedepth=12)))


########################
## Compare predictions
########################
## Source deterministic model
source("Predictions_ProductivityModel3_Ricker.R")

## Change the following between output_Ricker and compare to output_Ricker_fixedinit
stan_model_output_Ricker <- output_Ricker

## Combine data with model output
names(dat); names(stan_model_output_Ricker)
Ricker_list <- Map(c, stan_model_output_Ricker, dat[(names(stan_model_output_Ricker))])

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

Ricker_sim <- lapply(Ricker_list, function(x) Ricker_sim_fxn(x)) ## this will take some time

## Save simulation
saveRDS(Ricker_sim, "Sim_9riv_Ricker.rds")
## If previously simulated
simmat3_list <- readRDS("Sim_9riv_Ricker.rds")
simmat3_list <- Ricker_sim

# For every day extract median and CI
median_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
lower_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
upper_simmat3 <- ldply(lapply(simmat3_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)

## Plot simulated GPP
dat3 <- ldply(dat[names(Ricker_sim)], data.frame)
df_sim3 <- as.data.frame(cbind(dat3$site_name, as.character(dat3$date),
                               dat3$GPP, median_simmat3[,2], lower_simmat3[,2], upper_simmat3[,2]))
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
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM3: Ricker")+
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










