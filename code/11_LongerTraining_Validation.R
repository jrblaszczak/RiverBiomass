##==============================================================================
## Script for fitting stan models to more years of training data for predictions
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Read in data files
Pot_TS <- readRDS("./rds files/SBPotomac_longTS.rds")

## Read in streamlight
Pot_SL <- readRDS("./rds files/SBPotomac_SL.rds")

#######################################################################
## Longer time series data prep - South Branch Potomac River, WV
#######################################################################

## Join TS data and StreamLight
colnames(Pot_SL)[colnames(Pot_SL) == "Date"] <- "date"
data <- left_join(Pot_TS, Pot_SL, by=c("site_name", "date"))

## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
data_siteyears$year
# SB Potomac = 2008 2010 2012 2013 2014 2015 2016 (previous out-of-sample year was 2013 predicted from 2012)

## Set any GPP < 0 to a small value between 0.05 to 0.13 g O2 m-2 d-1
data[which(data$GPP < 0),]$GPP <- sample(exp(-3):exp(-2), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

# remove 2013
Pot_longtrain <- data[which(data$year %in% c("2012","2014","2015","2016")),]
Pot_shorttrain <- data[which(data$year == "2012"),]

## Relativize light and discharge
rel_LQT <- function(x){
  x$light_rel_PPFD <- x$light/max(x$light)
  x$light_rel_PAR <- x$PAR_surface/max(x$PAR_surface)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  x<-x[order(x$date),]
  return(x)
}

Pot_longtrain <- rel_LQT(Pot_longtrain)
Pot_shorttrain <- rel_LQT(Pot_shorttrain)

## visualize years included
ggplot(Pot_longtrain, aes(date, GPP))+
  geom_point()+
  geom_point(data = Pot_shorttrain, aes(date, GPP), color = "blue")+
  geom_point(data = data[which(data$year == "2013"),], aes(date, GPP), color = "purple")

plot_grid(
  ggplot(Pot_longtrain, aes(date, light_rel_PAR))+
  geom_point(),
  ggplot(Pot_longtrain, aes(date, tQ))+
    geom_point(),
  ncol= 1)

ggplot(Pot_longtrain, aes(tQ))+
  geom_histogram()+
  geom_histogram(data = Pot_shorttrain, aes(tQ), fill = "blue")

####################
## Stan data prep ##
####################
rstan_options(auto_write=TRUE)
options(mc.cores=8)

stan_data_compile <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_Pot_long <- stan_data_compile(Pot_longtrain)
stan_data_Pot_short <- stan_data_compile(Pot_shorttrain)

##############################################################
## Run Stan to get parameter estimates - initial tests
##############################################################

## Initial tests
#AR
test_ar <- stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
             data=stan_data_Pot_short,
             chains=4,iter=5000, 
             control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_ar)

#Ricker - P reparameterized
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

test_ricker <- stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                    data=stan_data_Pot_short,
                    init = init_Ricker,chains=4,iter=5000,
                    control=list(max_treedepth=12, adapt_delta=0.95))
launch_shinystan(test_ricker)

rm(test_ar, test_ricker)


###################################################
## Run Stan to get parameter estimates - all sites
###################################################
stan_data_l <- list(stan_data_Pot_long, stan_data_Pot_short)

## PM 1 - Standard time series
PM_outputlist_AR <- lapply(stan_data_l,
                           function(x) rstan::stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                                   data=x, chains=4, iter=5000,
                                                   control=list(max_treedepth=12, adapt_delta=0.95)))
names(PM_outputlist_AR) <- c("Pot_long","Pot_short")
saveRDS(PM_outputlist_AR, "./rds files/stan_Pot_output_AR_2023_03_12.rds")


## PM 2 - Latent Biomass (Ricker)
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}
PM_outputlist_Ricker <- lapply(stan_data_l,
                               function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                                data=x, init = init_Ricker, chains=4, iter=10000,
                                                control=list(max_treedepth=12, adapt_delta=0.99)))
names(PM_outputlist_Ricker) <- c("Pot_long","Pot_short")
launch_shinystan(PM_outputlist_Ricker$Pot_long)
launch_shinystan(PM_outputlist_Ricker$Pot_short)

saveRDS(PM_outputlist_Ricker, "./rds files/stan_Pot_output_Ricker_2023_03_15.rds")


########################################################################
## Compare out-of-sample predictions from long versus short TS
########################################################################

## Import model fit if needed
PM_outputlist_AR <- readRDS("./rds files/stan_Pot_output_AR_2023_03_12.rds")
PM_outputlist_Ricker <- readRDS("./rds files/stan_Pot_output_Ricker_2023_03_15.rds")

## Prep out-of-sample data for short versus long TS
# define oos year
Pot_oos <- data[which(data$year == "2013"),]
# extract longtrain and shorttrain max light and discharge by site
long_max <- Pot_longtrain %>% 
  summarise_at(.vars = c("Q","PAR_surface"), .funs = max)
short_max <- Pot_shorttrain %>% 
  summarise_at(.vars = c("Q","PAR_surface"), .funs = max)
# relativize Q and light by the max specific to the long and short training data sets
oos_relativize <- function(prev_max, post_dat){
  
  max.vals <- prev_max
  dat <- post_dat
  
  dat$light_rel_PAR <- dat$PAR_surface/max.vals$PAR_surface
  dat$tQ <- dat$Q/max.vals$Q
  
  dat <- dat[order(dat$date),]
  return(dat)
  
}

Pot_oos_long <- oos_relativize(long_max, Pot_oos)
Pot_oos_short <- oos_relativize(short_max, Pot_oos)

# visualize difference
ggplot(Pot_oos_long, aes(date, tQ))+
  geom_line()+
  geom_line(data = Pot_oos_short, aes(date, tQ), color="purple")

ggplot(Pot_longtrain, aes(date, tQ))+
  geom_line()+
  geom_line(data = Pot_shorttrain, aes(date, tQ), color="purple")


## source simulation models
source("Predicted_ProductivityModel_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p

## Ricker simulations using long and short posteriors
Ricker_sim_fxn <- function(x, oos_dat){
  #separate data
  output <- x
  df <- oos_dat
  
  # extract
  pars3<-extract(output, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  biomat3<-matrix(NA,length(df$GPP),length(unlist(pars3$sig_p)))
  rmsemat3<-matrix(NA,length(df$GPP),1)
  #Simulated
  for (i in 1:length(pars3$r)){
    simmat3[,i]<-PM_Ricker(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    biomat3[,i]<-PM_Ricker_B(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],pars3$sig_o[i],df)
    rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  }
  
  l <- list(simmat3, rmsemat3, biomat3)
  return(l)
  
}

Ricker_sim_Pot_long <- Ricker_sim_fxn(PM_outputlist_Ricker$Pot_long, Pot_oos_long)
names(Ricker_sim_Pot_long) <- c("GPP_sim","LB_sim","RMSE_sim")
Ricker_sim_Pot_short <- Ricker_sim_fxn(PM_outputlist_Ricker$Pot_short, Pot_oos_short)
names(Ricker_sim_Pot_short) <- c("GPP_sim","LB_sim","RMSE_sim")

## Save simulation
saveRDS(Ricker_sim_Pot_long, "./rds files/sim_Pot_long_Ricker_2023_03_15b.rds")
saveRDS(Ricker_sim_Pot_short, "./rds files/sim_Pot_short_Ricker_2023_03_15b.rds")

################################
## Predicted time series - compile original GPP data with simulated GPP based on median parameter estimates
################################

# import simulation if needed
Ricker_sim_Pot_long <- readRDS("./rds files/sim_Pot_long_Ricker_2023_03_15b.rds")
Ricker_sim_Pot_short <- readRDS("./rds files/sim_Pot_short_Ricker_2023_03_15b.rds")

# pair mean, lower, and upper CI of GPP with data
GPP_oos_preds_ts <- function(preds, dat){

  # For every day extract median and CI
  mean_simmat <- apply(preds[[1]], 1, function(x) mean(x))
  lower_simmat <- apply(preds[[1]], 1, function(x) quantile(x, probs = 0.025))
  upper_simmat <- apply(preds[[1]], 1, function(x) quantile(x, probs = 0.975))
  
  ## Plot simulated GPP
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, mean_simmat, lower_simmat, upper_simmat))
  colnames(df_sim) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  
  return(df_sim)
  
}

long_Ricker_simdat <- GPP_oos_preds_ts(Ricker_sim_Pot_long, Pot_oos_long)
short_Ricker_simdat <- GPP_oos_preds_ts(Ricker_sim_Pot_short, Pot_oos_short)


## Visualize including original used in manuscript
short_Ricker_simdat_orig <- readRDS("./rds files/Potomac_orig_oos_simdat.rds")

# colors
PM_long.col <- "#d95f02"
PM_short.col <- "#7570b3"#"red"

## Compare short simulations to check that they are the same (yes)
ggplot(short_Ricker_simdat, aes(Date, GPP))+
  geom_point(size=1.5, color="black")+
  #long training
  geom_line(aes(Date, sim_GPP), color=PM_long.col, size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_long.col, alpha=0.2, show.legend = FALSE)+
  geom_line(data=short_Ricker_simdat_orig, aes(Date, sim_GPP), color="blue", size=1.2)+
  geom_ribbon(data=short_Ricker_simdat_orig, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill="blue", alpha=0.2, show.legend = FALSE)

## Compare short to long simulations
ggplot(long_Ricker_simdat, aes(Date, GPP))+
  geom_point(size=1.5, color="black")+
  #long training
  geom_line(aes(Date, sim_GPP), color=PM_long.col, size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_long.col, alpha=0.2, show.legend = FALSE)+
  #short training
  geom_line(data=short_Ricker_simdat, aes(Date, sim_GPP), color=PM_short.col, size=1.2)+
  geom_ribbon(data=short_Ricker_simdat, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_short.col, alpha=0.2, show.legend = FALSE)+
  #previous short training
  #geom_line(data=short_Ricker_simdat_orig, aes(Date, sim_GPP), color="blue", size=1.2)+
  #geom_ribbon(data=short_Ricker_simdat_orig, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
  #            fill="blue", alpha=0.2, show.legend = FALSE)+
  #theme
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.text.x = element_text(angle=25, hjust = 1),
        axis.title.y = element_text(size = 14),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))


## Calculate metrics - coverage, NRMSE (accuracy), MRE (bias)
calc_gof_metrics <- function(x, mod){
  
  ts_sub <- x
  
  ##RMSE
  rmse <- sqrt(sum((ts_sub$sim_GPP-ts_sub$GPP)^2)/length(ts_sub$GPP))
  
  ##NRMSE
  nrmse <- rmse/(max(ts_sub$GPP)-min(ts_sub$GPP))
  
  ## RRMSE
  rrmse_pre <- sum(((ts_sub$sim_GPP-ts_sub$GPP)/ts_sub$GPP)^2)
  rrmse <- sqrt(rrmse_pre*(100/length(ts_sub$GPP)))
  
  ## MRE
  mre_pre <- sum((ts_sub$sim_GPP-ts_sub$GPP)/ts_sub$GPP)
  mre <- mre_pre*(100/length(ts_sub$GPP))
  
  ## Coverage
  ts_sub$c_yn <- ifelse(ts_sub$GPP >= ts_sub$sim_GPP_lower & ts_sub$GPP <= ts_sub$sim_GPP_upper,
                        yes=1, no=0)
  cov_pct <- (sum(ts_sub$c_yn)/length(ts_sub$c_yn))*100
  
  ## Compile
  metrics <- as.data.frame(cbind(rmse, nrmse, rrmse, mre, cov_pct,mod))
  return(metrics)
  
}

calc_gof_metrics(long_Ricker_simdat, "LB-TS long")
calc_gof_metrics(short_Ricker_simdat, "LB-TS short")
calc_gof_metrics(short_Ricker_simdat_orig, "LB-TS short orig")

## Compare 1:1 line predictions
plot_grid(
  ggplot(long_Ricker_simdat, aes(GPP, sim_GPP))+
    geom_point(color = PM_long.col)+
    scale_x_continuous(limits = c(0,15))+
    scale_y_continuous(limits = c(0,15))+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw()+
    labs(title = "Long training set (n = 4 years)"),
  ggplot(short_Ricker_simdat, aes(GPP, sim_GPP))+
    geom_point(color = PM_short.col)+
    scale_x_continuous(limits = c(0,15))+
    scale_y_continuous(limits = c(0,15))+
    geom_abline(slope = 1, intercept = 0)+
    theme_bw()+
    labs(title = "Short training set (n = 1 year)"),
  ncol = 2)
  

## fitting the model over more years doesn't substantially change the dynamics
## but different times of the year are predicted better by short versus long

################################
## Compare parameter estimates
################################

par_LBTS <- lapply(PM_outputlist_Ricker, function(x) extract(x, c("r","lambda","s","c","sig_p","sig_o")))
mean_LBTS <- lapply(par_LBTS, function(x) lapply(x, function(y) mean(y)))
sd_LBTS <- lapply(par_LBTS, function(x) lapply(x, function(y) sd(y)))

# Parameter posteriors
par_post_df <- ldply(par_LBTS, data.frame)
par_post_df2 <- gather(par_post_df, parameters, value, r:sig_o)


ggplot(par_post_df2, aes(value, fill = .id))+
  geom_density(alpha=0.2)+
  facet_wrap(~parameters, scales = "free")+
  labs(x="Value",y="Density")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=14),
        axis.text.x = element_text(size=12, angle = 45, hjust = 1),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=14), title = element_text(size=14))




