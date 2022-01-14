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

########################################################
## Stan data prep for initial fit
########################################################

stan_data_compile_PAR <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}
stan_data_compile_PPFD <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PPFD, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_PAR <- lapply(list("Potomac" = ex, "Paint Branch" = ex2), function(x) stan_data_compile_PAR(x))
stan_data_PPFD <- lapply(list("Potomac" = ex, "Paint Branch" = ex2), function(x) stan_data_compile_PPFD(x))


###########################
## Fit data to models
##########################
## Stan prep
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

Ricker_output_PAR <- lapply(stan_data_PAR,
                            function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                             data=x,chains=3,iter=5000,init = init_Ricker,
                                             control=list(max_treedepth=12)))
Ricker_output_PPFD <- lapply(stan_data_PPFD,
                             function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                              data=x,chains=3,iter=5000,init = init_Ricker,
                                              control=list(max_treedepth=12)))
sim_Ricker_output <- list(Ricker_output_PAR, Ricker_output_PPFD)

saveRDS(sim_Ricker_output, "./rds files/sim_Ricker_output_2021_12_28.rds")


################################################
## Simulate data using deterministic function
################################################

## Extract parameter estimates from simulation
#sim_Ricker_output <- readRDS("./rds files/sim_Ricker_output_2021_12_28.rds")
Ricker_output_PAR <- sim_Ricker_output[[1]]
Ricker_output_PPFD <- sim_Ricker_output[[2]]

## Bring in simulation code
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p

## Pair up data
names(stan_data_PAR); names(Ricker_output_PAR)
names(stan_data_PPFD); names(Ricker_output_PPFD)
PAR_list <- Map(c, Ricker_output_PAR, stan_data_PAR)
PPFD_list <- Map(c, Ricker_output_PPFD, stan_data_PPFD)

Ricker_sim_fxn <- function(x){
  
  # separate data
  output <- x[[1]]
  df <- x

  # extract
  pars<-extract(output, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o"))
  simmat<-matrix(NA,length(df$GPP),length(unlist(pars$sig_p)))
  biomat<-matrix(NA,length(df$GPP),length(unlist(pars$sig_p)))
  rmsemat<-matrix(NA,length(df$GPP),1)
  # simulate
  for (i in 1:length(pars$r)){
    simmat[,i]<-PM_Ricker(pars$r[i],pars$lambda[i],pars$s[i],pars$c[i],pars$sig_p[i],pars$sig_o[i],df)
    }
  
  return(simmat)
  
}

Ricker_sim_PAR <- lapply(PAR_list, function(x) Ricker_sim_fxn(x))
Ricker_sim_PPFD <- lapply(PPFD_list, function(x) Ricker_sim_fxn(x))

# For every day extract mean and sd of GPP
mean_sim_PAR <- ldply(lapply(Ricker_sim_PAR, function(z) apply(z, 1, function(x) mean(x))), data.frame)
sd_sim_PAR <- ldply(lapply(Ricker_sim_PAR, function(z) apply(z, 1, function(x) sd(x))), data.frame)
mean_sim_PPFD <- ldply(lapply(Ricker_sim_PPFD, function(z) apply(z, 1, function(x) mean(x))), data.frame)
sd_sim_PPFD <- ldply(lapply(Ricker_sim_PPFD, function(z) apply(z, 1, function(x) sd(x))), data.frame)
#rename PAR colnames
sim_PAR <- cbind(mean_sim_PAR, sd_sim_PAR[,2])
colnames(sim_PAR) <- c("SiteID", "mean_GPP_PAR","sd_GPP_PAR")
#rename PPFD colnames
sim_PPFD <- cbind(mean_sim_PPFD, sd_sim_PPFD[,2])
colnames(sim_PPFD) <- c("SiteID", "mean_GPP_PPFD","sd_GPP_PPFD")

## Visualize
sim_PAR$order <- as.numeric(rownames(sim_PAR)); sim_PPFD$order <- as.numeric(rownames(sim_PPFD))
ggplot(sim_PAR, aes(order, mean_GPP_PAR))+
  geom_line()+facet_wrap(~SiteID, nrow=1, scales="free_x")
ggplot(sim_PPFD, aes(order, mean_GPP_PPFD))+
  geom_line()+facet_wrap(~SiteID, nrow=1, scales="free_x")


#############################################################################
## Replace original GPP with predicted mean and sd GPP in data list and fit model again
#############################################################################
#create copies of data lists
stan_data_PAR_v2 <- stan_data_PAR
stan_data_PPFD_v2 <- stan_data_PPFD
#replace GPP for PAR
stan_data_PAR_v2$Potomac$GPP <- sim_PAR[which(sim_PAR$SiteID == "Potomac"),]$mean_GPP_PAR
stan_data_PAR_v2$`Paint Branch`$GPP <- sim_PAR[which(sim_PAR$SiteID == "Paint Branch"),]$mean_GPP_PAR
stan_data_PAR_v2$Potomac$GPP_sd <- sim_PAR[which(sim_PAR$SiteID == "Potomac"),]$sd_GPP_PAR
stan_data_PAR_v2$`Paint Branch`$GPP_sd <- sim_PAR[which(sim_PAR$SiteID == "Paint Branch"),]$sd_GPP_PAR
#replace GPP for PPFD
stan_data_PPFD_v2$Potomac$GPP <- sim_PPFD[which(sim_PPFD$SiteID == "Potomac"),]$mean_GPP_PPFD
stan_data_PPFD_v2$`Paint Branch`$GPP <- sim_PPFD[which(sim_PPFD$SiteID == "Paint Branch"),]$mean_GPP_PPFD
stan_data_PPFD_v2$Potomac$GPP_sd <- sim_PPFD[which(sim_PPFD$SiteID == "Potomac"),]$sd_GPP_PPFD
stan_data_PPFD_v2$`Paint Branch`$GPP_sd <- sim_PPFD[which(sim_PPFD$SiteID == "Paint Branch"),]$sd_GPP_PPFD


#Fit models again
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

Ricker_output_PAR_v2 <- lapply(stan_data_PAR_v2,
                            function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                             data=x,chains=3,iter=5000,init = init_Ricker,
                                             control=list(max_treedepth=12)))
Ricker_output_PPFD_v2 <- lapply(stan_data_PPFD_v2,
                             function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                              data=x,chains=3,iter=5000,init = init_Ricker,
                                              control=list(max_treedepth=12)))
sim_Ricker_output_v2 <- list(Ricker_output_PAR_v2, Ricker_output_PPFD_v2)

saveRDS(sim_Ricker_output_v2, "./rds files/sim_Ricker_output_v2_2022_01_03.rds")

#################################################################################
## Compare parameters from simulation to posterior distributions
################################################################################




samp_param <- get_sampler_params(Ricker_output_PAR$Potomac, inc_warmup = FALSE)





















