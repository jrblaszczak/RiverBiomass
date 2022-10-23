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
## (1) Source & visualize data
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

## Visualize light_rel differences - using light_rel_PAR
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
## (2) Stan data prep for initial fit
########################################################
stan_data_compile_PAR <- function(x){
  data <- list(Ndays=length(x$GPP), light = x$light_rel_PAR, GPP = x$GPP,
               GPP_sd = x$GPP_sd, tQ = x$tQ)
  return(data)
}

stan_data_PAR <- lapply(list("Potomac" = ex, "Paint Branch" = ex2), function(x) stan_data_compile_PAR(x))

###########################
## (3) Fit models to data
##########################
## Stan prep
rstan_options(auto_write=TRUE)
options(mc.cores= 8) #parallel::detectCores())

## S-TS model
STS_output_PAR <- lapply(stan_data_PAR,
                        function(x) stan("Stan_ProductivityModel1_Autoregressive_obserr.stan",
                                         data=x,chains=4,iter=5000,
                                         control=list(max_treedepth=12, adapt_delta = 0.95)))

## LB-TS model
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

LBTS_output_PAR <- lapply(stan_data_PAR,
                            function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2.stan",
                                             data=x,chains=4,iter=5000,init = init_Ricker,
                                             control=list(max_treedepth=12, adapt_delta = 0.95)))

## Save initial model fit
init_model_output <- list(STS_output_PAR, LBTS_output_PAR)
saveRDS(init_model_output, "./rds files/Param_Rec_Test_init_output_2022_07_10.rds")
rm(init_model_output); rm(STS_output_PAR); rm(LBTS_output_PAR)

################################################################
## (4) Extract parameter estimates from initial model fit
###############################################################

## Extract parameter estimates from simulation
init_model_output <- readRDS("./rds files/Param_Rec_Test_init_output_2022_07_10.rds")
STS_output_PAR <- init_model_output[[1]]
LBTS_output_PAR <- init_model_output[[2]]

#extract
p_STS <- lapply(STS_output_PAR, function(x) extract(x, c("phi","alpha","beta","sig_p","sig_o")))
p_LBTS<- lapply(LBTS_output_PAR, function(x) extract(x, c("r","lambda","s","c","sig_p","sig_o")))

#mean and sd
mean_STS <- lapply(p_STS, function(x) lapply(x, function(y) mean(y)))
sd_STS <- lapply(p_STS, function(x) lapply(x, function(y) sd(y)))
mean_LBTS <- lapply(p_LBTS, function(x) lapply(x, function(y) mean(y)))
sd_LBTS <- lapply(p_LBTS, function(x) lapply(x, function(y) sd(y)))

##############################################################
## (5) Simulate GPP ts using extracted parameter estimates
##############################################################

## Bring in simulation code
source("Predicted_ProductivityModel_Autoregressive.R") 
source("Predicted_ProductivityModel_Ricker.R")

## STS - simulate GPP again for recovery test data set
Pot_STS_GPP <- PM_AR(phi=mean_STS$Potomac$phi,
                    alpha=mean_STS$Potomac$alpha,
                    beta=mean_STS$Potomac$beta,
                    sig_p=mean_STS$Potomac$sig_p,
                    sig_o=mean_STS$Potomac$sig_o, df=ex)
Paint_STS_GPP <- PM_AR(phi=mean_STS$`Paint Branch`$phi,
                      alpha=mean_STS$`Paint Branch`$alpha,
                      beta=mean_STS$`Paint Branch`$beta,
                      sig_p=mean_STS$`Paint Branch`$sig_p,
                      sig_o=mean_STS$`Paint Branch`$sig_o, df=ex2)
plot(1:length(Pot_STS_GPP), Pot_STS_GPP, type="l")
points(ex$GPP) #good
plot(1:length(Paint_STS_GPP), Paint_STS_GPP, type="l")
points(ex2$GPP) #divergent fall peak in predicted (line) GPP

## LBTS - simulate GPP again for recovery test data set
Pot_LBTS_GPP <- PM_Ricker(r = mean_LBTS$Potomac$r,
                            lambda = mean_LBTS$Potomac$lambda,
                            s = mean_LBTS$Potomac$s,
                            c = mean_LBTS$Potomac$c, 
                            sig_p = mean_LBTS$Potomac$sig_p,
                            sig_o = mean_LBTS$Potomac$sig_o, df = ex)
Paint_LBTS_GPP <- PM_Ricker(r = mean_LBTS$`Paint Branch`$r,
                              lambda = mean_LBTS$`Paint Branch`$lambda,
                              s = mean_LBTS$`Paint Branch`$s,
                              c = mean_LBTS$`Paint Branch`$c, 
                              sig_p = mean_LBTS$`Paint Branch`$sig_p,
                              sig_o = mean_LBTS$`Paint Branch`$sig_o, df = ex2)
plot(1:length(Pot_LBTS_GPP), Pot_LBTS_GPP, type="l")
points(ex$GPP) #good
plot(1:length(Paint_LBTS_GPP), Paint_LBTS_GPP, type="l")
points(ex2$GPP) #much better than S-TS prediction, just delayed spring peak


#############################################################################
## (6) Replace original GPP with predicted GPP in stan data list
#############################################################################

## Compile for Stan again
# STS
stan_simPot_STS <- list(Ndays=length(Pot_STS_GPP), light=ex$light_rel_PAR, GPP = Pot_STS_GPP,
                        prior_sig_o_mean = mean_STS$Potomac$sig_o,
                        prior_sig_o_sd = sd_STS$Potomac$sig_o, tQ = ex$tQ)
stan_simPaint_STS <- list(Ndays=length(Paint_STS_GPP), light=ex2$light_rel_PAR, GPP = Paint_STS_GPP,
                          prior_sig_o_mean = mean_STS$`Paint Branch`$sig_o,
                          prior_sig_o_sd = sd_STS$`Paint Branch`$sig_o, tQ = ex2$tQ)
stan_STS_sim_list <- list("Potomac_STS"=stan_simPot_STS,
                         "PaintBranch_STS"=stan_simPaint_STS)
#LBTS
stan_simPot_LBTS <- list(Ndays=length(Pot_LBTS_GPP), light=ex$light_rel_PAR, GPP = Pot_LBTS_GPP,
                        prior_sig_o_mean = mean_LBTS$Potomac$sig_o,
                        prior_sig_o_sd = sd_LBTS$Potomac$sig_o, tQ = ex$tQ)
stan_simPaint_LBTS <- list(Ndays=length(Paint_LBTS_GPP), light=ex2$light_rel_PAR, GPP = Paint_LBTS_GPP,
                          prior_sig_o_mean = mean_LBTS$`Paint Branch`$sig_o,
                          prior_sig_o_sd = sd_LBTS$`Paint Branch`$sig_o, tQ = ex2$tQ)
stan_LBTS_sim_list <- list("Potomac_LBTS"=stan_simPot_LBTS,
                      "PaintBranch_LBTS"=stan_simPaint_LBTS)

###########################################################################
## (7) Fit models to simulated data from initial parameter estimates
###########################################################################
rstan_options(auto_write=TRUE)
options(mc.cores=8)#parallel::detectCores())
#STS
recov_STS_output <- lapply(stan_STS_sim_list,
                       function(x) stan("Stan_ProductivityModel1_Autoregressive_obserr_simulation.stan",
                                        data=x,chains=4,iter=5000,
                                        control=list(max_treedepth=12, adapt_delta = 0.95)))

#LBTS
init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}
recov_LBTS_output <- lapply(stan_LBTS_sim_list,
                            function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2_simulation.stan",
                                             data=x,chains=4,iter=5000,init = init_Ricker,
                                             control=list(max_treedepth=12)))

## Save parameter recovery model fit
recovery_model_output <- list(recov_STS_output, recov_LBTS_output)
saveRDS(recovery_model_output, "./rds files/Param_Rec_Test_recovery_output_2022_07_10.rds")
rm(recovery_model_output); rm(recov_STS_output); rm(recov_LBTS_output)


#################################################################################
## (8) Compare parameters from simulation to posterior distributions
################################################################################

recovery_model_output <- readRDS("./rds files/Param_Rec_Test_recovery_output_2022_07_10.rds")
recov_STS_output <- recovery_model_output[[1]]
recov_LBTS_output <- recovery_model_output[[2]]

## STS - vis_recovery
STS_sim_recp <- lapply(recov_STS_output, function(x) ldply(extract(x, c("phi","alpha","beta","sig_p","sig_o")),data.frame))
# need STS_sim_recp, and mean_STS from part 4
STS_vis_recovery <- function(final_distributions, used_parameters, plot.title){
  
  ## distributions of most recent fit
  recpars <- final_distributions
  colnames(recpars) <- c("parameter","value")
  
  ## params used to simulate
  orig_meanpars <- ldply(used_parameters, data.frame)
  colnames(orig_meanpars) <- c("parameter","value")
  
  
  ggplot(recpars, aes(value))+
    geom_density(fill="red", alpha=0.2)+
    facet_wrap(~parameter, scales = "free")+
    geom_vline(data=orig_meanpars,aes(xintercept = value))+
    labs(x="Value",y="Density",title=plot.title)+
    theme(legend.position = "none",
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=14),
          axis.text.x = element_text(size=12, angle = 45, hjust = 0.5),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14), title = element_text(size=14))
}

STS_vis_recovery(STS_sim_recp$Potomac_STS, mean_STS$Potomac, "S-TS Potomac Parameter Recovery")
STS_vis_recovery(STS_sim_recp$PaintBranch_STS, mean_STS$`Paint Branch`, "S-TS Paint Branch Parameter Recovery")


## LBTS - vis_recovery
LBTS_sim_recp <- lapply(recov_LBTS_output, function(x) ldply(extract(x, c("r","lambda","s","c","sig_p","sig_o")),data.frame))
# need LBTS_sim_recp, and mean_LBTS from part 4
LBTS_vis_recovery <- function(final_distributions, used_parameters, plot.title){
  
  ## distributions of most recent fit
  recpars <- final_distributions
  colnames(recpars) <- c("parameter","value")
  
  ## params used to simulate
  orig_meanpars <- ldply(used_parameters, data.frame)
  colnames(orig_meanpars) <- c("parameter","value")
  
  
  ggplot(recpars, aes(value))+
    geom_density(fill="red", alpha=0.2)+
    facet_wrap(~parameter, scales = "free")+
    geom_vline(data=orig_meanpars,aes(xintercept = value))+
    labs(x="Value",y="Density",title=plot.title)+
    theme(legend.position = "none",
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=14),
          axis.text.x = element_text(size=12, angle = 45, hjust = 0.5),
          axis.text.y = element_text(size=12),
          axis.title = element_text(size=14), title = element_text(size=14))
}

LBTS_vis_recovery(LBTS_sim_recp$Potomac_LBTS, mean_LBTS$Potomac, "LB-TS Potomac Parameter Recovery")
LBTS_vis_recovery(LBTS_sim_recp$PaintBranch_LBTS, mean_LBTS$`Paint Branch`, "LB-TS Paint Branch Parameter Recovery")





