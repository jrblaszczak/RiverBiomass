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

saveRDS(sim_Ricker_output, "./rds files/sim_Ricker_output_2021_01_14.rds")


################################################
## Extract parameter estimates
################################################

## Extract parameter estimates from simulation
sim_Ricker_output <- readRDS("./rds files/sim_Ricker_output_2021_01_14.rds")
Ricker_output_PAR <- sim_Ricker_output[[1]]
Ricker_output_PPFD <- sim_Ricker_output[[2]]

#extract
p_PAR<- lapply(Ricker_output_PAR, function(x) extract(x, c("r","lambda","s","c","sig_p","sig_o")))
p_PPFD<- lapply(Ricker_output_PPFD, function(x) extract(x, c("r","lambda","s","c","sig_p","sig_o")))
#mean and sd
mean_PAR <- lapply(p_PAR, function(x) lapply(x, function(y) mean(y)))
sd_PAR <- lapply(p_PAR, function(x) lapply(x, function(y) sd(y)))
mean_PPFD <- lapply(p_PPFD, function(x) lapply(x, function(y) mean(y)))
sd_PPFD <- lapply(p_PPFD, function(x) lapply(x, function(y) sd(y)))

############################################
## Simulate final GPP ts using extracted parameter estimates
############################################

## Bring in simulation code
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p

## simulate GPP again for final data set
Pot_GPP_PAR <- PM_Ricker_lv(r = mean_PAR$Potomac$r,
                         lambda = mean_PAR$Potomac$lambda,
                         s = mean_PAR$Potomac$s,
                         c = mean_PAR$Potomac$c, 
                         sig_p = mean_PAR$Potomac$sig_p,
                         sig_o = mean_PAR$Potomac$sig_o, 
                         df = ex, light_version = ex$light_rel_PAR)
Pot_GPP_PPFD <- PM_Ricker_lv(r = mean_PPFD$Potomac$r,
                          lambda = mean_PPFD$Potomac$lambda,
                          s = mean_PPFD$Potomac$s,
                          c = mean_PPFD$Potomac$c, 
                          sig_p = mean_PPFD$Potomac$sig_p,
                          sig_o = mean_PPFD$Potomac$sig_o, 
                          df = ex, light_version = ex$light_rel_PPFD)
Paint_GPP_PAR <- PM_Ricker_lv(r = mean_PAR$`Paint Branch`$r,
                           lambda = mean_PAR$`Paint Branch`$lambda,
                           s = mean_PAR$`Paint Branch`$s,
                           c = mean_PAR$`Paint Branch`$c, 
                           sig_p = mean_PAR$`Paint Branch`$sig_p,
                           sig_o = mean_PAR$`Paint Branch`$sig_o, 
                           df = ex2, light_version = ex2$light_rel_PAR)
Paint_GPP_PPFD <- PM_Ricker_lv(r = mean_PPFD$`Paint Branch`$r,
                            lambda = mean_PPFD$`Paint Branch`$lambda,
                            s = mean_PPFD$`Paint Branch`$s,
                            c = mean_PPFD$`Paint Branch`$c, 
                            sig_p = mean_PPFD$`Paint Branch`$sig_p,
                            sig_o = mean_PPFD$`Paint Branch`$sig_o, 
                            df = ex2, light_version = ex2$light_rel_PPFD)

plot(1:352, Pot_GPP_PAR, type="l")
lines(1:352, Pot_GPP_PPFD, type="l", col="blue")

plot(1:312, Paint_GPP_PAR, type="l")
lines(1:312, Paint_GPP_PPFD, type="l", col="blue")


#############################################################################
## Replace original GPP with predicted GPP in data list and fit model again
#############################################################################

stan_simPot_PAR <- list(Ndays=length(Pot_GPP_PAR), light=ex$light_rel_PAR, GPP = Pot_GPP_PAR,
                        prior_sig_o_mean = mean_PAR$Potomac$sig_o,
                        prior_sig_o_sd = sd_PAR$Potomac$sig_o, tQ = ex$tQ)
stan_simPot_PPFD <- list(Ndays=length(Pot_GPP_PPFD), light=ex$light_rel_PPFD, GPP = Pot_GPP_PPFD,
                        prior_sig_o_mean = mean_PPFD$Potomac$sig_o,
                        prior_sig_o_sd = sd_PPFD$Potomac$sig_o, tQ = ex$tQ)

stan_simPaint_PAR <- list(Ndays=length(Paint_GPP_PAR), light=ex2$light_rel_PAR, GPP = Paint_GPP_PAR,
                        prior_sig_o_mean = mean_PAR$`Paint Branch`$sig_o,
                        prior_sig_o_sd = sd_PAR$`Paint Branch`$sig_o, tQ = ex2$tQ)
stan_simPaint_PPFD <- list(Ndays=length(Paint_GPP_PPFD), light=ex2$light_rel_PPFD, GPP = Paint_GPP_PPFD,
                         prior_sig_o_mean = mean_PPFD$`Paint Branch`$sig_o,
                         prior_sig_o_sd = sd_PPFD$`Paint Branch`$sig_o, tQ = ex2$tQ)

stan_sim_list <- list("Potomac_PAR"=stan_simPot_PAR,
                      "Potomac_PPFD"=stan_simPot_PPFD,
                      "PaintBranch_PAR"=stan_simPaint_PAR,
                      "PaintBranch_PPFD"=stan_simPaint_PPFD)

#Fit models again
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

init_Ricker <- function(...) {
  list(c = 0.5, s = 1.5)
}

Ricker_sim_output <- lapply(stan_sim_list,
                            function(x) stan("Stan_ProductivityModel2_Ricker_s_mod2_simulation.stan",
                                             data=x,chains=3,iter=5000,init = init_Ricker,
                                             control=list(max_treedepth=12)))

saveRDS(Ricker_sim_output, "./rds files/Ricker_recover_output.rds")


#################################################################################
## Compare parameters from simulation to posterior distributions
################################################################################

#Ricker_sim_output <- readRDS("./rds files/Ricker_recover_output.rds")

## extract parameters from fit to simulated GPP
sim_recparams <- lapply(Ricker_sim_output, function(x) ldply(extract(x, c("r","lambda","s","c","sig_p","sig_o")),data.frame))

vis_recovery <- function(final_distributions, used_parameters){
  
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
    theme_bw()+
    theme()
}

vis_recovery(sim_recparams$Potomac_PAR, mean_PAR$Potomac)
vis_recovery(sim_recparams$Potomac_PPFD, mean_PPFD$Potomac)

vis_recovery(sim_recparams$PaintBranch_PPFD, mean_PPFD$`Paint Branch`)
vis_recovery(sim_recparams$PaintBranch_PAR, mean_PAR$`Paint Branch`)



test <- Ricker_sim_output$PaintBranch_PAR
launch_shinystan(Ricker_output_PPFD$`Paint Branch`)











