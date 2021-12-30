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
## Fit data to see parameters estimates
########################################################

## Stan data prep
rstan_options(auto_write=TRUE)
options(mc.cores=6)#parallel::detectCores())

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
#extract
p_PAR<- lapply(Ricker_output_PAR, function(x) extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))
p_PPFD<- lapply(Ricker_output_PPFD, function(x) extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))
#median
med_PAR <- lapply(p_PAR, function(x) lapply(x, function(y) median(y)))
med_PPFD <- lapply(p_PPFD, function(x) lapply(x, function(y) median(y)))

## Ricker
PM_Ricker <- function(r, lambda, s, c, sig_p, sig_o, df, light_version) {
  
  ## Data
  Ndays<-length(df$GPP)
  GPP <- df$GPP
  GPP_sd <- df$GPP_sd
  light <- light_version
  tQ <- df$tQ # discharge standardized to max value
  
  ## Vectors for model output of P, B, pred_GPP
  P <- numeric(Ndays)
  P[1] <- 1
  for(i in 2:length(tQ)){
    P[i] = exp(-exp(s*100*(tQ[i] - c)))
  }
  
  B<-numeric(Ndays)
  B[1] <- log(GPP[1]/light[1])
  pred_GPP<-numeric(Ndays)
  pred_GPP[1] <- light[1]*exp(B[1])
  
  ## Process Model
  for (j in 2:Ndays){
    B[j] <- MCMCglmm::rtnorm(1, mean = (B[j-1] + r + lambda*exp(B[j-1]))*P[j], sd = sig_p, upper = 5)
  }
  
  for (i in 2:Ndays){
    pred_GPP[i] <- MCMCglmm::rtnorm(1, mean = light[i]*exp(B[i]), sd = sig_o, lower=0.01)
  }
  
  return(pred_GPP)
}

## simulate GPP
Pot_GPP_PAR <- PM_Ricker(r = med_PAR$Potomac$r,
                         lambda = med_PAR$Potomac$lambda,
                         s = med_PAR$Potomac$s,
                         c = med_PAR$Potomac$c, 
                         sig_p = med_PAR$Potomac$sig_p,
                         sig_o = med_PAR$Potomac$sig_o, 
                         df = ex, light_version = ex$light_rel_PAR)
Pot_GPP_PPFD <- PM_Ricker(r = med_PPFD$Potomac$r,
                         lambda = med_PPFD$Potomac$lambda,
                         s = med_PPFD$Potomac$s,
                         c = med_PPFD$Potomac$c, 
                         sig_p = med_PPFD$Potomac$sig_p,
                         sig_o = med_PPFD$Potomac$sig_o, 
                         df = ex, light_version = ex$light_rel_PPFD)
Paint_GPP_PAR <- PM_Ricker(r = med_PAR$`Paint Branch`$r,
                         lambda = med_PAR$`Paint Branch`$lambda,
                         s = med_PAR$`Paint Branch`$s,
                         c = med_PAR$`Paint Branch`$c, 
                         sig_p = med_PAR$`Paint Branch`$sig_p,
                         sig_o = med_PAR$`Paint Branch`$sig_o, 
                         df = ex2, light_version = ex2$light_rel_PAR)
Paint_GPP_PPFD <- PM_Ricker(r = med_PPFD$`Paint Branch`$r,
                          lambda = med_PPFD$`Paint Branch`$lambda,
                          s = med_PPFD$`Paint Branch`$s,
                          c = med_PPFD$`Paint Branch`$c, 
                          sig_p = med_PPFD$`Paint Branch`$sig_p,
                          sig_o = med_PPFD$`Paint Branch`$sig_o, 
                          df = ex2, light_version = ex2$light_rel_PPFD)


## Plot comparison (better agreement)
plot(ex$GPP, pch=19)
lines(Pot_GPP_PAR, col="blue")
lines(Pot_GPP_PPFD, col="red")


## Plot comparison (worse agreement)
plot(ex2$GPP, pch=19)
lines(Paint_GPP_PAR, col="blue")
lines(Paint_GPP_PPFD, col="red") #better

#############################################################################
## Replace GPP with predicted with known params and fit model again
#############################################################################





#################################################################################
## Compare parameters from simulation to posterior distributions
################################################################################




samp_param <- get_sampler_params(Ricker_output_PAR$Potomac, inc_warmup = FALSE)





















