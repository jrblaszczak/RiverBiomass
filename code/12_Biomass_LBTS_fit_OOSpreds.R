##==============================================================================
## Script for new Figure 2 - rmax distribution, LB-TS fit, LB, & out-of-sample LB-TS prediction and 1:1 plots
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2","ggExtra","patchwork","grid","gridExtra"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")
## Source functions to extract and summarize parameter estimates
source("StanParameterExtraction_Source.R")

## Import Stan model fit
stan_model_output_LBTS <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")

##############################################
## A) rmax distribution comparisons
##############################################
## Extract parameters
##############################################
par_LBTS <- lapply(stan_model_output_LBTS, function(x) rstan::extract(x, c("r","lambda","s","c","B","P","pred_GPP","sig_p","sig_o")))

r_K_func <- function(x) {
  rx <- x$r
  rm <- as.data.frame(cbind(rx[1:2500],rx[2501:5000],rx[5001:7500]))
  rm <- rm %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  lx <- x$lambda
  lambda <- as.data.frame(cbind(lx[1:2500],lx[2501:5000],lx[5001:7500]))
  lambda <- lambda %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  new <- as.data.frame(cbind(rm$Avg, lambda$Avg))
  colnames(new) <- c("r","lambda")
  new$K <- (-1*new$r)/new$lambda
  
  long <- gather(new)
  
  return(long)
}


rK_ws <- ldply(lapply(par_LBTS, function(x) r_K_func(x)), data.frame)

rk_formatting <- function(x){
  x$short_name <- revalue(as.character(x$.id), replace = c("nwis_02336526"="Proctor Creek, GA",
                                                           "nwis_01649190"="Paint Branch, MD",
                                                           "nwis_07191222"="Beaty Creek, OK",
                                                           "nwis_01608500"="S. Br. Potomac River, WV",
                                                           "nwis_11044000"="Santa Margarita River, CA",
                                                           "nwis_08447300"="Pecos River, TX"))
  x$short_name_SO <- revalue(as.character(x$.id), replace = c("nwis_02336526"="Proctor Creek, GA (2)",
                                                              "nwis_01649190"="Paint Branch, MD (2)",
                                                              "nwis_07191222"="Beaty Creek, OK (3)",
                                                              "nwis_01608500"="S. Br. Potomac River, WV (5)",
                                                              "nwis_11044000"="Santa Margarita River, CA (6)",
                                                              "nwis_08447300"="Pecos River, TX (7)"))
  x$short_name_SO <- factor(x$short_name_SO, levels= c("Proctor Creek, GA (2)",
                                                       "Paint Branch, MD (2)",
                                                       "Beaty Creek, OK (3)",
                                                       "S. Br. Potomac River, WV (5)", 
                                                       "Santa Margarita River, CA (6)",
                                                       "Pecos River, TX (7)"))
  #x$key <- factor(x$key, levels = c("r","lambda","K"))
  return(x)
}

rK_vals <- rk_formatting(rK_ws)  


## r plot
r_vals <- rK_vals[which(rK_vals$key == "r"),]
r_vals$max_inc <- exp(r_vals$value) - 1


r_max_colors <- hcl.colors(n=6, palette = "Spectral")

r_plot <- ggplot(r_vals, aes(max_inc, fill = short_name_SO))+
  geom_density(color="black", alpha=0.75)+
  #colors scaled by rmax value
  scale_fill_manual(values = c("Pecos River, TX (7)" = r_max_colors[1],
                               "Paint Branch, MD (2)" = r_max_colors[2],
                               "Santa Margarita River, CA (6)" = r_max_colors[3],
                               "Proctor Creek, GA (2)" = r_max_colors[4],
                               "S. Br. Potomac River, WV (5)" = r_max_colors[5],
                               "Beaty Creek, OK (3)" = r_max_colors[6]))+
  scale_x_continuous(limits = c(0,0.6), labels = scales::percent)+
  labs(x = "Daily Maximum % Increase in Biomass",
       y = "Density")+
  theme_bw()
r_plot



## K plot
K_vals <- rK_vals[which(rK_vals$key == "K"),]

K_plot <- ggplot(K_vals, aes(value, fill = short_name_SO))+
  geom_density(color="black", alpha=0.5)+
  #scale_x_continuous(limits = c(0,0.6), labels = scales::percent)+
  labs(x = "Latent Biomass Carry Capacity",
       y = "Density")+
  theme_bw()



## Inspect r vs K covariance
r_K_func2 <- function(x) {
  rx <- x$r
  rm <- as.data.frame(cbind(rx[1:2500],rx[2501:5000],rx[5001:7500]))
  rm <- rm %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  lx <- x$lambda
  lambda <- as.data.frame(cbind(lx[1:2500],lx[2501:5000],lx[5001:7500]))
  lambda <- lambda %>%
    rowwise() %>% mutate(Avg=mean(c(V1, V2, V3))) 
  
  new <- as.data.frame(cbind(rm$Avg, lambda$Avg))
  colnames(new) <- c("r","lambda")
  new$K <- (-1*new$r)/new$lambda

  return(new)
}


rK_covar_df <- ldply(lapply(par_LBTS, function(x) r_K_func2(x)), data.frame)
rK_covar_df <- rk_formatting(rK_covar_df)
rK_covar_df$short_name <- factor(rK_covar_df$short_name, levels= site_order_list)

ggplot(rK_covar_df, aes(r, K))+
  geom_point()+
  facet_wrap(~short_name, ncol = 2, scales = "free")+
  theme_bw()+
  theme(legend.position = "none",
          panel.background = element_rect(color = "black", fill=NA, size=1),
          axis.text = element_text(size=12),
          axis.title = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
          strip.background = element_rect(fill="white", color="black"),
          strip.text = element_text(size=15))

#investigate correlation between r and K
rK_covar_list <- split(rK_covar_df, rK_covar_df$short_name)

lapply(rK_covar_list, function(x) cor.test(x = x$r, y = x$K))

##############################################
## B) WS model fit
##############################################
## extract estimates
LBTS_medpar <- ldply(lapply(stan_model_output_LBTS,
                           function(x) LBTS_extract_medians(rstan::extract(x,c("B","P","pred_GPP")))),
                    data.frame)

## Create data frame
dat <- ldply(df, data.frame)

## GPP fit
lbts_GPP <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, LBTS_medpar$pred_GPP, LBTS_medpar$pred_GPP_Q.025, LBTS_medpar$pred_GPP_Q.975))
colnames(lbts_GPP) <- c("site_name","Date","GPP","pred_GPP","pred_GPP_lower","pred_GPP_upper")
lbts_GPP$Date <- as.POSIXct(as.character(lbts_GPP$Date), format="%Y-%m-%d")
lbts_GPP[,3:6] <- apply(lbts_GPP[,3:6],2,function(x) as.numeric(as.character(x)))
## Arrange rivers by river order
lbts_GPP <- left_join(lbts_GPP, site_info[,c("site_name","short_name")])
lbts_GPP$short_name <- factor(lbts_GPP$short_name, levels=site_order_list)

## Latent Biomass fit
df_modB3 <- as.data.frame(cbind(dat$site_name, as.character(dat$date), PM3_medpar$B, PM3_medpar$B_Q.025, PM3_medpar$B_Q.975))
colnames(df_modB3) <- c("site_name","Date","B","B_lower","B_upper")
df_modB3$Date <- as.POSIXct(as.character(df_modB3$Date), format="%Y-%m-%d")
df_modB3[,3:5] <- apply(df_modB3[,3:5],2,function(x) as.numeric(as.character(x)))
## Arrange rivers by river order
df_modB3 <- left_join(df_modB3, site_info[,c("site_name","short_name")])
df_modB3$short_name <- factor(df_modB3$short_name, levels=site_order_list)

## Subset to only Potomac River
Pot_GPP <- lbts_GPP[which(lbts_GPP$short_name == "S. Br. Potomac River, WV"),]
Pot_LB <- df_modB3[which(df_modB3$short_name == "S. Br. Potomac River, WV"),]

plot_grid(
## Plot GPP pred
ggplot(Pot_GPP, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, pred_GPP), color=PM_Ricker.col, size=0.75)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="Within-Sample LB-TS Model Fit")+
  geom_ribbon(aes(ymin=pred_GPP_lower,ymax=pred_GPP_upper),
              fill=PM_Ricker.col, alpha=0.5, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15)),

## Plot latent biomass predictions
ggplot(Pot_LB, aes(Date, exp(B)))+
  geom_line(size=0.75, color="chartreuse4")+
  labs(y="Latent Biomass",title="Within-Sample LB-TS Model Estimates")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.5, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15)),
ncol = 1, align = "hv")


##############################################
## B) OOS Model estimates
##############################################

source("DataSource_6rivers_oos_StreamLight.R")
source("Predicted_ProductivityModel_Ricker.R") # parameters: r, lambda, s, c, sig_p

## Generate predictions
GPP_oos_preds_ts <- function(preds, df){
  
  simmat_list <- readRDS(preds)
  simmat_list <- simmat_list[names(df)]
  
  # For every day extract median and CI
  median_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) median(x))), data.frame)
  lower_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.025))), data.frame)
  upper_simmat <- ldply(lapply(simmat_list, function(z) apply(z[[1]], 1, function(x) quantile(x, probs = 0.975))), data.frame)
  
  ## Plot simulated GPP
  dat <- ldply(df, data.frame)
  df_sim <- as.data.frame(cbind(dat$site_name, as.character(dat$date), dat$GPP, median_simmat$X..i.., lower_simmat$X..i.., upper_simmat$X..i..))
  colnames(df_sim) <- c("site_name","Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
  df_sim$Date <- as.POSIXct(as.character(df_sim$Date), format="%Y-%m-%d")
  df_sim[,3:6] <- apply(df_sim[,3:6],2,function(x) as.numeric(as.character(x)))
  
  ## Arrange rivers by river order
  df_sim <- left_join(df_sim, site_info[,c("site_name","short_name")])
  df_sim$short_name <- factor(df_sim$short_name, levels=site_order_list)
  
  return(df_sim)
  
}

LB_simdat <- GPP_oos_preds_ts("./rds files/Sim_6riv_Ricker_oos_2022_02_27.rds", dat_oos)

## Subset to Potomac
scaleFUN <- function(x) sprintf("%.1f", x)

LB_simdat_site <- LB_simdat[which(LB_simdat$short_name == "S. Br. Potomac River, WV"),]

xy.limits <- range(c(LB_simdat_site$GPP,
                     LB_simdat_site$sim_GPP))

inset_plot <- ggplot(LB_simdat_site, aes(GPP, sim_GPP))+
  geom_point(color = PM_Ricker.col, size=0.5)+
  scale_x_continuous(limits=c(0,xy.limits[2])) + 
  scale_y_continuous(limits=c(0,xy.limits[2])) +
  labs(x="GPP Data",y="Predicted\n GPP")+
  geom_abline(slope = 1, intercept = 0)+
  theme_classic()+
  theme(axis.title = element_text(size = 8),
        axis.text = element_text(size = 7))


main_plot <- ggplot(LB_simdat_site, aes(Date, GPP))+
  geom_point(size=1, color="black")+
  
  geom_line(data=LB_simdat_site, aes(Date, sim_GPP), color=PM_Ricker.col, size=1)+
  geom_ribbon(data=LB_simdat_site, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM_Ricker.col, alpha=0.4, show.legend = FALSE)+
  
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  coord_cartesian(ylim=c(0,max(LB_simdat_site$sim_GPP_upper)*2))+
  scale_y_continuous(labels=scaleFUN)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
       title = "Out-of-Sample GPP Predictions")

plot.with.inset <- ggdraw() +
  draw_plot(main_plot) +
  draw_plot(inset_plot, x = 0.58, y = 0.52, width = .3, height = .3)
plot.with.inset

comb.plot <- ggdraw()+
  draw_plot(main_plot)+
  draw_plot(inset_plot, x=0.58, y=0.52, width=0.35, height = 0.35)
comb.plot

######################################
## Together
######################################

r_max_colors <- hcl.colors(n=6, palette = "Spectral")
PM_Ricker.col <- r_max_colors[5] #to match D and E

plot_grid(
  
  plot_grid(
    ## Plot GPP estimates
    ggplot(Pot_GPP, aes(Date, GPP))+
      geom_point(size=2, color="black")+
      geom_line(aes(Date, pred_GPP), color=PM_Ricker.col, size=0.75)+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
           title="2012 Within-Sample LB-TS Model GPP Estimates")+
      geom_ribbon(aes(ymin=pred_GPP_lower,ymax=pred_GPP_upper),
                  fill=PM_Ricker.col, alpha=0.5, show.legend = FALSE)+
      theme_classic()+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle=25, hjust = 1)),
    
    ## Plot latent biomass predictions
    ggplot(Pot_LB, aes(Date, exp(B)))+
      geom_line(size=0.75, color="chartreuse4")+
      labs(y="Latent Biomass",title="2012 Within-Sample LB-TS Model Biomass Estimates")+
      geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
                  fill="chartreuse4", alpha=0.5, show.legend = FALSE)+
      theme_classic()+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle=25, hjust = 1)),
    
    ## Plot OOS predictions
    ggplot(LB_simdat_site, aes(Date, GPP))+
      geom_point(size=1, color="black")+
      geom_line(data=LB_simdat_site, aes(Date, sim_GPP), color=PM_Ricker.col, size=0.75)+
      geom_ribbon(data=LB_simdat_site, aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
                  fill=PM_Ricker.col, alpha=0.4, show.legend = FALSE)+
      theme_classic()+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(),
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(angle=25, hjust = 1))+
      coord_cartesian(ylim=c(0,max(LB_simdat_site$sim_GPP_upper)*1.75))+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),
           title = "2013 Out-of-Sample LB-TS Model GPP Predictions"),
    
    ncol = 1, align = "hv", labels = c("A","B","C")),
  
  plot_grid(
    ## r estimate distributions
    ggplot(r_vals, aes(max_inc, fill = short_name_SO))+
      geom_density(color="black", alpha=0.5)+
      scale_x_continuous(limits = c(0,0.6), labels = scales::percent)+
      labs(x = "Daily Maximum % Increase in Biomass",
           y = "Density")+
      scale_fill_manual(values = c("Pecos River, TX (7)" = r_max_colors[1],
                                   "Paint Branch, MD (2)" = r_max_colors[2],
                                   "Santa Margarita River, CA (6)" = r_max_colors[3],
                                   "Proctor Creek, GA (2)" = r_max_colors[4],
                                   "S. Br. Potomac River, WV (5)" = r_max_colors[5],
                                   "Beaty Creek, OK (3)" = r_max_colors[6]))+
      theme_classic()+
      theme(legend.position = "none",
            axis.text = element_text(size = 10),
            panel.background = element_rect(color = "black", fill=NA, size=1)),
    
    ## K estimate distributions
    ggplot(K_vals, aes(value, fill = short_name_SO))+
      geom_density(color="black", alpha=0.6)+
      labs(x = "Latent Biomass Carry Capacity (K)",
           y = "Density")+
      scale_fill_manual(values = c("Pecos River, TX (7)" = r_max_colors[1],
                                   "Paint Branch, MD (2)" = r_max_colors[2],
                                   "Santa Margarita River, CA (6)" = r_max_colors[3],
                                   "Proctor Creek, GA (2)" = r_max_colors[4],
                                   "S. Br. Potomac River, WV (5)" = r_max_colors[5],
                                   "Beaty Creek, OK (3)" = r_max_colors[6]))+
      theme_classic()+
      theme(legend.position = "none",
            axis.text = element_text(size = 10),
            panel.background = element_rect(color = "black", fill=NA, size=1)),
    
    ncol = 1, align = "hv", labels = c("D","E")),
  
  ncol = 2, rel_widths = c(1, 0.9))

## save legend
plot_grid(get_legend(r_plot+theme(legend.position = "bottom")))

## save inset plot for Potomac and add in illustrator to Panel C (OOS preds)
ggplot(LB_simdat_site, aes(GPP, sim_GPP))+
  geom_point(color = PM_Ricker.col, size=1)+
  scale_x_continuous(limits=c(0,xy.limits[2]), expand = c(0,0.5)) + 
  scale_y_continuous(limits=c(0,xy.limits[2]), expand = c(0,0.5)) +
  labs(x="GPP Data",y="Predicted GPP")+
  geom_abline(slope = 1, intercept = 0)+
  theme_classic()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12))






  
}







