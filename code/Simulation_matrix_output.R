# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse",
         "rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("Oregon_ClackamasData_Source.R")
df <- dat$Clackamas_OR
df <- df[which(df$year == "2010"),] ## Subset to 2010

# source simulation models
# input variables: GPP, GPP_sd, light, tQ
source("Simulated_ProductivityModel1_Autoregressive.R") # estimated parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # estimated parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # estimated parameters: r, beta_0, s, c, sig_p

# for parameter extraction
source("StanParameterExtraction_Source.R")

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"



##################
## Model 1 Output
##################
PM1_DataOutput <- readRDS("PM1_DataOutput.rds")

pars1<-extract(PM1_DataOutput, c("phi","alpha","beta","sig_p"))

simmat1<-matrix(NA,length(df[,1]),length(unlist(pars1$phi)))
rmsemat1<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars1$phi)){
  
  simmat1[,i]<-PM1(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p,df)
  rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
  
}

## View RMSE distribution
hist(rmsemat1, xlim = c(0,20))

## For every day extract median and CI
mean_simmat1 <- apply(simmat1, 1, function(x) mean(x))
lower_simmat1 <- apply(simmat1, 1, function(x) quantile(x, probs = 0.025))
upper_simmat1 <- apply(simmat1, 1, function(x) quantile(x, probs = 0.975))

## Plot simulated GPP
df_sim1 <- as.data.frame(cbind(as.character(df$date), df$GPP, mean_simmat1, lower_simmat1, upper_simmat1))
colnames(df_sim1) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
df_sim1[,2:5] <- apply(df_sim1[,2:5],2,function(x) as.numeric(as.character(x)))

df_sim1_plot <- ggplot(df_sim1, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM1.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: GPP")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM1.col, alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,30))
df_sim1_plot


##################
## Model 2 Output
##################
PM2_DataOutput <- readRDS("PM2_DataOutput.rds")

pars2<-extract(PM2_DataOutput, c("r","K","s","c","B","P","pred_GPP","sig_p"))
   
simmat2<-matrix(NA,length(df[,1]),length(unlist(pars2$r)))
rmsemat2<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars2$r)){
  
  simmat2[,i]<-PM2(pars2$r[i],pars2$K[i],pars2$s[i],pars2$c[i],pars2$sig_p,df)
  rmsemat2[i]<-sqrt(sum((simmat2[,i]-df$GPP)^2)/length(df$GPP))
  
}

## View RMSE distribution
hist(rmsemat2)

## For every day extract median and CI
mean_simmat2 <- apply(simmat2, 1, function(x) mean(x))
lower_simmat2 <- apply(simmat2, 1, function(x) quantile(x, probs = 0.025))
upper_simmat2 <- apply(simmat2, 1, function(x) quantile(x, probs = 0.975))

## Plot simulated GPP
df_sim2 <- as.data.frame(cbind(as.character(df$date), df$GPP, mean_simmat2, lower_simmat2, upper_simmat2))
colnames(df_sim2) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim2$Date <- as.POSIXct(as.character(df_sim2$Date), format="%Y-%m-%d")
df_sim2[,2:5] <- apply(df_sim2[,2:5],2,function(x) as.numeric(as.character(x)))

df_sim2_plot <- ggplot(df_sim2, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM2.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM2: Logistic")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM2.col, alpha=0.2, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,30))
df_sim2_plot


## Plot latent B
PM2_medpar <- mechB_extract_medians(rstan::extract(PM2_DataOutput,c("r","K","s","c","B","P","pred_GPP","sig_p")))
df_modB2 <- as.data.frame(cbind(as.character(df$date), PM2_medpar$B, PM2_medpar$B_Q.025, PM2_medpar$B_Q.975))
colnames(df_modB2) <- c("Date","B","B_lower","B_upper")
df_modB2$Date <- as.POSIXct(as.character(df_modB2$Date), format="%Y-%m-%d")
df_modB2[,2:4] <- apply(df_modB2[,2:4],2,function(x) as.numeric(as.character(x)))

df_modB2_plot <- ggplot(df_modB2, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM2: Logistic")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  scale_y_continuous(limits=c(0,35))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))
df_modB2_plot

##################
## Model 3 Output
##################
PM3_DataOutput <- readRDS("PM3_DataOutput.rds")

pars3<-extract(PM3_DataOutput, c("r","beta_0","s","c","B","P","pred_GPP","sig_p"))

simmat3<-matrix(NA,length(df[,1]),length(unlist(pars3$r)))
rmsemat3<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars3$r)){
  
  simmat3[,i]<-PM3(pars3$r[i],pars3$beta_0[i],pars3$s[i],pars3$c[i],pars3$sig_p,df)
  rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^3)/length(df$GPP))
  
}

## View RMSE distribution
hist(rmsemat3)

## For every day extract median and CI
mean_simmat3 <- apply(simmat3, 1, function(x) mean(x))
lower_simmat3 <- apply(simmat3, 1, function(x) quantile(x, probs = 0.025))
upper_simmat3 <- apply(simmat3, 1, function(x) quantile(x, probs = 0.975))

## Plot simulated GPP
df_sim3 <- as.data.frame(cbind(as.character(df$date), df$GPP, mean_simmat3, lower_simmat3, upper_simmat3))
colnames(df_sim3) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper")
df_sim3$Date <- as.POSIXct(as.character(df_sim3$Date), format="%Y-%m-%d")
df_sim3[,2:5] <- apply(df_sim3[,2:5],2,function(x) as.numeric(as.character(x)))

df_sim3_plot <- ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM3.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM3: Ricker")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM3.col, alpha=0.4, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,30))
df_sim3_plot


## Plot latent B
PM3_medpar <- mechB_extract_medians(rstan::extract(PM3_DataOutput, c("r","beta_0","s","c","B","P","pred_GPP","sig_p")))

df_modB3 <- as.data.frame(cbind(as.character(df$date), PM3_medpar$B, PM3_medpar$B_Q.025, PM3_medpar$B_Q.975))
colnames(df_modB3) <- c("Date","B","B_lower","B_upper")
df_modB3$Date <- as.POSIXct(as.character(df_modB3$Date), format="%Y-%m-%d")
df_modB3[,2:4] <- apply(df_modB3[,2:4],2,function(x) as.numeric(as.character(x)))

df_modB3_plot <- ggplot(df_modB3, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM3: Ricker")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  scale_y_continuous(limits=c(0,35))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))
df_modB3_plot





##################################################
## Compare simulations and RMSE of PM1 and PM2
###################################################
## simulation comparison
plot_grid(
  df_sim1_plot,
  df_sim2_plot,
  df_sim3_plot,
  ncol=1
)

plot_grid(
  df_modB2_plot,
  df_modB3_plot, ncol=1)

## RMSE comparison
rmse_comp <- as.data.frame(as.matrix(cbind(rmsemat1, rmsemat2, rmsemat3)))
colnames(rmse_comp) <- c("PM1 RMSE","PM2 RMSE", "PM3 RMSE")
rmse_comp_long <- gather(rmse_comp)

rmse_comp_mean <- rmse_comp_long %>%
  group_by(key) %>%
  summarise(rating.mean = mean(na.omit(value)))

ggplot(rmse_comp_long, aes(value, fill=key))+
  geom_density(alpha=0.3)+
  scale_fill_manual("",values=c("PM1 RMSE" = PM1.col,"PM2 RMSE" = PM2.col,"PM3 RMSE" = PM3.col))+
  scale_x_continuous(trans="log", limits=c(0.35,35),breaks = c(0.5,1,3,5,10,30), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15),
        legend.position = c(.25, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=14))+
  labs(y="Density","RMSE")+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[1],
             color=PM1.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[2],
             color=PM2.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[3],
             color=PM3.col, linetype = "dashed", size = 1)






#######################################
## figures modified for presentations
############################################

        