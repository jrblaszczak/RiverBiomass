## Figures for model fit and RMSE comparison

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("Oregon_ClackamasData_Source.R")
df <- dat$Clackamas_OR
df <- df[which(df$year == "2010"),] ## Subset to 2010
# temporary solution to turbidity data gaps
df$mean_daily_turb <- na.approx(df$mean_daily_turb, maxgap = 3)
colnames(df)[which(colnames(df) == "mean_daily_turb")] <- "turb"

# source simulation models
source("Simulated_ProductivityModel1_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Simulated_ProductivityModel4_Ricker_lightadj.R") # parameters: alpha_1, lambda, s, c, sig_p
source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, beta_2, s, c, sig_p

# for parameter extraction
source("StanParameterExtraction_Source.R")

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"
PM4.col <- "#598C8D"
PM5.col <- "#743731"

## Import stan fits
stan_model_output_AR <- readRDS("stan_model_output_AR.rds")
stan_model_output_Logistic <- readRDS("stan_model_output_Logistic.rds")
stan_model_output_Ricker <- readRDS("stan_model_output_Ricker.rds")
stan_model_output_Ricker_Ladj <- readRDS("stan_model_output_Ricker_Ladj.rds")
stan_model_output_Gompertz <- readRDS("stan_model_output_Gompertz.rds")

##########################
## Model 1 Output - AR
#########################

## No light modification ##
pars1 <- extract(stan_model_output_AR[[1]], c("phi","alpha","beta","sig_p"))
simmat1<-matrix(NA,length(df[,1]),length(unlist(pars1$phi)))
rmsemat1<-matrix(NA,length(df[,1]),1)
# Simulate
for (i in 1:length(pars1$phi)){
  simmat1[,i]<-PM1(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p[i],df)
  rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
}

## Benthic light modification ##
pars1_BL <- extract(stan_model_output_AR[[2]], c("phi","alpha","beta","sig_p","a"))
simmat1_BL<-matrix(NA,length(df[,1]),length(unlist(pars1_BL$phi)))
rmsemat1_BL<-matrix(NA,length(df[,1]),1)
# Simulate
for (i in 1:length(pars1_BL$phi)){
  simmat1_BL[,i]<-PM1_BL(pars1_BL$phi[i],pars1_BL$alpha[i],pars1_BL$beta[i],pars1_BL$a[i],pars1_BL$sig_p,df)
  rmsemat1_BL[i]<-sqrt(sum((simmat1_BL[,i]-df$GPP)^2)/length(df$GPP))
}


## Remove and save simulation
rm(stan_model_output_AR)
simmat1_list <- list(simmat1, simmat1_BL)
saveRDS(simmat1_list, "Sim_matrix_AR.rds")
saveRDS(rmsemat1, "RMSE_matrix_AR.rds")

## If previously simulated
simmat1_list <- readRDS("Sim_matrix_AR.rds")
simmat1 <- simmat1_list[[1]]
simmat1_BL <- simmat1_list[[2]]

# For every day extract median and CI
median_simmat1 <- apply(simmat1, 1, function(x) median(x))
lower_simmat1 <- apply(simmat1, 1, function(x) quantile(x, probs = 0.025))
upper_simmat1 <- apply(simmat1, 1, function(x) quantile(x, probs = 0.975))

# For every day extract median and CI
median_simmat1_BL <- apply(simmat1_BL, 1, function(x) median(x))
lower_simmat1_BL <- apply(simmat1_BL, 1, function(x) quantile(x, probs = 0.025))
upper_simmat1_BL <- apply(simmat1_BL, 1, function(x) quantile(x, probs = 0.975))

## Plot simulated GPP
df_sim1 <- as.data.frame(cbind(as.character(df$date), df$GPP, median_simmat1,lower_simmat1, upper_simmat1,
                               median_simmat1_BL,lower_simmat1_BL, upper_simmat1_BL))
colnames(df_sim1) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper",
                       "sim_GPP_BL","sim_GPP_BL_lower","sim_GPP_BL_upper")
df_sim1$Date <- as.POSIXct(as.character(df_sim1$Date), format="%Y-%m-%d")
df_sim1[,2:8] <- apply(df_sim1[,2:8],2,function(x) as.numeric(as.character(x)))

df_sim1_plot <- ggplot(df_sim1, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM1.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: GPP")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM1.col, alpha=0.3, show.legend = FALSE)+
  #geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  #geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
  #            fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,30))
df_sim1_plot


df_sim1_BL_plot <- ggplot(df_sim1, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM1.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM1: GPP with benthic light in blue")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM1.col, alpha=0.3, show.legend = FALSE)+
  geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
              fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,30))
df_sim1_BL_plot

plot_grid(df_sim1_plot,
          df_sim1_BL_plot,
          ncol = 1)

## Save image



##############################
## Model 2 Output - Logistic
##############################
## No light modification ##
pars2<-extract(stan_model_output_Logistic[[1]], c("r","K","s","c","B","P","pred_GPP","sig_p"))
simmat2<-matrix(NA,length(df[,1]),length(unlist(pars2$sig_p)))
rmsemat2<-matrix(NA,length(df[,1]),1)
#Simulate
for (i in 1:length(pars2$sig_p)){
  simmat2[,i]<-PM2(r=pars2$r[i],K=pars2$K[i],s=pars2$s[i],c=pars2$c[i],sig_p=pars2$sig_p[i],df)
  rmsemat2[i]<-sqrt(sum((simmat2[,i]-df$GPP)^2)/length(df$GPP))
}

## Benthic light modification ##
pars2_BL<-extract(stan_model_output_Logistic[[2]], c("r","K","s","c","a","B","P","pred_GPP","sig_p"))
simmat2_BL<-matrix(NA,length(df[,1]),length(unlist(pars2_BL$sig_p)))
rmsemat2_BL<-matrix(NA,length(df[,1]),1)
# Simulation
for (i in 1:length(pars2_BL$sig_p)){
  simmat2_BL[,i]<-PM2_BL(r=pars2_BL$r[i],K=pars2_BL$K[i],s=pars2_BL$s[i],c=pars2_BL$c[i],a=pars2_BL$a[i],sig_p=pars2_BL$sig_p[i],df)
  rmsemat2_BL[i]<-sqrt(sum((simmat2_BL[,i]-df$GPP)^2)/length(df$GPP))
}

## Save simulation
simmat2_list <- list(simmat2, simmat2_BL)
saveRDS(simmat2_list, "Sim_matrix_Logistic.rds")
saveRDS(rmsemat2, "RMSE_matrix_Logistic.rds")
# If already previously simulated
simmat2_list <- readRDS("Sim_matrix_Logistic.rds")
simmat2 <- simmat2_list[[1]]
simmat2_BL <- simmat2_list[[2]]

## For every day extract median and CI
median_simmat2 <- apply(simmat2, 1, function(x) median(x, na.rm = TRUE))
lower_simmat2 <- apply(simmat2, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upper_simmat2 <- apply(simmat2, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

## For every day extract median and CI
median_simmat2_BL <- apply(simmat2_BL, 1, function(x) median(x, na.rm = TRUE))
lower_simmat2_BL <- apply(simmat2_BL, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upper_simmat2_BL <- apply(simmat2_BL, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))




## Plot simulated GPP
df_sim2 <- as.data.frame(cbind(as.character(df$date), df$GPP, median_simmat2, lower_simmat2, upper_simmat2,
                               median_simmat2_BL, lower_simmat2_BL, upper_simmat2_BL))
colnames(df_sim2) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper",
                       "sim_GPP_BL","sim_GPP_BL_lower","sim_GPP_BL_upper")
df_sim2$Date <- as.POSIXct(as.character(df_sim2$Date), format="%Y-%m-%d")
df_sim2[,2:8] <- apply(df_sim2[,2:8],2,function(x) as.numeric(as.character(x)))

df_sim2_plot <- ggplot(df_sim2, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM2.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM2: Logistic")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM2.col, alpha=0.2, show.legend = FALSE)+
  #geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  #geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
  #            fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim2_plot


df_sim2_BL_plot <- ggplot(df_sim2, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM2.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM2: Logistic with benthic light in blue")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM2.col, alpha=0.2, show.legend = FALSE)+
  geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
              fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim2_BL_plot

plot_grid(df_sim2_plot,
          df_sim2_BL_plot,
          ncol = 1)


## Plot latent B
## no light modification ##
PM2_medpar <- mechB_extract_medians(rstan::extract(stan_model_output_Logistic[[1]],c("r","K","s","c","B","P","pred_GPP","sig_p")))
df_modB2 <- as.data.frame(cbind(as.character(df$date), PM2_medpar$B, PM2_medpar$B_Q.025, PM2_medpar$B_Q.975))
colnames(df_modB2) <- c("Date","B","B_lower","B_upper")
df_modB2$Date <- as.POSIXct(as.character(df_modB2$Date), format="%Y-%m-%d")
df_modB2[,2:4] <- apply(df_modB2[,2:4],2,function(x) as.numeric(as.character(x)))

## benthic light ##
PM2_BL_medpar <- mechB_extract_medians(rstan::extract(stan_model_output_Logistic[[2]],c("r","K","s","c","a","B","P","pred_GPP","sig_p")))
df_modB2_BL <- as.data.frame(cbind(as.character(df$date), PM2_BL_medpar$B, PM2_BL_medpar$B_Q.025, PM2_BL_medpar$B_Q.975))
colnames(df_modB2_BL) <- c("Date","B","B_lower","B_upper")
df_modB2_BL$Date <- as.POSIXct(as.character(df_modB2_BL$Date), format="%Y-%m-%d")
df_modB2_BL[,2:4] <- apply(df_modB2_BL[,2:4],2,function(x) as.numeric(as.character(x)))

df_modB2_plot <- ggplot(df_modB2, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM2: Logistic")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  geom_line(data=df_modB2_BL, aes(Date, exp(B)), size=1.2, color="grey45")+
  geom_ribbon(data=df_modB2_BL, aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="grey45", alpha=0.3, show.legend = FALSE)+
  scale_y_continuous(limits=c(0,35))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))
df_modB2_plot

rm(stan_model_output_Logistic)

###############################
## Model 3 Output - Ricker
###############################
## no light modification ##
pars3<-extract(stan_model_output_Ricker[[1]], c("r","lambda","s","c","B","P","pred_GPP","sig_p"))
simmat3<-matrix(NA,length(df[,1]),length(unlist(pars3$r)))
rmsemat3<-matrix(NA,length(df[,1]),1)
#Simulated
for (i in 1:length(pars3$r)){
  simmat3[,i]<-PM3(pars3$r[i],pars3$lambda[i],pars3$s[i],pars3$c[i],pars3$sig_p[i],df)
  rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
}

## benthic light modification ##
pars3_BL<-extract(stan_model_output_Ricker[[2]], c("r","lambda","s","c","a","B","P","pred_GPP","sig_p"))
simmat3_BL<-matrix(NA,length(df[,1]),length(unlist(pars3_BL$sig_p)))
rmsemat3_BL<-matrix(NA,length(df[,1]),1)
# Simulation
for (i in 1:length(pars3_BL$sig_p)){
  simmat3_BL[,i]<-PM3_BL(r=pars3_BL$r[i],lambda=pars3_BL$lambda[i],s=pars3_BL$s[i],c=pars3_BL$c[i],a=pars3_BL$a[i],sig_p=pars3_BL$sig_p[i],df)
  rmsemat3_BL[i]<-sqrt(sum((simmat3_BL[,i]-df$GPP)^2)/length(df$GPP))
}


## Save simulation
simmat3_list <- list(simmat3, simmat3_BL)
saveRDS(simmat3_list, "Sim_matrix_Ricker.rds")
saveRDS(rmsemat3, "RMSE_matrix_Ricker.rds")
# If already previously simulated
simmat3_list <- readRDS("Sim_matrix_Ricker.rds")
simmat3 <- simmat3_list[[1]]
simmat3_BL <- simmat3_list[[2]]

## For every day extract median and CI
median_simmat3 <- apply(simmat3, 1, function(x) median(x))
lower_simmat3 <- apply(simmat3, 1, function(x) quantile(x, probs = 0.025))
upper_simmat3 <- apply(simmat3, 1, function(x) quantile(x, probs = 0.975))

# For every day extract median and CI
median_simmat3_BL <- apply(simmat3_BL, 1, function(x) median(x, na.rm = TRUE))
lower_simmat3_BL <- apply(simmat3_BL, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upper_simmat3_BL <- apply(simmat3_BL, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

## Plot simulated GPP
df_sim3 <- as.data.frame(cbind(as.character(df$date), df$GPP, median_simmat3, lower_simmat3, upper_simmat3,
                               median_simmat3_BL, lower_simmat3_BL, upper_simmat3_BL))
colnames(df_sim3) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper",
                       "sim_GPP_BL","sim_GPP_BL_lower","sim_GPP_BL_upper")
df_sim3$Date <- as.POSIXct(as.character(df_sim3$Date), format="%Y-%m-%d")
df_sim3[,2:8] <- apply(df_sim3[,2:8],2,function(x) as.numeric(as.character(x)))

df_sim3_plot <- ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM3.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM3: Ricker")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM3.col, alpha=0.4, show.legend = FALSE)+
  #geom_line(aes(Date, sim_GPP_BL), color="#480711", size=1.2)+
  #geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
  #            fill="#480711", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim3_plot

df_sim3_BL_plot <- ggplot(df_sim3, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM3.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM3: Ricker with benthic light in blue")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM3.col, alpha=0.4, show.legend = FALSE)+
  geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
              fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim3_BL_plot

plot_grid(df_sim3_plot,
          df_sim3_BL_plot,
          ncol=1)


## Plot latent B
PM3_medpar <- mechB_extract_medians(rstan::extract(stan_model_output_Ricker[[1]], c("r","lambda","s","c","B","P","pred_GPP","sig_p")))

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


############################################
## Model 4 Output - Ricker, adjusted light
#############################################
## no light modifications ##
pars4<-extract(stan_model_output_Ricker_Ladj[[1]], c("alpha_1","lambda","s","c","r","B","P","pred_GPP","sig_p"))
simmat4<-matrix(NA,length(df[,1]),length(unlist(pars4$alpha_1)))
rmsemat4<-matrix(NA,length(df[,1]),1)
#Simulate
for (i in 1:length(pars4$alpha_1)){
  simmat4[,i]<-PM4(pars4$alpha_1[i],pars4$lambda[i],pars4$s[i],pars4$c[i],pars4$sig_p[i],df)
  rmsemat4[i]<-sqrt(sum((simmat4[,i]-df$GPP)^2)/length(df$GPP))
}

## benthiclight modifications ##
pars4_BL<-extract(stan_model_output_Ricker_Ladj[[2]], c("alpha_1","lambda","s","c","r","a","B","P","pred_GPP","sig_p"))
simmat4_BL<-matrix(NA,length(df[,1]),length(unlist(pars4_BL$alpha_1)))
rmsemat4_BL<-matrix(NA,length(df[,1]),1)
#Simulate
for (i in 1:length(pars4_BL$alpha_1)){
  simmat4_BL[,i]<-PM4_BL(pars4_BL$alpha_1[i],pars4_BL$lambda[i],pars4_BL$s[i],pars4_BL$c[i],pars4_BL$a[i],pars4_BL$sig_p[i],df)
  rmsemat4_BL[i]<-sqrt(sum((simmat4_BL[,i]-df$GPP)^2)/length(df$GPP))
}


## For every day extract median and CI
median_simmat4 <- apply(simmat4, 1, function(x) median(x, na.rm = TRUE))
lower_simmat4 <- apply(simmat4, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upper_simmat4 <- apply(simmat4, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))

## For every day extract median and CI
median_simmat4_BL <- apply(simmat4_BL, 1, function(x) median(x, na.rm = TRUE))
lower_simmat4_BL <- apply(simmat4_BL, 1, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
upper_simmat4_BL <- apply(simmat4_BL, 1, function(x) quantile(x, probs = 0.975, na.rm = TRUE))


## Save simulation
simmat4_list <- list(simmat4, simmat4_BL)
saveRDS(simmat4_list, "Sim_matrix_Ricker_lightadj.rds")
saveRDS(rmsemat4, "RMSE_matrix_Ricker_lightadj.rds")
# If already previously simulated
simmat4_list <- readRDS("Sim_matrix_Ricker_lightadj.rds")
simmat4 <- simmat4_list[[1]]
simmat4_BL <- simmat4_list[[2]]



## Plot simulated GPP
df_sim4 <- as.data.frame(cbind(as.character(df$date), df$GPP, median_simmat4, lower_simmat4, upper_simmat4,
                               median_simmat4_BL, lower_simmat4_BL, upper_simmat4_BL))
colnames(df_sim4) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper",
                       "sim_GPP_BL","sim_GPP_BL_lower","sim_GPP_BL_upper")
df_sim4$Date <- as.POSIXct(as.character(df_sim4$Date), format="%Y-%m-%d")
df_sim4[,2:8] <- apply(df_sim4[,2:8],2,function(x) as.numeric(as.character(x)))

df_sim4_plot <- ggplot(df_sim4, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM4.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM4: Ricker with light adjustment")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM4.col, alpha=0.4, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim4_plot


df_sim4_BL_plot <- ggplot(df_sim4, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM4.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM4: Ricker with light adjustment and benthic light in blue")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM4.col, alpha=0.4, show.legend = FALSE)+
  geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
              fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,20))
df_sim4_BL_plot

plot_grid(df_sim4_plot,
          df_sim4_BL_plot,
          ncol=1)


## Plot latent B
PM4_medpar <- mechB_extract_medians(rstan::extract(stan_model_output_Ricker_Ladj[[1]], c("alpha_0","alpha_1","lambda","s","c","r","B","P","pred_GPP","sig_p")))

df_modB4 <- as.data.frame(cbind(as.character(df$date), PM4_medpar$B, PM4_medpar$B_Q.025, PM4_medpar$B_Q.975))
colnames(df_modB4) <- c("Date","B","B_lower","B_upper")
df_modB4$Date <- as.POSIXct(as.character(df_modB4$Date), format="%Y-%m-%d")
df_modB4[,2:4] <- apply(df_modB4[,2:4],2,function(x) as.numeric(as.character(x)))

df_modB4_plot <- ggplot(df_modB4, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM4: Ricker with light adjustment")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  scale_y_continuous(limits=c(0,35))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))
df_modB4_plot


###############################
## Model 5 Output - Gompertz
###############################
## no light modification ##
pars5<-extract(stan_model_output_Gompertz[[1]], c("beta_0","beta_1","beta_2","s","c","B","P","pred_GPP","sig_p"))
simmat5<-matrix(NA,length(df[,1]),length(unlist(pars5$beta_0)))
rmsemat5<-matrix(NA,length(df[,1]),1)
#Simulate
for (i in 1:length(pars5$beta_0)){
  simmat5[,i]<-PM5(pars5$beta_0[i],pars5$beta_1[i],pars5$beta_2[i],pars5$s[i],pars5$c[i],pars5$sig_p[i],df)
  rmsemat5[i]<-sqrt(sum((simmat5[,i]-df$GPP)^2)/length(df$GPP))
}

## benthic light modification ##
pars5_BL<-extract(stan_model_output_Gompertz[[2]], c("beta_0","beta_1","beta_2","s","c","a","B","P","pred_GPP","sig_p"))
simmat5_BL<-matrix(NA,length(df[,1]),length(unlist(pars5_BL$beta_0)))
rmsemat5_BL<-matrix(NA,length(df[,1]),1)
#Simulate
for (i in 1:length(pars5_BL$beta_0)){
  simmat5_BL[,i]<-PM5_BL(pars5_BL$beta_0[i],pars5_BL$beta_1[i],pars5_BL$beta_2[i],pars5_BL$s[i],pars5_BL$c[i],pars5_BL$a[i],pars5_BL$sig_p[i],df)
  rmsemat5_BL[i]<-sqrt(sum((simmat5_BL[,i]-df$GPP)^2)/length(df$GPP))
}

## Save simulation
simmat5_list <- list(simmat5, simmat5_BL)
saveRDS(simmat5_list, "Sim_matrix_Gompertz.rds")
saveRDS(rmsemat5, "RMSE_matrix_Gompertz.rds")
# If already previously simulated
simmat5_list <- readRDS("Sim_matrix_Gompertz.rds")
simmat5 <- simmat5_list[[1]]
simmat5_BL <- simmat5_list[[2]]

## For every day extract median and CI
median_simmat5 <- apply(simmat5, 1, function(x) median(x, na.rm=TRUE))
lower_simmat5 <- apply(simmat5, 1, function(x) quantile(x, probs = 0.025, na.rm=TRUE))
upper_simmat5 <- apply(simmat5, 1, function(x) quantile(x, probs = 0.975, na.rm=TRUE))
## For every day extract median and CI
median_simmat5_BL <- apply(simmat5_BL, 1, function(x) median(x, na.rm=TRUE))
lower_simmat5_BL <- apply(simmat5_BL, 1, function(x) quantile(x, probs = 0.025, na.rm=TRUE))
upper_simmat5_BL <- apply(simmat5_BL, 1, function(x) quantile(x, probs = 0.975, na.rm=TRUE))


## Plot simulated GPP
df_sim5 <- as.data.frame(cbind(as.character(df$date), df$GPP, median_simmat5, lower_simmat5, upper_simmat5,
                               median_simmat5_BL, lower_simmat5_BL, upper_simmat5_BL))
colnames(df_sim5) <- c("Date","GPP","sim_GPP","sim_GPP_lower","sim_GPP_upper",
                       "sim_GPP_BL","sim_GPP_BL_lower","sim_GPP_BL_upper")
df_sim5$Date <- as.POSIXct(as.character(df_sim5$Date), format="%Y-%m-%d")
df_sim5[,2:8] <- apply(df_sim5[,2:8],2,function(x) as.numeric(as.character(x)))


df_sim5_plot <- ggplot(df_sim5, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM5.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM5: Gompertz")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM5.col, alpha=0.4, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,35))
df_sim5_plot

df_sim5_BL_plot <- ggplot(df_sim5, aes(Date, GPP))+
  geom_point(size=2, color="black")+
  geom_line(aes(Date, sim_GPP), color=PM5.col, size=1.2)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'),title="PM5: Gompertz benthic light in blue")+
  geom_ribbon(aes(ymin=sim_GPP_lower,ymax=sim_GPP_upper),
              fill=PM5.col, alpha=0.4, show.legend = FALSE)+
  geom_line(aes(Date, sim_GPP_BL), color="blue", size=1.2)+
  geom_ribbon(aes(ymin=sim_GPP_BL_lower,ymax=sim_GPP_BL_upper),
              fill="blue", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))+
  scale_y_continuous(limits=c(0,35))
df_sim5_BL_plot

plot_grid(df_sim5_plot,
          df_sim5_BL_plot,
          ncol=1)




## Plot latent B
PM5_medpar <- mechB_extract_medians(rstan::extract(stan_model_output_Gompertz[[1]], c("beta_0","beta_1","beta_2","s","c","B","P","pred_GPP","sig_p")))

df_modB5 <- as.data.frame(cbind(as.character(df$date), PM5_medpar$B, PM5_medpar$B_Q.025, PM5_medpar$B_Q.975))
colnames(df_modB5) <- c("Date","B","B_lower","B_upper")
df_modB5$Date <- as.POSIXct(as.character(df_modB5$Date), format="%Y-%m-%d")
df_modB5[,2:4] <- apply(df_modB5[,2:4],2,function(x) as.numeric(as.character(x)))

df_modB5_plot <- ggplot(df_modB5, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="PM5: Gompertz")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  scale_y_continuous(limits=c(0,35))+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15))
df_modB5_plot















##################################################
## Compare simulations and RMSE of all models
###################################################
## simulation comparison
plot_grid(
  df_sim1_plot,
  df_sim2_plot,
  df_sim3_plot,
  df_sim4_plot,
  df_sim5_plot,
  ncol=1
)

plot_grid(
  df_sim3_plot+scale_y_continuous(limits=c(0,15)),
  df_modB3_plot, ncol=1)

## RMSE comparison
#rmsemat_list <- readRDS("rmsemat_list.rds")
rmse_comp <- as.data.frame(as.matrix(cbind(rmsemat1, rmsemat2, rmsemat3, rmsemat4, rmsemat5)))
colnames(rmse_comp) <- c("PM1: GPP RMSE","PM2: Logistic RMSE",
                         "PM3: Ricker RMSE","PM4: Ricker Light Adj. RMSE","PM5: Gompertz RMSE")
rmse_comp_long <- gather(rmse_comp)

rmse_comp_mean <- rmse_comp_long %>%
  group_by(key) %>%
  summarise(rating.mean = mean(na.omit(value)))

ggplot(rmse_comp_long, aes(value, fill=key))+
  geom_density(alpha=0.5)+
  scale_fill_manual("",values=c("PM1: GPP RMSE" = PM1.col,"PM2: Logistic RMSE" = PM2.col,
                                "PM3: Ricker RMSE" = PM3.col,"PM4: Ricker Light Adj. RMSE" = PM4.col,
                                "PM5: Gompertz RMSE" = PM5.col))+
  scale_x_continuous(trans="log", limits=c(1,40), breaks = c(1,3,5,10,30), expand = c(0.01,0.01))+
  scale_y_continuous(expand = c(0,0.01))+
  theme(panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=15), axis.text = element_text(size=13),
        axis.title.y = element_text(size=15),
        legend.position = c(.5, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.text = element_text(size=14))+
  labs(y="Density",x="Daily RMSE")+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[1],
             color=PM1.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[2],
             color=PM2.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[3],
             color=PM3.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[4],
             color=PM4.col, linetype = "dashed", size = 1)+
  geom_vline(xintercept = rmse_comp_mean$rating.mean[5],
             color=PM5.col, linetype = "dashed", size = 1)





        