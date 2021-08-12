## Within-sample daily latent biomass fit

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")
## Source parameter extraction code
source("StanParameterExtraction_Source.R")

## Import stan Ricker model fit
stan_model_output_Ricker <- readRDS("./rds files/stan_6riv_output_Ricker_2021_06_01.rds")

##################
## Plot biomass
##################
dat <- ldply(df, data.frame)

PM3_medpar <- ldply(lapply(stan_model_output_Ricker,
                           function(x) mechB_extract_medians(rstan::extract(x,c("B","P","pred_GPP")))),
                    data.frame)

df_modB3 <- as.data.frame(cbind(dat$site_name, as.character(dat$date), PM3_medpar$B, PM3_medpar$B_Q.025, PM3_medpar$B_Q.975))
colnames(df_modB3) <- c("site_name","Date","B","B_lower","B_upper")
df_modB3$Date <- as.POSIXct(as.character(df_modB3$Date), format="%Y-%m-%d")
df_modB3[,3:5] <- apply(df_modB3[,3:5],2,function(x) as.numeric(as.character(x)))

## Arrange rivers by river order
df_modB3 <- left_join(df_modB3, site_info[,c("site_name","short_name")])
df_modB3$short_name <- factor(df_modB3$short_name, levels=site_order_list)


## Plot latent biomass predictions
ggplot(df_modB3, aes(Date, exp(B)))+
  geom_line(size=1.2, color="chartreuse4")+
  labs(y="Latent Biomass",title="LB-TS Year 1")+
  geom_ribbon(aes(ymin=exp(B_lower),ymax=exp(B_upper)),
              fill="chartreuse4", alpha=0.3, show.legend = FALSE)+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=15))+
  facet_wrap(~short_name, scales = "free_x", ncol = 2)

