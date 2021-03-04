## RMSE comparison figure

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "reshape2"), require, character.only=T)

# source simulation models
source("Simulated_ProductivityModel1_Autoregressive.R") # parameters: phi, alpha, beta, sig_p
source("Simulated_ProductivityModel2_Logistic.R") # parameters: r, K, s, c, sig_p
source("Simulated_ProductivityModel3_Ricker.R") # parameters: r, lambda, s, c, sig_p
source("Simulated_ProductivityModel5_Gompertz.R") # parameters: beta_0, beta_1, s, c, sig_p


##############################
## Data Import & Processing
## (modified 'DataSource_6rivers.R')
##############################
data <- readRDS("../data/NWIS_6site_subset.rds")
data$date <- as.POSIXct(as.character(data$date), format="%Y-%m-%d")

site_info <- readRDS("../data/NWIS_6siteinfo_subset.rds")
## Change river names to short names
site_info[,c("site_name","long_name","NHD_STREAMORDE")]
site_info$short_name <- revalue(as.character(site_info$site_name), replace = c("nwis_01649500"="Anacostia River, MD",
                                                                               "nwis_02234000"="St. John's River, FL",
                                                                               "nwis_03058000"="West Fork River, WV",
                                                                               "nwis_08180700"="Medina River, TX",
                                                                               "nwis_10129900"="Silver Creek, UT",
                                                                               "nwis_14211010"="Clackamas River, OR"))


## How many days of data per site per year
data$year <- year(data$date)
data_siteyears <- data %>%
  group_by(site_name, year) %>%
  tally()
## Select the *second* of each site year pair
data <- rbind(data[which(data$site_name == "nwis_08180700" & data$year %in% c(2011)),],
              data[which(data$site_name == "nwis_10129900" & data$year %in% c(2016)),],
              data[which(data$site_name == "nwis_03058000" & data$year %in% c(2015)),],
              data[which(data$site_name == "nwis_01649500" & data$year %in% c(2013)),],
              data[which(data$site_name == "nwis_14211010" & data$year %in% c(2013)),],
              data[which(data$site_name == "nwis_02234000" & data$year %in% c(2014)),])

## Set any GPP < 0 to a small value close to 0
data[which(data$GPP < 0),]$GPP <- sample(exp(-6):exp(-4), 1)

## Create a GPP SD; SD = (CI - mean)/1.96
data$GPP_sd <- (((data$GPP.upper - data$GPP)/1.96) + ((data$GPP.lower - data$GPP)/-1.96))/2

## visualize
ggplot(data, aes(date, GPP))+
  geom_point()+geom_line()+
  facet_wrap(~site_name,scales = "free_x")

## split list by ID
l <- split(data, data$site_name)

rel_LQT <- function(x){
  x$light_rel <- x$light/max(x$light)
  x$temp_rel <- x$temp/max(x$temp)
  x$tQ <- x$Q/max(x$Q)
  
  #x$std_light <- (x$light-mean(x$light))/sd(x$light)
  #x$std_temp <- (x$temp-mean(x$temp))/sd(x$temp)
  #x$tQ <- (x$Q-mean(x$Q))/sd(x$Q)
  x<-x[order(x$date),]
  return(x)
}

dat_oos <- lapply(l, function(x) rel_LQT(x))

rm(data,l, data_siteyears)

# colors
PM1.col <- "#d95f02"
PM2.col <- "#7570b3"
PM3.col <- "#1C474D"
PM4.col <- "#743731"


##################################################
## Compare RMSE of all models
###################################################
## extract RMSE
simmat1_list <- readRDS("Sim_6riv_AR_oos.rds")
simmat2_list <- readRDS("Sim_6riv_Logistic_oos.rds")
simmat3_list <- readRDS("Sim_6riv_Ricker_oos.rds")
simmat4_list <- readRDS("Sim_6riv_Gompertz_oos.rds")

## extract only daily rmse
rmsemat1 <- ldply(lapply(simmat1_list, function(x) return(x[[2]])), data.frame)
rmsemat2 <- ldply(lapply(simmat2_list, function(x) return(x[[2]])), data.frame)
rmsemat3 <- ldply(lapply(simmat3_list, function(x) return(x[[2]])), data.frame)
rmsemat4 <- ldply(lapply(simmat4_list, function(x) return(x[[2]])), data.frame)


## Combine
rmse_comp <- as.data.frame(as.matrix(cbind(rmsemat1, rmsemat2$X..i.., rmsemat3$X..i.., rmsemat4$X..i..)))
colnames(rmse_comp) <- c("site_name","PM1: GPP","PM2: Logistic",
                         "PM3: Ricker","PM4: Gompertz")
rmse_comp$site_name <- as.factor(rmse_comp$site_name)
rmse_comp_long <- melt(rmse_comp, id.vars = "site_name")
rmse_comp_long$value <- as.numeric(as.character(rmse_comp_long$value))

rmse_comp_mean <- rmse_comp_long %>%
  group_by(site_name, variable) %>%
  summarise(rating.mean = mean(na.omit(value)))

## Arrange rivers by river order
rmse_comp_long <- left_join(rmse_comp_long, site_info[,c("site_name","short_name")])
rmse_comp_long$short_name <- factor(rmse_comp_long$short_name, levels=c("Silver Creek, UT",
                                                            "Medina River, TX",
                                                            "Anacostia River, MD",
                                                            "West Fork River, WV",
                                                            "St. John's River, FL",
                                                            "Clackamas River, OR"))


## Visualize
ggplot(rmse_comp_long, aes(x=variable, y=value, fill=variable, group=variable))+
  geom_jitter(size=0.2, color="gray75", width=0.3)+
  geom_boxplot(alpha=0.5, outlier.shape = NA)+ 
  facet_wrap(~short_name, ncol = 2)+
  scale_fill_manual("",values=c("PM1: GPP" = PM1.col,"PM2: Logistic" = PM2.col,
                                "PM3: Ricker" = PM3.col,"PM4: Gompertz" = PM4.col))+
  labs(y="Out-of-Sample Daily RMSE")+
  theme(legend.position = "right",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1),
        strip.background = element_rect(fill="white", color="black"),
        strip.text = element_text(size=12))+
  scale_y_continuous(trans="log", limits=c(1,40), breaks = c(1,3,10,30), expand = c(0.01,0.01))


## Anova
rmse_list <- split(rmse_comp_long, rmse_comp_long$site_name)
rmse_anova <- lapply(rmse_list, function(x) aov(value ~ variable, data=x))
lapply(rmse_anova, function(x) summary(x))
rmse_tukey <- lapply(rmse_anova, function(x) TukeyHSD(x))
lapply(rmse_tukey, function(x) print(x))

plot(rmse_tukey$nwis_01649500, las=1)


