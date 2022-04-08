##==============================================================================
## Script for extracting information for results
## Code author: J.R. Blaszczak
##==============================================================================

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate","parallel",
         "tidyverse","rstan","bayesplot","shinystan","Metrics","MCMCglmm",
         "wesanderson"), require, character.only=T)

## Source data
source("DataSource_6rivers_StreamLight.R")

########################################
## Parameter summaries
########################################
## Import stan model fit
stan_model_output_STS <- readRDS("./rds files/stan_6riv_output_AR_2022_02_22.rds")
stan_model_output_LBTS <- readRDS("./rds files/stan_6riv_output_Ricker_2022_02_27.rds")

stan_psum <- function(x){
  
  s_init <- as.data.frame(summary(x)$summary)
  s_init$pars <- row.names(s_init)
  s_init$n_eff_pct <- s_init$n_eff/10000  ## effective samples are the number of independent samples with the same estimation power as the N autocorrelated samples
  s_init$n_eff_less10pct <- ifelse(s_init$n_eff_pct < 0.10, yes = "true", no = "false") # 10% is often used as a threshold, below which the chains for a parameter did not properly converge
  s <- s_init[,c("pars","50%", "2.5%", "97.5%","Rhat","n_eff","n_eff_less10pct")]

  return(s)
  
}

##
## YEAR 1 ##
##
pSTS <- ldply(lapply(stan_model_output_STS, function(z) stan_psum(z)), data.frame)
#STS params of interest: c("phi","alpha","beta","sig_p","sig_o")
write.csv(pSTS, "./tables/STS_ws_posterior_sum.csv")
pSTS_sub <- pSTS[which(pSTS$pars %in% c("phi","alpha","beta","sig_p","sig_o")),]
write.csv(pSTS_sub, "./tables/STS_ws_posteriorsubset_sum.csv")

pLBTS <- ldply(lapply(stan_model_output_LBTS, function(z) stan_psum(z)), data.frame)
#LB-TS params of interest: c("r","lambda","s","c","sig_p","sig_o")
write.csv(pLBTS, "./tables/LBTS_ws_posterior_sum.csv")
pLBTS_sub <- pLBTS[which(pLBTS$pars %in% c("r","lambda","s","c","sig_p","sig_o")),]
write.csv(pLBTS_sub, "./tables/LBTS_ws_posteriorsubset_sum.csv")

##
## YEAR 2 ##
##
Yr2_output_STS <- readRDS("./rds files/stan_6riv_2ndYr_output_AR_2022_03_06.rds")
Yr2_output_LBTS <- readRDS("./rds files/stan_6riv_2ndYr_output_Ricker_2022_03_06.rds")
## need stan_psum function from above

pSTS2 <- ldply(lapply(Yr2_output_STS, function(z) stan_psum(z)), data.frame)
#STS params of interest: c("phi","alpha","beta","sig_p","sig_o")
write.csv(pSTS2, "./tables/STS_ws_posterior_sum_Yr2.csv")
pSTS2_sub <- pSTS2[which(pSTS2$pars %in% c("phi","alpha","beta","sig_p","sig_o")),]
write.csv(pSTS2_sub, "./tables/STS_ws_posteriorsubset_sum_Yr2.csv")

pLBTS2 <- ldply(lapply(Yr2_output_LBTS, function(z) stan_psum(z)), data.frame)
#LB-TS params of interest: c("r","lambda","s","c","sig_p","sig_o")
write.csv(pLBTS2, "./tables/LBTS_ws_posterior_sum_Yr2.csv")
pLBTS2_sub <- pLBTS2[which(pLBTS2$pars %in% c("r","lambda","s","c","sig_p","sig_o")),]
write.csv(pLBTS2_sub, "./tables/LBTS_ws_posteriorsubset_sum_Yr2.csv")

######################
## S-TS description
######################


######################
## LB-TS description
######################


#############################################
## Critical flow thresholds and sensitivity
#############################################
##
## 1 - Magnitude of c relative to other flow metrics
##

## Import and merge bankfull discharge with site info
RI_2 <- read.csv("../data/RI_2yr_flood_6riv.csv", header=T)
sapply(RI_2, class)
site_info <- merge(site_info, RI_2, by="site_name")

## Reimport c estimates (within-sample) if not already loaded
pLBTS_sub <- read.csv("./tables/LBTS_ws_posteriorsubset_sum.csv", header=T)
c_sites <- pLBTS_sub[which(pLBTS_sub$pars == "c"),]
colnames(c_sites)[2] <- "site_name"

## Convert c estimate to discharge
c_Qc <- function(site, df, site_info, c_sites){
  
  Q_sub <- df[[site]]
  
  ## critical Q based on velocity
  crit_Q <- site_info[which(site_info$site_name == site),]$RI_2yr_Q
  
  ## convert c to critical Q based on GPP (Qc)
  Qc <- c_sites[which(c_sites$site_name == site),]$X50.*max(Q_sub$Q, na.rm = T)
  Qc_lower <- c_sites[which(c_sites$site_name == site),]$X2.5.*max(Q_sub$Q, na.rm = T)
  Qc_upper <- c_sites[which(c_sites$site_name == site),]$X97.5.*max(Q_sub$Q, na.rm = T)
  Qc_lmu <- as.data.frame(cbind(Qc,Qc_lower, Qc_upper, site))
  
  return(Qc_lmu)

}

## convert across all sites
site_list <- levels(as.factor(site_info$site_name))
Qc_all <- ldply(lapply(site_list, function(x) c_Qc(x,df,site_info,c_sites)), data.frame)
colnames(Qc_all) <- c("Qc_cms","Qc_lower_cms","Qc_upper_cms", "site_name")

## merge with site_info, convert 2 year flood to cms (divide by 35.314666212661), and calculate differences
crits <- merge(site_info[,c("site_name","RI_2yr_Q","short_name")],Qc_all, by="site_name")
crits$Q_2yrRI_cms <- crits$RI_2yr_Q/35.314666212661
## what is the max Q observed in the site
maxQ_obs <- ldply(lapply(df, function(x) max(x$Q)), data.frame)
colnames(maxQ_obs) <- c("site_name","Q_maxobs_cms")
crits <- merge(crits,maxQ_obs,by="site_name")
sapply(crits,class)
crits[,c("Qc_cms","Qc_lower_cms","Qc_upper_cms")] <- apply(crits[,c("Qc_cms","Qc_lower_cms","Qc_upper_cms")], 2, function(x) as.numeric(x))

## visualize
Qc_plot_df <- gather(crits[,c("short_name","Qc_cms","Q_2yrRI_cms","Q_maxobs_cms")], Q_type, Q_cms, Qc_cms:Q_maxobs_cms)
## Add in uncertainty estimates
cUI <- crits[,c("short_name","Qc_lower_cms","Qc_upper_cms")]
cUI$Q_type <- "Qc_cms"
Qc_plot_df <- merge(Qc_plot_df, cUI, by=c("short_name","Q_type"), all = TRUE)

##order
Qc_plot_df$short_name <- factor(Qc_plot_df$short_name, levels= site_order_list)
Qc_plot_df$Q_type <- factor(Qc_plot_df$Q_type, levels= c("Qc_cms","Q_maxobs_cms","Q_2yrRI_cms"))
## colors
Qcol <- c(wes_palette("Moonrise2")[2],wes_palette("Moonrise2")[1],wes_palette("Moonrise2")[4])

## plot
ggplot(Qc_plot_df, aes(color=Q_type, y=Q_cms+1, x=short_name)) + 
  geom_point(size=4)+
  geom_errorbar(aes(ymin = (Qc_lower_cms+1), ymax = (Qc_upper_cms+1)),
                width=0.2, size=0.8)+
  scale_y_continuous(trans="log", breaks=c(0.1+1, 10+1, 100+1, 1000+1),
                     labels=c("0.1", "10", "100","1,000"),
                     limits = c(0.1+1, 1000+1))+
  xlab("River") + ylab("Discharge (cms)")+
  scale_color_manual("", values = c("Qc_cms" = Qcol[1],
                                        "Q_maxobs_cms" = Qcol[2],
                                        "Q_2yrRI_cms" = Qcol[3]),
                    labels = c("Qc_cms" = expression(paste(Q[c])),
                               "Q_maxobs_cms" = expression(paste("Observed ",Q[max])),
                               "Q_2yrRI_cms" = expression(paste("Estimated ",Q["2yr"]))))+
    theme(panel.background = element_rect(fill = "white", color="black"),
        panel.grid = element_line(color = "gray85", linetype = "dashed", size = 0.2),
        axis.text.x = element_text(size=15, angle=45, hjust=1), 
        axis.text.y = element_text(size=15),
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        legend.position = c(0.85,0.85),
        legend.text.align = 0,
        legend.box.background = element_rect(color = "black"),
        axis.title = element_text(size=20))



## Quantify differences
crits$Qc_Qmax_pct <- crits$Qc_cms/crits$Q_maxobs_cms
crits$Qc_Q2yr_pct <- crits$Qc_cms/crits$Q_2yrRI_cms
## ranges
range(crits$Qc_Qmax_pct)*100; median(crits$Qc_Qmax_pct)*100 # range: 14 - 101%; median: 36%
range(crits$Qc_Q2yr_pct)*100; median(crits$Qc_Q2yr_pct)*100 # range: 3 - 72%; median: 16%
## are smaller or larger streams more sensitive? which sites have trouble converging c
c_eval <- merge(crits, c_sites, by="site_name")


## visualize
Qc_pct_df <- gather(crits[,c("short_name","Qc_Qmax_pct","Qc_Q2yr_pct")], pct_type, Pct, Qc_Qmax_pct:Qc_Q2yr_pct)
Qc_pct_df$short_name <- factor(Qc_pct_df$short_name, levels= site_order_list)
Qc_pct_df$pct_type <- factor(Qc_pct_df$pct_type, levels= c("Qc_Qmax_pct","Qc_Q2yr_pct"))

Qpct_col <- c(wes_palette("Moonrise2")[1],wes_palette("Moonrise2")[4])

ggplot(Qc_pct_df, aes(x=short_name,y=Pct, color = pct_type)) +
  geom_point(size=4)+#, position = position_dodge(width=0.2))+
  scale_y_continuous(labels = function(x) paste0(x * 100, '%'))+
  xlab("River") + ylab("Percent")+
  scale_color_manual("", values = c("Qc_Qmax_pct" = Qpct_col[1],
                                   "Qc_Q2yr_pct" = Qpct_col[2]),
                    labels = c("Qc_Qmax_pct" = expression(paste("% ",Q[c],"/",Q[max])),
                               "Qc_Q2yr_pct" = expression(paste("%",Q[c],"/",Q["2yr"]))))+
  theme(panel.background = element_rect(fill = "white", color="black"),
        panel.grid = element_line(color = "gray85", linetype = "dashed", size = 0.2),
        axis.text.x = element_text(size=15, angle=45, hjust=1), 
        axis.text.y = element_text(size=15),
        legend.text = element_text(size=18), legend.position = c(0.15,0.85),
        axis.title = element_text(size=20))

##
## 2 - Patterns in s
##

## evaluate patterns in s
s_sites <- pLBTS_sub[which(pLBTS_sub$pars == "s"),]
colnames(s_sites)[2] <- "site_name"
s_sites <- merge(s_sites, site_info[,c("site_name","short_name")], by="site_name")
range(s_sites$X50.); median(s_sites$X50.) # range: 1.26 - 1.74; median: 1.5

## visualize
s_sites$short_name <- factor(s_sites$short_name, levels= site_order_list)
ggplot(s_sites, aes(x = short_name, y = X50.))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin=X2.5., ymax=X97.5.), width=0.2)+
  theme_bw()

##
## 3 - Change in c between first and second year
##

## Import stan_model_output_LBTS (Year 1) and Yr2_output_LBTS (Year 2) model fits
## Calculate median values and convert to discharge

## Year 1
Qc_all_Yr1 <- Qc_all
Qc_all_Yr1$Year <- "Year 1"
## Year 2
## Source data - will overwrite year 1, re-import year 1
source("DataSource_6rivers_2ndYr_StreamLight.R")
df_yr2 <- df
source("DataSource_6rivers_StreamLight.R")
df_yr1 <- df

## Use c_Qc function to convert Year 2 data
Qc_all_Yr2 <- ldply(lapply(site_list, function(x) c_Qc(x,df_yr2,site_info,c_sites)), data.frame)
colnames(Qc_all_Yr2) <- c("Qc_cms","Qc_lower_cms","Qc_upper_cms", "site_name")
Qc_all_Yr2$Year <- "Year 2"

# Combine and restructure
Qc_all_Yrs <- rbind(Qc_all_Yr1, Qc_all_Yr2)
Qc_all_Yrs <- Qc_all_Yrs[,c("site_name", "Year","Qc_cms")]
Qc_wide <- spread(Qc_all_Yrs, Year, Qc_cms)
Qc_wide[,c("Year 1","Year 2")] <- apply(Qc_wide[,c("Year 1","Year 2")], 2, function(x) as.numeric(x))

Qc_wide$diff <- Qc_wide$`Year 2` - Qc_wide$`Year 1`
Qc_wide$pct.diff <- (Qc_wide$diff/Qc_wide$`Year 1`)*100
Qc_wide <- merge(Qc_wide, site_info, by="site_name")

ggplot(Qc_wide, aes(short_name, pct.diff))+geom_point()

## Plot differences in posterior c distributions

Yr1_c_ppd <- lapply(stan_model_output_LBTS, function(x) extract(x, c("c")))
Yr2_c_ppd <- lapply(Yr2_output_LBTS, function(x) extract(x, c("c")))



######################
## Hysteresis
######################
## See hysteresis code





















