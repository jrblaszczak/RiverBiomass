## Conceptual Figure

# load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","patchwork"), require, character.only=T)

## Source data
source("DataSource_9rivers.R")

dat <- df$nwis_08180700


## Time series
GPP_TS <- ggplot(dat, aes(date, GPP))+
  geom_point(size=2, color="chartreuse4")+
  geom_errorbar(aes(ymin=GPP.lower, ymax=GPP.upper), color="chartreuse4")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1))+
  coord_cartesian(ylim=c(0,8))

ratio_QL <- max(dat$light)/max(dat$Q)

external_TS <- ggplot(dat, aes(date, Q*ratio_QL))+
  geom_point(data=dat, aes(date, light), size=1.5, color="darkgoldenrod3")+
  geom_line(size=1, color="deepskyblue4")+
  scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
  labs(y=expression('Daily PPFD'), title="External Drivers")+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_blank(), axis.text = element_text(size=12),
        axis.title.y.left = element_text(size=14, color="darkgoldenrod3"),
        axis.title.y.right = element_text(size=14, color="deepskyblue4"),
        axis.text.x = element_text(angle=25, hjust = 1))

(GPP_TS + theme(axis.text.x = element_blank())) +
    external_TS + plot_layout(ncol=1)

## GPP versus light
ggplot(dat, aes(light, GPP))+
  geom_point(size=2, color="black")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily Light")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=15), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1))+
  geom_quantile(quantiles = 0.5, size=1.5)+
  coord_cartesian(ylim=c(0,8))
summary(lm(GPP ~ light, data=dat)) # Medina: slope = 0.003, intercept = 2.5


## predict daily GPP based on linear relationship with light
coeffs <- summary(lm(GPP ~ light, data=dat))$coeff
dat$GPP_linpred <- dat$GPP*coeffs[2,1] + coeffs[1,1]

ggplot(dat, aes(light, GPP_linpred))+
  geom_point(size=2, color="black")+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily Light")+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.title.x = element_text(size=15), axis.text = element_text(size=12),
        axis.title.y = element_text(size=15), axis.text.x = element_text(angle=25, hjust = 1))

sum(dat$GPP_linpred)
sum(dat$GPP)

## productivity hysteresis
dat[which(dat$Q > 30),]$date

storm_light_plot <- function(dat, start.date, end.date, GPP.upper.lim){
  
  storm1 <- subset(dat, date >= start.date & date <= end.date)
  ratio_GPPL <- max(storm1$GPP)/max(storm1$light)
  
  #storm1$roll3 <- zoo::rollmean(storm1$GPP, k=3)
  
  plot_grid(ggplot(storm1, aes(date, light*ratio_GPPL))+
              geom_polygon(data=storm1, aes(date, Q/25), fill="deepskyblue4", alpha=0.5)+
              geom_point(color="darkgoldenrod3", size=1.5) + geom_line(color="darkgoldenrod3", size=1)+
              geom_point(data=storm1, aes(date, GPP), size=1.5)+
              geom_line(data=storm1, aes(date, GPP), size=1)+
              scale_y_continuous(sec.axis = sec_axis(~./ratio_GPPL, name=expression("PPFD")))+
              coord_cartesian(ylim=c(0,GPP.upper.lim))+
              labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
              theme_bw(),
            
            ggplot(storm1, aes(light, GPP, color=date))+
              geom_point() + geom_path() +
              labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="PPFD")+
              theme_bw())

}

storm_light_plot(dat, "2010-05-10", "2010-06-10",5)
storm_light_plot(dat, "2010-09-05", "2010-10-05",3.5)

