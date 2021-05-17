## Subset data from hypoxia database that is already linked to NHD
## JRB

## Load packages
lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "tidyverse","data.table","patchwork"), require, character.only=T)

############################
## To create linked file
############################

## Import site data from Appling
site <- fread("../data/site_data.csv")
names(site)

## Import hypoxia dataset and subset to Appling
hyp <- fread("../data/GRDO_GEE_HA_NHD_2021_02_07.csv")
PC <- hyp[which(hyp$DB_Source == "PC"),]
names(PC)

## Merge both based on siteID
head(site$site_name); head(PC$SiteID)
colnames(PC)[which(colnames(PC) == "SiteID")] <- "site_name"

df <- left_join(site, PC, by="site_name")

## Export
write.csv(df, "../data/PC_site_attribs.csv")

##########################################
## Re-import for site selection
##########################################
df <- fread("../data/PC_site_attribs.csv")
colnames(df)

## subset
sub <- df[,c("site_name","long_name","NHD_STREAMORDE","US_state",
           "site_type","struct.canal_flag","struct.dam_flag","struct.npdes_flag")]

##################################################
## Select sites based on data quality
##################################################

# Import and subset model diagnostics
diagnostics <- read.table("../data/diagnostics.tsv",sep = "\t", header=T)
diagnostics <- diagnostics[which(diagnostics$site %in% sub$site_name),]
high_sites <- unique(diagnostics[which(diagnostics$model_confidence == "H"),]$site) ## 254

# Subset s based on high sites and site type and flags
s <- sub[which(sub$site_name %in% high_sites),] ## 254
s <- s[which(s$site_type == "ST"),] ## 249
s <- s[which(s$struct.dam_flag %in% c(NA,"95")),]

# Import time series
NWIS <- read.table("../data/daily_predictions.tsv", sep='\t', header = TRUE)
NWIS$date <- as.POSIXct(as.character(NWIS$date), format="%Y-%m-%d")
head(NWIS)

## Subset columns and sites
NWIS_sub <- NWIS[,c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                    "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                    "temp.water","discharge","shortwave","velocity")]
colnames(NWIS_sub) <- c("site_name","date","GPP","GPP.lower","GPP.upper", "GPP.Rhat",
                        "ER","ER.lower","ER.upper","K600","K600.lower","K600.upper",
                        "temp","Q","light","velocity")

## Subset to sites in high_sites (sites with high confidence rating)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% high_sites),]
# Confirm
length(levels(as.factor(NWIS_sub$site_name))) #254 sites

## Identify which sites have the most continuous data
NWIS_sub$doy <- yday(NWIS_sub$date)
NWIS_sub$year <- year(NWIS_sub$date)

## count days per year
dat_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  count()
## identify the max day gap per year
gap_per_year <- NWIS_sub %>%
  group_by(site_name, year) %>%
  mutate(gap = doy - lag(doy, default=doy[1]))
maxgap <- gap_per_year %>%
  group_by(site_name, year) %>%
  summarize_at(.vars = "gap", .funs = max)
## subset for sites with a max gap of 7 days
sub_by_gap <- maxgap[which(maxgap$gap <= 7),]
length(levels(as.factor(sub_by_gap$site_name))) ## 210
## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))
## at least 300 days per year
sub_by_gap <- sub_by_gap[which(sub_by_gap$n >= 300),]
sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),] # 59 observations

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),] ## only sites with two or more years

## Subset to years that meet criteria
sub_by_gap$site_year <- paste(sub_by_gap$site_name,sub_by_gap$year,sep = "_")
TS$site_year <- paste(TS$site_name, TS$year,sep = "_")
TS <- TS[which(TS$site_year %in% sub_by_gap$site_year),]
TS_site <- s[which(s$site_name %in% high_q$site_name),]

## Attach the median GPP
TS$GPP_temp <- TS$GPP
TS[which(TS$GPP < 0),]$GPP_temp <- sample(exp(-6):exp(-4), 1)
TS_gpp <- TS %>%
  group_by(site_name) %>%
  summarise_at(.vars = "GPP_temp", .funs = c(mean, max))
colnames(TS_gpp) <- c("site_name","GPP_mean","GPP_max")
TS_site <- left_join(TS_site, TS_gpp, by="site_name")

## Assign a stream order classification
TS_site$order_group <- "NA"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(1,2)),]$order_group <- "small"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(3,4)),]$order_group <- "mid"
TS_site[which(TS_site$NHD_STREAMORDE %in% c(5,6)),]$order_group <- "large"

###########################################################################
## Choose two consecutive river years from small, mid, and large rivers
###########################################################################

## choose sites from different groups
View(TS_site[which(TS_site$order_group == "small"),])
View(TS_site[which(TS_site$order_group == "mid"),])
View(TS_site[which(TS_site$order_group == "large"),])

## plot
sid <- "nwis_11273400"
TS_site[which(TS_site$site_name == sid),]

plot_grid(
  ggplot(TS[which(TS$site_name == sid),], aes(date, GPP_temp))+
    geom_line()+labs(title=TS_site[which(TS_site$site_name == sid),]$long_name),
  ggplot(TS[which(TS$site_name == sid & TS$year %in% c(2015,2016)),], aes(date, GPP_temp))+
  geom_line(),
  ncol = 1)

## small: nwis_05406457 2015,2016 (Order 1; BLACK EARTH CREEK NR BREWERY RD AT CROSS PLAINS,WI)
## small: nwis_01656903 2013,2014 (Order 2; FLATLICK BRANCH ABOVE FROG BRANCH AT CHANTILLY, VA)
## mid: nwis_14206950 2009,2010 (Order 3; FANNO CREEK AT DURHAM, OR)
## mid: nwis_07191222 2009,2010 (Order 3; Beaty Creek near Jay, OK)
## large: nwis_01608500 2012,2013 (Order 5; SOUTH BRANCH POTOMAC RIVER NEAR SPRINGFIELD, WV)
## large: nwis_11273400 2015,2016 (Order 6; SAN JOAQUIN R AB MERCED R NR NEWMAN CA)

site_subset <- rbind(TS[which(TS$site_name == "nwis_05406457" & TS$year %in% c(2015,2016)),],
               TS[which(TS$site_name == "nwis_01656903" & TS$year %in% c(2013,2014)),],
               TS[which(TS$site_name == "nwis_14206950" & TS$year %in% c(2009,2010)),],
               TS[which(TS$site_name == "nwis_07191222" & TS$year %in% c(2009,2010)),],
               TS[which(TS$site_name == "nwis_01608500" & TS$year %in% c(2012,2013)),],
               TS[which(TS$site_name == "nwis_11273400" & TS$year %in% c(2015,2016)),])

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),]

###################################################
## Check other covariate data quality
###################################################
site_sub_list <- split(site_subset, site_subset$site_name)

plotting_covar <- function(x) {
  
  df <- x
  p <- plot_grid(

    ggplot(df, aes(date, GPP))+
      geom_point(color="chartreuse4", size=2)+
      geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="darkolivegreen4")+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), title=df$site_name[1])+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, Q))+
      geom_line(size=1.5, color="deepskyblue4")+
      labs(y="Q (cms)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, temp))+
      geom_line(size=1.5, color="#A11F22")+
      labs(y="Water Temp (C)")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12)),
    
    ggplot(df, aes(date, light))+
      geom_point(size=2, color="darkgoldenrod3")+
      labs(y="Incoming Light", x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.text = element_text(size=12),
            axis.title = element_text(size=12)),
    ncol=1, align="hv")
  
  return(p)
  
}

plotting_covar(site_sub_list$nwis_14206950)

setwd("~/GitHub/RiverBiomass/figures/Site Covariate Plots")
lapply(site_sub_list, function(x) ggsave(plot = plotting_covar(x),filename = paste(x$site_name[1],"covar.jpg",sep = "")))


###########################
## Export
###########################

## NWIS site subset
setwd("~/GitHub/RiverBiomass/code")
saveRDS(site_subset, "./rds files/NWIS_6site_subset.rds")
saveRDS(TS_site_subset, "./rds files/NWIS_6siteinfo_subset.rds")



##############################
## Plot for talks
###########################

df <- site_sub_list$nwis_14211010
ratio_QL <- max(df$light)/max(df$Q)
GPP_plot <- ggplot(df, aes(date, GPP))+
  geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
      geom_point(color="chartreuse4", size=2)+geom_line(color="chartreuse4", size=1)+
      labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'))+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text.x = element_blank(),
            axis.text.y = element_text(size=12),
            axis.title.y = element_text(size=12))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
    
data_plot <- ggplot(df, aes(date, Q*ratio_QL))+
      geom_point(data=df, aes(date, light), size=1.5, color="darkgoldenrod3")+
      geom_line(size=1, color="deepskyblue4")+
      scale_y_continuous(sec.axis = sec_axis(~./ratio_QL, name=expression("Daily Q (cms)")))+
      labs(y=expression('Daily PPFD'))+# ('*~mu~mol~ m^-2~d^-1*')'), x="Date")+
      theme(legend.position = "none",
            panel.background = element_rect(color = "black", fill=NA, size=1),
            axis.title.x = element_blank(), axis.text = element_text(size=12),
            axis.title.y.left = element_text(size=12, color="darkgoldenrod3"),
            axis.title.y.right = element_text(size=12, color="deepskyblue4"),
            axis.text.x = element_text(angle=25, hjust = 1),
            strip.background = element_rect(fill="white", color="black"),
            strip.text = element_text(size=15))+
  geom_vline(xintercept = as.POSIXct("2013-01-01", format="%Y-%m-%d"),size=1, linetype="dashed")
    
GPP_plot + data_plot + plot_layout(ncol = 1)

ggplot(df, aes(light, GPP))+
  geom_point(size=1.5)+
  labs(y=expression('GPP (g '*~O[2]~ m^-2~d^-1*')'), x="Daily PPFD")+
  theme_bw()+
  theme(legend.position = "none",
        panel.background = element_rect(color = "black", fill=NA, size=1),
        axis.text = element_text(size=14),
        axis.title = element_text(size=18))


