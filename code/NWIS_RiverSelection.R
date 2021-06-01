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
colnames(site)

## Import StreamLight from Savoy
# https://www.sciencebase.gov/catalog/item/5f974adfd34e198cb77db168
SL <- read.table("../data/StreamLight_site information and parameters.txt", header=T)
colnames(SL)[colnames(SL) == "Site_ID"] <- "site_name"

## Merge
df <- left_join(site, SL, "site_name")
colnames(df)

## subset
sub <- df[,c("site_name","long_name","StreamOrde",
           "site_type","struct.canal_flag","struct.dam_flag","struct.npdes_flag")]

## Secondary stream order source from hypoxia data set
#https://www.sciencebase.gov/catalog/item/606f60afd34ef99870188ee5
hyp <- fread("../data/GRDO_GEE_HA_NHD.csv")
hyp <- hyp[which(hyp$DB_Source == "PC"), c("SiteID","ORD_STRA","NHD_STREAMORDE")]
colnames(hyp)[which(colnames(hyp) == "SiteID")] <- "site_name"

#Merge
sub <- left_join(sub, hyp, by="site_name")

##################################################
## Select sites based on data quality
##################################################

# Import and subset model diagnostics
diagnostics <- read.table("../data/diagnostics.tsv",sep = "\t", header=T)
diagnostics <- diagnostics[which(diagnostics$site %in% sub$site_name),]
highq_sites <- diagnostics[which(diagnostics$K600_daily_sigma_Rhat < 1.05 & 
                                   diagnostics$err_obs_iid_sigma_Rhat < 1.05 &
                                   diagnostics$err_proc_iid_sigma_Rhat < 1.05 &
                                   diagnostics$neg_GPP < 15 & diagnostics$pos_ER < 15),] #229
highq_site_names <- unique(highq_sites$site) ## 208


# Subset s based on high sites and site type and flags
s <- sub[which(sub$site_name %in% highq_site_names),] ## 208
s <- s[which(s$site_type == "ST"),] ## 204
s <- s[which(s$struct.dam_flag %in% c(NA,"95")),] ## 97
# which have light from Phil
#s_l <- s[!is.na(s$StreamOrde),] # 41

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

## Subset to sites in high_sites (sites with high confidence rating and limited dam interference)
NWIS_sub <- NWIS_sub[which(NWIS_sub$site_name %in% s$site_name),]
# Confirm
length(levels(as.factor(NWIS_sub$site_name))) ## 97 when subsetting for s

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
## subset for sites with a max gap of 14 days
sub_by_gap <- maxgap[which(maxgap$gap <= 14),]
length(levels(as.factor(sub_by_gap$site_name))) ## 92
## merge with number of days per year
sub_by_gap <- merge(sub_by_gap, dat_per_year, by=c("site_name","year"))
## at least 275 days per year
sub_by_gap <- sub_by_gap[which(sub_by_gap$n >= 275),]
sub_by_gap_sum <- sub_by_gap %>% group_by(site_name) %>% count()
high_q <- sub_by_gap_sum[which(sub_by_gap_sum$n >= 2),] # 41

## Subset NWIS_sub
TS <- NWIS_sub[which(NWIS_sub$site_name %in% high_q$site_name),] ## only sites with two or more years

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
TS_site[which(TS_site$NHD_STREAMORDE %in% c(3,4,5)),]$order_group <- "mid"
TS_site[which(TS_site$NHD_STREAMORDE >= 6),]$order_group <- "large"

###########################################################################
## Choose two consecutive river years from small, mid, and large rivers
###########################################################################

## choose sites from different groups
View(TS_site[which(TS_site$order_group == "small"),])
View(TS_site[which(TS_site$order_group == "mid"),])
View(TS_site[which(TS_site$order_group == "large"),])

## plot
sid <- "nwis_08447300"
two_years <- c(2012,2013)
TS_site[which(TS_site$site_name == sid),]

plot_grid(
  ggplot(TS[which(TS$site_name == sid),], aes(date, GPP_temp))+
    geom_line()+labs(title=TS_site[which(TS_site$site_name == sid),]$long_name),
  ggplot(TS[which(TS$site_name == sid & TS$year %in% two_years),], aes(date, GPP_temp))+
  geom_line(),
  ncol = 1)

## small: nwis_02336526 2015,2016 (Order 2; PROCTOR CREEK AT JACKSON PARKWAY, AT ATLANTA, GA) - light
## small: nwis_01649190 2010,2011 (Order 2; PAINT BRANCH NEAR COLLEGE PARK, MD) - light
## mid: nwis_07191222 2009,2010 (Order 3; Beaty Creek near Jay, OK) - light
## mid: nwis_01608500 2012,2013 (Order 5; SOUTH BRANCH POTOMAC RIVER NEAR SPRINGFIELD, WV) - light
## large: nwis_11044000 2015,2016 (Order 6; SANTA MARGARITA R NR TEMECULA CA) - no light
## large: nwis_08447300 2012,2013 (Order 7: Pecos Rv at Brotherton Rh nr Pandale, TX) - no light

site_subset <- rbind(TS[which(TS$site_name == "nwis_02336526" & TS$year %in% c(2015,2016)),],
               TS[which(TS$site_name == "nwis_01649190" & TS$year %in% c(2010,2011)),],
               TS[which(TS$site_name == "nwis_07191222" & TS$year %in% c(2009,2010)),],
               TS[which(TS$site_name == "nwis_01608500" & TS$year %in% c(2012,2013)),],
               TS[which(TS$site_name == "nwis_11044000" & TS$year %in% c(2015,2016)),],
               TS[which(TS$site_name == "nwis_08447300" & TS$year %in% c(2012,2013)),])

TS_site_subset <- df[which(df$site_name %in% site_subset$site_name),]

###########################
## Export
###########################

## NWIS site subset
setwd("~/GitHub/RiverBiomass/code")
saveRDS(site_subset, "./rds files/NWIS_6site_subset_SL.rds")
saveRDS(TS_site_subset, "./rds files/NWIS_6siteinfo_subset_SL.rds")



##############################
## Plot for talks
###########################
site_subset_list <- split(site_subset, site_subset$site_name)


df <- site_subset_list$nwis_01649190
ratio_QL <- max(df$light)/max(df$Q)
GPP_plot <- ggplot(df, aes(date, GPP_temp))+
  #geom_errorbar(aes(ymin = GPP.lower, ymax = GPP.upper), width=0.2,color="chartreuse4")+
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


