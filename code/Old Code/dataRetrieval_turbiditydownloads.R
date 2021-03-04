## Downloading turbidity data for sites

lapply(c("plyr","dplyr","ggplot2","cowplot","lubridate",
         "dataRetrieval","tidyverse"), require, character.only=T)

###############################
## Clackamas River, OR
##############################
## Discharge and turbidity
OR_site_do <- readNWISuv(siteNumbers="14211010", parameterCd=c("00060", "63680"),
                         startDate="2010-01-01", endDate="2010-12-31")
nrow(OR_site_do)
#vis
ggplot(OR_site_do, aes(dateTime, X_63680_00000))+geom_point()
ggplot(OR_site_do, aes(dateTime, X_00060_00000))+geom_point()
#calc mean daily
OR_site_do$doy <- yday(OR_site_do$dateTime)

OR_daily_turb <- OR_site_do %>%
  mutate(date = floor_date(dateTime, unit = "day")) %>%
  group_by(date) %>%
  summarise(mean_daily_turb = mean(X_63680_00000))
#vis
ggplot(OR_daily_turb, aes(date, mean_daily_turb))+geom_point()

write.csv(OR_daily_turb, "../data/Clackamas_daily_turb.csv")

