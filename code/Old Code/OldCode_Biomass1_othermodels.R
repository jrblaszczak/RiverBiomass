#################
## Other Models
#################

## PM 2 - Latent Biomass (Ricker) - time varying r
init_Ricker <- function(...) {
  list(c = 0.5, s = 100)
}

PM_outputlist_Ricker_tvr <- lapply(stan_data_l,
                                   function(x) stan("Stan_ProductivityModel2_Ricker_fixedinit_r.stan",
                                                    data=x,chains=3,iter=5000,init = init_Ricker,
                                                    control=list(max_treedepth=12)))
PM_Ricker_elapsedtime_tvr <- lapply(PM_outputlist_Ricker_tvr, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Ricker_tvr, "./rds files/stan_6riv_output_Ricker_tvr_2021_06_01.rds")
saveRDS(PM_Ricker_elapsedtime_tvr, "./rds files/stan_6riv_Ricker_tvr_time_2021_06_01.rds")





## PM 3 - Latent Biomass (Gompertz)
init_Gompertz <- function(...) {
  list(c = 0.5, s = 200)
}

PM_outputlist_Gompertz <- lapply(stan_data_l,
                                 function(x) stan("Stan_ProductivityModel3_Gompertz_fixedinit_obserr.stan",
                                                  data=x,chains=3,iter=5000, 
                                                  control=list(max_treedepth=12)))
PM_Gompertz_elapsedtime <- lapply(PM_outputlist_Gompertz, function(x) return(get_elapsed_time(x)))
saveRDS(PM_outputlist_Gompertz, "./rds files/stan_6riv_output_Gompertz_2021_06_01.rds")
saveRDS(PM_Gompertz_elapsedtime, "./rds files/stan_6riv_Gompertz_time_2021_06_01.rds")


## View
launch_shinystan(PM_outputlist_AR$nwis_08447300)
launch_shinystan(PM_outputlist_Ricker$nwis_07191222)
launch_shinystan(PM_outputlist_Ricker_tvr$nwis_11044000)

PM_outputlist_AR$nwis_01608500$lp__
