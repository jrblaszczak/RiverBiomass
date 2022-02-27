########################################
## Extract Stan Parameter Estimates 
#######################################

## For S-TS

STS_extract_medians <- function(PM_par){
  
  ## Summarize parameters
  med_par_PM <- function(par) {
    
    ## Find the med of all
    med_par <- lapply(par, function(x) median(x))
    
    ## med of ts parameters
    med_pred_GPP_ts <- apply(exp(par$l_pred_GPP),2,median)
    
    ## sd of ts parameters
    Q.025_pred_GPP_ts <- apply(exp(par$l_pred_GPP),2, function(x) quantile(x, probs = 0.025))
    Q.975_pred_GPP_ts <- apply(exp(par$l_pred_GPP),2, function(x) quantile(x, probs = 0.975))
    
    ## Compile in list and return
    med_par_ts <- list(med_par, med_pred_GPP_ts,
                       Q.025_pred_GPP_ts,Q.975_pred_GPP_ts)
    names(med_par_ts) <- c("par","pred_GPP","pred_GPP_Q.025","pred_GPP_Q.975")
    return(med_par_ts)
  }
  
  PM_medpar <- med_par_PM(PM_par)
}


## For LB-TS model

LBTS_extract_medians <- function(PM_par){
  
  ## Summarize parameters
  med_par_PM <- function(par) {
    
    ## Find the med of all
    med_par <- lapply(par, function(x) median(x))
    
    ## med of ts parameters
    med_pred_GPP_ts <- apply(par$pred_GPP,2,median)
    med_B_ts <- apply(par$B,2,median)
    med_P_ts <- apply(par$P,2,median)
    
    ## sd of ts parameters
    Q.025_pred_GPP_ts <- apply(par$pred_GPP,2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    Q.975_pred_GPP_ts <- apply(par$pred_GPP,2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
    Q.025_B_ts <- apply(par$B,2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    Q.975_B_ts <- apply(par$B,2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
    Q.025_P_ts <- apply(par$P,2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
    Q.975_P_ts <- apply(par$P,2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
    
    ## Compile in list and return
    med_par_ts <- list(med_par, med_pred_GPP_ts,med_B_ts, med_P_ts,
                     Q.025_pred_GPP_ts,Q.975_pred_GPP_ts,
                     Q.025_B_ts, Q.975_B_ts,
                     Q.025_P_ts, Q.975_P_ts)
  names(med_par_ts) <- c("par","pred_GPP","B","P","pred_GPP_Q.025","pred_GPP_Q.975",
                         "B_Q.025","B_Q.975","P_Q.025","P_Q.975")
  return(med_par_ts)
}

  PM_medpar <- med_par_PM(PM_par)
}

## For mechanistic model 2 with latent biomass
# random parameter set to simulate data

LBTS_randpar <- function(PM_par){
  
  ## Summarize parameters
  rand_par_PM <- function(par) {
    i <- sample(size=1, x=1:length(par))
    randpar <- lapply(par, function(x) return(x[i]))
    return(randpar)
  }
  
  PM_randpar <- rand_par_PM(PM_par)
  return(PM_randpar)
}




