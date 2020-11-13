## Model-specific functions for creating figures

## Confidence interval extraction

CI_par_PM1 <- function(par) {
  ## Find CI of all
  medCI_par <- lapply(par, function(x) quantile(x, probs = 0.5))
  lowerCI_par <- lapply(par, function(x) quantile(x, probs = 0.025))
  upperCI_par <- lapply(par, function(x) quantile(x, probs = 0.975))
  ## CI of ts parameters
  medCI_l_pred_GPP_ts <- apply(par$l_pred_GPP,2,function(x) quantile(x, probs = 0.5))
  lowerCI_l_pred_GPP_ts <- apply(par$l_pred_GPP,2,function(x) quantile(x, probs = 0.025))
  upperCI_l_pred_GPP_ts <- apply(par$l_pred_GPP,2,function(x) quantile(x, probs = 0.975))
  ## Compile in list and return
  CI_par_l <- list(medCI_par, medCI_l_pred_GPP_ts,lowerCI_par, lowerCI_l_pred_GPP_ts,upperCI_par, upperCI_l_pred_GPP_ts)
  names(CI_par_l) <- c("medCI_par", "medCI_l_pred_GPP","lowCI_par","lowCI_l_pred_GPP","upCI_par","upCI_l_pred_GPP")
  return(CI_par_l)
}

CI_par_PM2 <- function(par) {
  ## Find CI of all
  medCI_par <- lapply(par, function(x) quantile(x, probs = 0.5))
  lowerCI_par <- lapply(par, function(x) quantile(x, probs = 0.025))
  upperCI_par <- lapply(par, function(x) quantile(x, probs = 0.975))
  ## CI of ts parameters
  medCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.5))
  lowerCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.025))
  upperCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.975))
  ## Compile in list and return
  CI_par_l <- list(lowerCI_par,medCI_par,upperCI_par,lowerCI_pred_GPP_ts,medCI_pred_GPP_ts,upperCI_pred_GPP_ts)
  names(CI_par_l) <- c("lowCI_par","medCI_par","upCI_par","lowCI_pred_GPP","medCI_pred_GPP","upCI_pred_GPP")
  return(CI_par_l)
}

CI_par_PM3 <- function(par) {
  ## Find CI of all
  medCI_par <- lapply(par, function(x) quantile(x, probs = 0.5))
  lowerCI_par <- lapply(par, function(x) quantile(x, probs = 0.025))
  upperCI_par <- lapply(par, function(x) quantile(x, probs = 0.975))
  ## CI of ts parameters
  medCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.5))
  lowerCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.025))
  upperCI_pred_GPP_ts <- apply(par$pred_GPP,2,function(x) quantile(x, probs = 0.975))
  ## Compile in list and return
  CI_par_l <- list(lowerCI_par,medCI_par,upperCI_par,lowerCI_pred_GPP_ts,medCI_pred_GPP_ts,upperCI_pred_GPP_ts)
  names(CI_par_l) <- c("lowCI_par","medCI_par","upCI_par","lowCI_pred_GPP","medCI_pred_GPP","upCI_pred_GPP")
  return(CI_par_l)
}

## Plot comparison of input versus output parameters
require(RColorBrewer)
ygb <- brewer.pal(10, "Paired")

in_vs_out <- function(plist, inputvaluestring){
  
  PM_inputs <- as.data.frame(as.matrix(inputvaluestring))
  PM_inputs$Param <- rownames(PM_inputs)
  colnames(PM_inputs)[1] <- "Value"
  
  lowCI <- ldply(plist$lowCI_par, data.frame)
  medCI <- ldply(plist$medCI_par, data.frame)
  upCI <- ldply(plist$upCI_par, data.frame)
  CI_list <- list(lowCI, medCI, upCI)
  CI_list <- lapply(CI_list, setNames, nm = c("Param","Value"))
  
  CI_df <- bind_cols(CI_list)
  CI_df <- CI_df[,c("Param","Value","Value1","Value2")]
  colnames(CI_df) <- c("Param","ParOutput_low","ParOutput_med","ParOutput_up")
  
  in_out <- merge(PM_inputs, CI_df, by="Param")
  
  ggplot(in_out, aes(log(abs(Value)), log(abs(ParOutput_med)), color=Param))+
    geom_abline(slope = 1, intercept = 0, color="grey", size=1) +
    geom_linerange(aes(ymin = log(abs(ParOutput_low)), ymax= log(abs(ParOutput_up))), color="black",size=1)+
    geom_point(size=4)+
    labs(x = "Input Parameter Values", y="Median Posterior Values")+
    theme(axis.text.x = element_text(angle = 45, hjust=1),
          axis.text = element_text(size=15),
          axis.title = element_text(size=20),
          legend.text = element_text(hjust = 0, size=15))+
    scale_color_manual(name="",values=c('alpha'=ygb[1], 'beta_r'=ygb[2],'bQ'=ygb[3],
    'critQ'=ygb[4], 'gamma'=ygb[5],'r'=ygb[6],'sig_p'=ygb[7],'beta_0'=ygb[8],
    'phi'=ygb[9],'beta'=ygb[10]), labels=c('alpha'=expression(alpha), 'beta_r'=expression(beta[r]),'bQ'="s",
                                'critQ'="critQ", 'gamma'=expression(gamma),'r'="r",'sig_p'=expression(sigma[p]),
                                'beta_0'=expression(theta),'phi'=expression(phi),'beta'=expression(beta[l])))+
    scale_x_continuous(limits=c(-6,6))+
    scale_y_continuous(limits=c(-6,6))
  
}

