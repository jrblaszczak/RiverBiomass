pars<-extract(PM2_DataOutput, c("r","K","s","c","B","P","pred_GPP","sig_p"))
   
simmat<-matrix(NA,length(df[,1]),length(unlist(pars$r)))
rmsemat<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars$r)){
  
  simmat[,i]<-PM2(pars$r[i],pars$K[i],pars$s[i],pars$c[i],pars$sig_p,df)
  rmsemat[i]<-sqrt(sum((simmat[,i]-df$GPP)^2)/length(df$GPP))
  
}


hist(rmsemat)
        