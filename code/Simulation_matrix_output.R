##################
## Model 1 Output
##################

pars1<-extract(PM1_DataOutput, c("phi","alpha","beta","sig_p"))

simmat1<-matrix(NA,length(df[,1]),length(unlist(pars1$r)))
rmsemat1<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars1$r)){
  
  simmat1[,i]<-PM1(pars1$phi[i],pars1$alpha[i],pars1$beta[i],pars1$sig_p,df)
  rmsemat1[i]<-sqrt(sum((simmat1[,i]-df$GPP)^2)/length(df$GPP))
  
}

hist(rmsemat1)

##################
## Model 2 Output
##################

pars2<-extract(PM2_DataOutput, c("r","K","s","c","B","P","pred_GPP","sig_p"))
   
simmat2<-matrix(NA,length(df[,1]),length(unlist(pars2$r)))
rmsemat2<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars2$r)){
  
  simmat2[,i]<-PM2(pars2$r[i],pars2$K[i],pars2$s[i],pars2$c[i],pars2$sig_p,df)
  rmsemat2[i]<-sqrt(sum((simmat2[,i]-df$GPP)^2)/length(df$GPP))
  
}

hist(rmsemat2)

##################
## Model 3 Output
##################

pars3<-extract(PM3_DataOutput, c("alpha","gamma","s","c","N","P","pred_GPP","sig_p"))

simmat3<-matrix(NA,length(df[,1]),length(unlist(pars3$r)))
rmsemat3<-matrix(NA,length(df[,1]),1)

for (i in 1:length(pars3$r)){
  
  simmat3[,i]<-PM3(pars3$alpha[i],pars3$gamma[i],pars3$s[i],pars3$c[i],pars3$sig_p,df)
  rmsemat3[i]<-sqrt(sum((simmat3[,i]-df$GPP)^2)/length(df$GPP))
  
}

hist(rmsemat3)


cbind(rmsemat1, rmsemat2, rmsemat3)







        