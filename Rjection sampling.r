######### Rjection sampling
#Feed the model with the obs input and samples from prior of paramters
#accept the sampeles with the highest likelihood-lowest ssq or distance
####################################
library(BioCro)
rm(list=ls())
data(obsBea)
###############################
#Number of parameter values
Npara<-10000
#Sampling from the prior distribution of the parameter R
fit<-data.frame(vmax=rnorm(Npara, 39,10),alpha=rnorm(Npara, 0.04,.01),
                predict=rep(NA,Npara),obs=rep(NA,Npara),dis=rep(NA,Npara))

#Loop running the model and calculating log likelihood values
for (i in 1:Npara) {
  #run the c4photo function for the proposed vmax at this iteration
  Simul.i<-c4photo(obsBea$Qp, obsBea$temp, obsBea$rh, vmax = fit[i,1], alpha =fit[i,2])
  
  #Calculation of log likelihood - P( Y=obs | Model(theta) ) - or the distance 
  # or the ssq
  deviation<-sqrt(sum((Simul.i$Assim-obsBea$A)^2))
  
  fit[i,3:5]<-c(mean(Simul.i$Assim),mean(obsBea$A),deviation)
}

plot(fit[fit$dis<80,1:2],col="lightgrey")
points(fit[fit$dis<60,1:2],col="grey",pch=18)
points(fit[fit$dis<30,1:2],col="red",pch=8)
points(fit[fit$dis<20,1:2],col="black",pch=18)
########## getting the summries and hist of accpeted sampeles
fit2<-fit%>%filter(dis<20)
hist(fit2$vmax, breaks=30)
hist(fit2$alpha, breaks=30)
