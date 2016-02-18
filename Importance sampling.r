######### Importance sampling
#Feed the model with the obs input and samples from prior of paramters
# find a weight for each set of vmax and alpha sample. Samples with higher weight (higher
#likelihood) is  taken
####################################
library(BioCro)
library(mc2d)
library(ggplot2)
library(dplyr)
data(obsBea)
############################
#Number of parameter values
Npara<-10000
#Sampling from the prior distribution of the parameter R
fit<-data.frame(vmax=rnorm(Npara, 39,10),alpha=rnorm(Npara, 0.04,.01),
                predict=rep(NA,Npara),obs=rep(NA,Npara),weight=rep(NA,Npara))
#Sampling from the prior distribution of the residual model standard deviation
Ssample<-runif(Npara,0,20)
#Initialisation of the vector of log likelihood
LogLike_vec<-rep(NA,Npara)
#Loop running the model and calculating log likelihood values
for (i in 1:Npara) {
  #run the c4photo function for the proposed vmax at this iteration
  Simul.i<-c4photo(obsBea$Qp, obsBea$temp, obsBea$rh, vmax = fit[i,1], alpha =fit[i,2])
  #Calculation of log likelihood - P( Y=obs | Model(theta) )
  # sd is unknown as well
  LogLikelihood.i<-sum(log(dnorm(obsBea$A,Simul.i$Assim,Ssample[i])))
  LogLike_vec[i]<-LogLikelihood.i
  fit[i,3:4]<-c(mean(Simul.i$Assim),mean(obsBea$A))
}
#Weight calculation
Weight<-exp(LogLike_vec)/sum(exp(LogLike_vec))
fit$weight<-Weight

fit%>%ggplot()+
  geom_point(aes(x=vmax,y=alpha,color=weight))+
  scale_colour_gradient(low = "grey", high = "black")

plot(fit$vmax,fit$weight)
plot(fit$alpha,fit$weight)
#### Expected value of Vmax
sum(Weight*fit$vmax)
#### Expected value of Vmax
sum(Weight*fit$alpha)
