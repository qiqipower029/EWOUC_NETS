##### Function for EWOC-NETS Simulation ######

likelihood<-function(rho0_p,gamma_p,theta,simdata){
  
  beta0<-logit(rho0_p)   #calculate beta0 for likelihood
  
  beta1<-(logit(theta)-logit(rho0_p))/gamma_p  #calculate beta1 for likelihood
  
  s<-simdata$S    #nets of recent patients
  
  prob<-exp(beta0+beta1*simdata$X)/(1+exp(beta0+beta1*simdata$X))
  
  likelihood<-prod((prob^s)*(1-prob)^(1-s))   #quasi likelihood of all recent patients
  
  return(likelihood)
}

### Update rho0 ###
rho0_up<-function(prerho0,posgamma,theta,simdata){
  rho0_new<-runif(1,0,theta) #randomly select new value for rho0
  
  lik_rho_new<-likelihood(rho0_new,posgamma,theta,simdata=simdata)
  lik_rho_old<-likelihood(prerho0,posgamma,theta,simdata=simdata)
  
  ratio_rho<-lik_rho_new/lik_rho_old
  
  if(ratio_rho=="NaN"){ratio_rho<-1}
  
  if(ratio_rho>=1){rho_next<-rho0_new}
  else{test<-runif(1)
  if(test<ratio_rho){rho_next<-rho0_new}
  else{rho_next<-prerho0}
  
  }
  return(rho_next)
}


# Update gamma #
gamma_up<-function(posrho0,pregamma,theta,simdata){
  gamma_new<-runif(1,0,1)
  
  lik_gamma_new<-likelihood(posrho0,gamma_new,theta,simdata=simdata)
  lik_gamma_old<-likelihood(posrho0,pregamma,theta,simdata=simdata)
  ratio_gamma<-lik_gamma_new/lik_gamma_old
  
  if(ratio_gamma=="NaN"){ratio_gamma<-1}
  if(ratio_gamma>=1){gamma_next<-gamma_new}
  else{
    test<-runif(1)
    if(test<ratio_gamma){gamma_next<-gamma_new}
    else{gamma_next<-pregamma}
  }
  return(gamma_next)
}


posterior<-function(orirho0,origamma,theta,simdata,int){
  # create variables to store updated rho0 and gamma values.
  rho0_s<-c(orirho0,rep(NA,int-1))
  gamma_s<-c(origamma,rep(NA,int-1))
  
  #First update rho0 based on input data and assume gamma is known 
  for(i1 in 2:int){
    rho0_s[i1]<-rho0_up(rho0_s[i1-1],origamma,theta,simdata=simdata)
  }
  
  #Second update gamma based on input data and 
  #assume rho0 (bayesian estimator) known
  b <- median(rho0_s[25000:30000])
  for(i2 in 2:int){
    gamma_s[i2]<-gamma_up(b,gamma_s[i2-1],theta,simdata=simdata)
  }
  
  posterior<-list(rho0_s,gamma_s)
  return(posterior)
  
}




