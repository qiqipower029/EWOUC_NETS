########## EWOC-NETS Simulation Code ########
# Toxiciy-Dose Model: logit(S_i) = beta_0+beta_1*X_i
# Assume standardized dose level from 0 to 1
library(boot)
library(truncnorm)
library(coda)
library(lattice)

sim<-100              # simulation iteration
iteration<-30000      # MCMC procedure iteration 
begin<-25000          # burn-in 
end<-30000

N <- 40          # Total sample size
fb <- 0.25       # feasibiity bound
theta<-0.476     # Target NETS  Pr(DLT)=0.33

#### Simulation Scenario Set-up ####
rho0 <- 0.05    # initial rho0
gamma <- 0.6
beta0<-logit(rho0)  #calculate beta0 for likelihood
beta1<-(logit(theta)-logit(rho0))/gamma  #calculate beta1 for likelihood

simdata <- list(S=0.05,X=10^-6,rho0 = 0.05,gamma=0.6,beta0=beta0,beta1=beta1)
w = 3 # Penalty term to calculate utility

for (i in 1:sim){
  
  if( (i!=1)&(fb<0.45) ){
    fb <- fb+0.05
  }
   
  out <- posterior(simdata$rho0[i],simdata$gamma[i],theta,simdata,iteration)
  
  estrho0 <- out[[1]][begin:end]
  estgamma <- out[[2]][begin:end]
  estbeta0 <- logit(estrho0)
  estbeta1 <- (logit(theta)-logit(estrho0))/estgamma
  
  simdata$rho0[i+1] <- median(estrho0)
  simdata$gamma[i+1] <- median(estgamma)
  simdata$beta0[i+1] <- median(estbeta0)
  simdata$beta1[i+1] <- median(estbeta1)
  simdata$X[i+1] <- quantile(estgammat,fb)
   
  x1<- beta0+beta1*simdata$X[i+1] 
  ANETS <- inv.logit(x1)
  simdata$S[i+1] <- rtruncnorm(1,a=0,b=1,mean=ANETS,sd=sqrt(ANETS*(1-ANETS)/N))
  
}
proc.time()-s
simdata
