
logit<-function(p){return(log(p/(1-p)))}
#---------------------------------------- Scenario 1 ----------------------------------------#
gammat = 1.0
gammae = 0.4
sim = 1000
tau=0.1 ######correlation between toxicity and efficacy######  
##### sencario setting######
scenario=sname     ##scenario name###
cohort=3           ##cohort size###
thetat<-0.476      ###target tolerated level###
thetae = 0.3
upe<-rep(NA,sim)    ####Utility###
res<-c(NA,NA,NA)     ####initial result####    
######model####
Xmin=0.2;Xmax=1;rhot=0.01;rhoe=0.02  #####initial value
b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
xi = exp(b0t + b1t*std.dose)/(1+exp(b0t + b1t*std.dose))
mu = exp(b0e + b1e*std.dose)/(1+exp(b0e + b1e*std.dose))
lambda = 1
sig2 = std.dose^lambda
w = 3
s1.pe = mu
s1.pt = xi
d.u<-s1.pe-w*s1.pt       ####true utility###
real.table=rbind(s1.pt,s1.pe,d.u)  ####combine true  toxicity efficacy and utility###


#---------------------------------------- Scenario 2 ----------------------------------------#
gammat = 0.8
gammae = 0.4
sim = 1000
tau=0.1 ######correlation between toxicity and efficacy######  
##### sencario setting######
scenario=sname     ##scenario name###
cohort=3           ##cohort size###
thetat<-0.476      ###target tolerated level###
thetae = 0.333
upe<-rep(NA,sim)    ####Utility###
res<-c(NA,NA,NA)     ####initial result####    
######model####
Xmin=0.2;Xmax=1;rhot=0.01;rhoe=0.02  #####initial value
b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
xi = exp(b0t + b1t*std.dose)/(1+exp(b0t + b1t*std.dose))
mu = exp(b0e + b1e*std.dose)/(1+exp(b0e + b1e*std.dose))
lambda = 1
sig2 = std.dose^lambda
w = 3
s1.pe = mu
s1.pt = xi
d.u<-s1.pe-w*s1.pt       ####true utility###
real.table=rbind(s1.pt,s1.pe,d.u)  ####combine true  toxicity efficacy and utility###
