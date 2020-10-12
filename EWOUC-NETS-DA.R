############EWOUC-NETS########
rm(list=ls())
load("C:/Users/Shiru/Desktop/EWOUC-NETS-DA/potential results/potent.txt")
############settings###########
dose_level=c(20,25,30,35,40,45)
dose_level=dose_level/45
Xmin=20/45
Xmax=45/45
N=12
Time=3
h1=1
h2=2
h3=3
reallambda=c(0.4,0.667,2)
####scenarios#########
path<-"C:/Users/Shiru/Desktop/?о?/scenario - NEST"
fileNames <- dir(path)  ##??ȡ??·???µ??ļ???
filePath <- sapply(fileNames, function(x){ 
                 paste(path,x,sep='/')})   ##???ɶ?ȡ?ļ?·??
data<- lapply(filePath, function(x){
             read.table(x)})  

###functions####
get_dose=function(x){
 q=NULL
 q=dose_level[x]
 return(q)
}
get_doselevel=function(x){
 q=NULL
 if(!is.na(x)){
 for(i in 1:length(dose_level)){
  if(dose_level[i]==x)q=i
 }}
 else q=NA
 return(q)
}
choose_dose1=function(m1,m2){
 q=NULL
 j=1
 if(m1<m2)return(NA)
 if(m1>=m2){
  for(i in 1:length(dose_level)){
   if(dose_level[i]>=m2&&dose_level[i]<=m1){
    q[j]=dose_level[i]
    j=j+1
   }
  }
  if(is.null(q))return(NA)
 }
 return(max(q))
}
choose_dose2=function(m1){
 q=NULL
 j=1
for(i in 1:length(dose_level)){
  if(dose_level[i]<=m1){
   q[j]=dose_level[i]
   j=j+1
  }
  if(is.null(q))return(NA)
 }
 return(max(q))
}
getprob=function(x,K,lambda,gammat,gammae,rhot,thetat,rhoe,thetae,phi){
 q=NULL
 q1=sum(lambda)
 beta0t=(gammat*logit(rhot)-Xmin*logit(thetat))/(gammat-Xmin)
 beta1t=(logit(thetat)-logit(rhot))/(gammat-Xmin)
 beta0e=(gammae*logit(rhoe)-Xmin*logit(thetae))/(gammae-Xmin)
 beta1e=(logit(thetae)-logit(rhoe))/(gammae-Xmin)
 q2=f1(x,beta0t,beta1t,beta0e,beta1e,phi)
 q3=pit(x,beta0t,beta1t)-q2
 q=q2*exp(-q1)/(q2*exp(-q1)+1-q3)
 return(q)
}
generate_Y=function(x,y,time,lambda,gammat,gammae,rhot,thetat,rhoe,thetae,phi){
 s=y
 K=get_K(time)
 s[s==0]=rbinom(length(s[s==0]),1,getprob(x,K,lambda,gammat,gammae,rhot,thetat,rhoe,thetae,phi))
 return(s)
}
get_K=function(time){
 s=matrix(0,ncol=3,nrow=length(time))
 for(i in 1:length(time)){
  if(time[i]>h1){
   s[i,1]=h1
   if(time[i]>h2){
    s[i,2]=h2-h1
    if(time[i]>h3){
     s[i,3]=h3-h2
    }
    else s[i,3]=time[i]-h2
   }
   else s[i,2]=time[i]-h1
  }
  else s[i,1]=time[i]
 }
 return(s)
}
get_sigma=function(time){
 s=matrix(0,ncol=3,nrow=length(time))
 for(i in 1:3){
  if(time[i]<=h1)s[i,1]=1
  if(time[i]<=h2&&time[i]>h1)s[i,2]=1
  if(time[i]>h2)s[i,3]=1
 }
 return(s)
}
logit<-function(p){return(log(p/(1-p)))}

####likelihood functions##########
pit=function(x,beta0t,beta1t) 1/(1+exp(-(beta0t+beta1t*x)))
pie=function(x,beta0e,beta1e) 1/(1+exp(-(beta0e+beta1e*x)))
f1=function(x,beta0t,beta1t,beta0e,beta1e,phi) pit(x,beta0t,beta1t)*pie(x,beta0e,beta1e)*(1+(exp(phi)-1)/(exp(phi)+1)*(1-pit(x,beta0t,beta1t))*(1-pie(x,beta0e,beta1e))) # pi11
likelihood=function(x,s,y,gammat,gammae,rhot,thetat,rhoe,thetae,phi){
 q=1
 beta0e=(gammae*logit(rhoe)-Xmin*logit(thetae))/(gammae-Xmin)
 beta1e=(logit(thetae)-logit(rhoe))/(gammae-Xmin)
 beta0t=(gammat*logit(rhot)-Xmin*logit(thetat))/(gammat-Xmin)
 beta1t=(logit(thetat)-logit(rhot))/(gammat-Xmin)
 for(j in 1:length(x)){
  z1=f1(x[j],beta0t,beta1t,beta0e,beta1e,phi)
  z2=pit(x[j],beta0t,beta1t)
  z3=pie(x[j],beta0e,beta1e)
  for(i in 1:3){
   q=q*(z1^(y[j,i]*s[j,i]))*(z2-z1)^((1-y[j,i])*s[j,i])*(z3-z1)^(y[j,i]*(1-s[j,i]))*(1-z2-z3+z1)^((1-y[j,i])*(1-s[j,i]))  
  }
 }
 return(q)
} # return the likelihood
likely=function(time,y,lambda){
 p=1
 for(i in 1:nrow(time)){
  sigma=get_sigma(time[i,])
  K=get_K(time[i,])
  for(j in 1:3){
   for(h in 1:3){
    p=p*lambda[h]^sigma[j,h]*exp(-y[i,j]*lambda[h]*sigma[j,h])
   }
  }
 }
return(p)
}
########mcmc##########
mcmc1=function(S,D,Y,gammat0,gammae0,rhot0,thetat,rhoe0,thetae,phi0,n){
  gammat=c(gammat0,NULL)
  gammae=c(gammae0,NULL)
  rhot=c(rhot0,NULL)
  rhoe=c(rhoe0,NULL)
  phi=c(phi0,NULL)
  for(i in 2:n){
   z=0
   z1=0
   while(z==0){
    phi[i]=runif(1,min=0,max=1)
    u=runif(1,0,1)
    # compare likelihood between phi i and phi i-1, compare the ratio with 1.
    a=min(1,likelihood(D,S,Y,gammat[i-1],gammae[i-1],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i])/likelihood(D,S,Y,gammat[i-1],gammae[i-1],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i-1]))
    if(!(is.na(a))){ 
     if(u<=a) z=1 # If a is not NA, compare a with u. if a > u, then z = 1. phi[i] is the new one.
    }
    if(z1==10){ # when the number of iteration is 11, phi[i] = phi[i-1]?
     phi[i]=phi[i-1]
     z=1 # What does z and z1 mean here? Guess: z means whether phi[i] has a value, z1 means the number of iterations.
    }
    z1=z1+1 # Calculate the number of iterations in the while loop.
   }
   z=0
   z1=0
   while(z==0){
    gammat[i]=runif(1,min=Xmin,max=Xmax)
    u=runif(1,0,1)
    a=min(1,likelihood(D,S,Y,gammat[i],gammae[i-1],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i])/likelihood(D,S,Y,gammat[i-1],gammae[i-1],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i]))
    if(!(is.na(a))){ 
     if(u<=a) z=1
    }
 if(z1==10){
     gammat[i]=gammat[i-1]
     z=1
    }
    z1=z1+1

   }
   z=0
   z1=0
   while(z==0){
    gammae[i]=runif(1,min=Xmin,max=Xmax)
    u=runif(1,0,1)
    a=min(1,likelihood(D,S,Y,gammat[i],gammae[i],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i])/likelihood(D,S,Y,gammat[i],gammae[i-1],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i]))
    if(!(is.na(a))){ 
     if(u<=a) z=1
    }
     if(z1==10){
     gammae[i]=gammae[i-1]
     z=1
    }
    z1=z1+1
   }
   z=0
   z1=0
   while(z==0){
    rhot[i]=runif(1,min=0,max=thetat)
    u=runif(1,0,1)
    a=min(1,likelihood(D,S,Y,gammat[i],gammae[i],rhot[i],thetat,rhoe[i-1],thetae,phi[i])/likelihood(D,S,Y,gammat[i],gammae[i],rhot[i-1],thetat,rhoe[i-1],thetae,phi[i]))
    if(!(is.na(a))){ 
     if(u<=a) z=1
    }
    if(z1==10){
     rhot[i]=rhot[i-1]
     z=1
    }
    z1=z1+1
   }
   z=0
   z1=0
   while(z==0){
    rhoe[i]=runif(1,min=0,max=thetae)
    u=runif(1,0,1)
    a=min(1,likelihood(D,S,Y,gammat[i],gammae[i],rhot[i],thetat,rhoe[i],thetae,phi[i])/likelihood(D,S,Y,gammat[i],gammae[i],rhot[i],thetat,rhoe[i-1],thetae,phi[i]))
    if(!(is.na(a))){ 
     if(u<=a) z=1
    }
     if(z1==10){
     rhoe[i]=rhoe[i-1]
     z=1
    }
    z1=z1+1
   }
   }
  return(cbind(gammat,gammae,rhot,rhoe,phi))
}

mcmc2=function(lambda,time,Y,n){
 lambda1=c(lambda[1],NULL)
 lambda2=c(lambda[2],NULL)
 lambda3=c(lambda[3],NULL)
 for(i in 2:n){ 
  z=0
  z1=0
  while(z==0){
   lambda1[i]= rgamma(1,0.2,2)
   u=runif(1,0,1)
   lambda_1=c(lambda1[i], lambda2[i-1], lambda3[i-1])
   lambda_2=c(lambda1[i-1], lambda2[i-1], lambda3[i-1])
   a=min(1,likely(time,Y,lambda_1)/likely(time,Y,lambda_2))
   if(!(is.na(a))){ 
     if(u<=a) z=1
   }
    if(z1==10){
     lambda1[i]=lambda1[i-1]
     z=1
    }
    z1=z1+1
  }
 z=0
 z1=0
 while(z==0){
   lambda2[i]= rgamma(1,0.33,2)
   u=runif(1,0,1)
   lambda_1=c(lambda1[i], lambda2[i], lambda3[i-1])
   lambda_2=c(lambda1[i], lambda2[i-1], lambda3[i-1])
   a=min(1,likely(time,Y,lambda_1)/likely(time,Y,lambda_2))
   if(!(is.na(a))){ 
     if(u<=a) z=1
   }
   if(z1==10){
     lambda2[i]=lambda2[i-1]
     z=1
    }
    z1=z1+1
  }
 z=0
 z1=0
 while(z==0){
   lambda3[i]= rgamma(1,1,2)
   u=runif(1,0,1)
   lambda_1=c(lambda1[i], lambda2[i], lambda3[i])
   lambda_2=c(lambda1[i], lambda2[i], lambda3[i-1])
   a=min(1,likely(time,Y,lambda_1)/likely(time,Y,lambda_2))
   if(!(is.na(a))){ 
     if(u<=a) z=1
   }
    if(z1==10){
     lambda3[i]=lambda3[i-1]
     z=1
    }
    z1=z1+1
  }
  }
 return(cbind(lambda1,lambda2,lambda3))
}
###for each scenarios####
R=list()
toxi=matrix(0,ncol=7,nrow=21)
toxi[1,]=c(0,0.092,0.25,0.417,0.583,0.75,0.917)
Does=matrix(0,ncol=13,nrow=20)
prob_=matrix(0,ncol=2,nrow=20)
gammat_=matrix(0,ncol=7,nrow=20)
gammae_=matrix(0,ncol=7,nrow=20)
for(sce in 1:2){
for(loop in 1:20){
A=data[[sce]]
#A[i,j]means the probability of ith grade toxicity when treated at the jth dose
#in scenatio 1
#the last nrow: the first means MTD, the second means MED
theta1=A[9,1]#TNETS
theta2=A[9,2]#P(Y=1|Dose=MED)
time=matrix(0,ncol=3)
D=NULL#dose
X=NULL#dose level
S=matrix(0,ncol=3)#3 patient one time
Y=matrix(0,ncol=3)
real_Y=matrix(0,ncol=3)
observe_Y=matrix(0,ncol=3)
lambda=NULL
D[1]=dose_level[1]
X[1]=get_doselevel(D[1])
S[1,]=potent[[loop]][[1]][(sce*6-6+X[1]),1:3]
real_Y[1,]=potent[[loop]][[1]][(sce*6-6+X[1]),7:9]
time[1,]=potent[[loop]][[1]][(sce*6-6+X[1]),10:12]
observe_Y[1,]=potent[[loop]][[1]][(sce*6-6+X[1]),4:6]
gammat=runif(1,min=Xmin,max=Xmax)
gammae=runif(1,min=Xmin,max=Xmax)
rhot=runif(1,min=0,max=theta1)
rhoe=runif(1,min=0,max=theta2)
phi=runif(1,min=0,max=1)
lambda[1]=rgamma(1,0.2,2)
lambda[2]=rgamma(1,0.33,2)
lambda[3]=rgamma(1,1,2)
Y[1,]=generate_Y(X[1],observe_Y[1,],time[1,],lambda,gammat,gammae,rhot,theta1,rhoe,theta2,phi)
fb=0.25
begin=15000
end=20000
for(i in 2:(N+1)){
 answer1=mcmc1(S,D,Y,gammat,gammae,rhot,theta1,rhoe,theta2,phi,20000)
 answer1=answer1[begin:end,]
 answer2=mcmc2(lambda,time,Y,20000)
 answer2=answer2[begin:end,]
 gammat=quantile(answer1[,1],min(fb+(i-2)*0.05,0.5))
 gammae=quantile(answer1[,2],min(fb+(i-2)*0.05,0.5))
 rhot=quantile(answer1[,3],0.5)
 rhoe=quantile(answer1[,4],0.5)
 phi=quantile(answer1[,5],0.5)
 lambda[1]=quantile(answer2[,1],0.5)
 lambda[2]=quantile(answer2[,2],0.5)
 lambda[3]=quantile(answer2[,3],0.5)
{ if(i==2) D[i]=dose_level[2]
 else{
  D[i]=choose_dose1(gammat,gammae)
  if(is.na(D[i])){
   D[i]=choose_dose2(gammat)
  }
 }
}
 X[i]=get_doselevel(D[i])
 if(i<=12){
 S=rbind(S,potent[[loop]][[i]][(sce*6-6+X[i]),1:3])
 real_Y=rbind(real_Y,potent[[loop]][[i]][(sce*6-6+X[i]),7:9])
 time=rbind(time,potent[[loop]][[i]][(sce*6-6+X[i]),10:12])
 observe_Y=rbind(observe_Y, potent[[loop]][[i]][(sce*6-6+X[i]),4:6])
 Y=rbind(Y, generate_Y(X[i],observe_Y[i,],time[i,],lambda,gammat,gammae,rhot,theta1,rhoe,theta2,phi))
}
}
for(j in 1:7) toxi[(loop+1),j]=length(S[S==toxi[1,j]])
prob_[loop,2]=length(real_Y[real_Y==1])/length(real_Y)
DOES=choose_dose1(gammat,gammae)
DOES=get_doselevel(DOES)
Does[loop,1]=DOES
Does[loop,2:13]=X[1:12]
prob_[loop,1]=(toxi[2,6]+toxi[2,7])/36
gammat_[loop,]=c(quantile(answer1[,1],0.05),quantile(answer1[,1],0.15),quantile(answer1[,1],0.25),quantile(answer1[,1],0.5),quantile(answer1[,1],0.75),quantile(answer1[,1],0.85),quantile(answer1[,1],0.95))
gammae_[loop,]=c(quantile(answer1[,2],0.05),quantile(answer1[,2],0.15),quantile(answer1[,2],0.25),quantile(answer1[,2],0.5),quantile(answer1[,2],0.75),quantile(answer1[,2],0.85),quantile(answer1[,2],0.95))
}
R=c(R,list(list(toxi,Does,prob_,gammat_,gammae_)))
}
R

