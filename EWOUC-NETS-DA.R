######EWOUC-NETS late on efficacy######
require(rjags) ####package for posterior distribution########
require(msm)
library(openxlsx)
####jags model specification####
model="model
{
for (i in 1:N)
{
# Logistic regression model for extensification
eff[i]~ dnorm(mu[i] + tau*(S[i]-xi[i]), 1/(sigma[i]))T(0,1)
S[i] ~ dnorm(xi[i], 1/sqrt((xi[i]*(1-xi[i])/N)))T(0,1)
logit(xi[i])<-(1/(gammat - Xmin))*(gammat*logit(rhot)- Xmin*logit(thetat)+(logit(thetat)-logit(rhot))*X[i])
logit(mu[i])<-(1/(gammae - Xmin))*(gammae*logit(rhoe)- Xmin*logit(thetae)+(logit(thetae)-logit(rhoe))*X[i])
sigma[i] = sqrt(X[i])
M[i]~dbern(pm[i])
pm[i]<-(1-eff[i]+eff[i]*exp(-lambda[1]*s[i,1]-lambda[2]*s[i,2]-lambda[3]*s[i,3]))*(1-C[i])
}

lambda[1]~dgamma(0.2,0.5)
lambda[2]~dgamma(0.33,0.5)
lambda[3]~dgamma(1,0.5)
tau~dunif(-1,1)
rhot~dunif(0,0.333)
rhoe~dunif(0,0.5)
gammat~dunif(Xmin,1.2)
gammae~dunif(Xmin,1.2)
}"
  ###############
  
  std.dose<-seq(0.2,1,0.2) ####dose levels######
  t.idx<-c(5,4,5,2,2) ####the index of MTD#####
  e.idx<-c(2,2,3,4,5)####the index of MED#######
  logit<-function(p){return(log(p/(1-p)))}
  
  
  simLate<-function(gammat,gammae,w,design,sname){
    N=12 ######sample size, 12 cohort each with 3 patients######
    
    #####simulation setting#######
    tau=0
    lambda=c(0.4,0.667,2)
    thetat<-0.476 ###
    thetae<-0.3
    Xmin=0.2;Xmax=1;rhot=0.01;rhoe=0.02  #####initial value
    b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
    b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
    b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
    b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
    xi = exp(b0t + b1t*std.dose)/(1+exp(b0t + b1t*std.dose))
    mu = exp(b0e + b1e*std.dose)/(1+exp(b0e + b1e*std.dose))
    lambda = 7
    sig2 = std.dose^lambda
    s1.pe = mu
    s1.pt = xi
    d.u<-s1.pe-w*s1.pt       ####true utility###
    real.table=rbind(s1.pt,s1.pe,d.u) ####scenario###
    
    d.table<-NULL
    all.table<-NULL
    
    #################################
    s.table<-NULL
    sim=1000
    set.seed(1234)
    r=round(runif(sim)*1000)###seed###
    a.T<-3########observation window#####
    K<-3#######3 knots######
    
    c.int<-seq(0,a.T,a.T/K)
    res=rep(0,3)
    final=rep(0,sim) ####final dose recommendation####
    
    ###########simulation begin#############
    for (k in 1:sim){
      init=list(rhot=rhot,rhoe=rhoe,gammae=gammat,gammat=gammae, tau=0, lambda=c(0.2,0.4,2),.RNG.name="base::Wichmann-Hill",.RNG.seed=r[k]) #####initial###
      S.int=NULL
      eff.int=NULL
      ####first cohort response of the data####
      newdata=list(eff=c(0,0,0),M=c(0,0,0),S=c(0,0,0),X=c(0.2,0.2,0.2),Xmin=0.2,thetat=0.476,thetae=0.3,N=3,C=c(1,1,1),s=matrix(c(1,1,1,1,1,1,1,1,1),ncol=3,nrow=3))
      #####################  
      
      #####initialize lists to restore posterior samples#######
      h.gammat<-list()
      h.gammae<-list()
      h.rhot<-list()
      h.rhoe<-list()
      h.tau<-list()
      est.lpe<-list()
      est.lpt<-list()
      est.pe = list()
      est.pt = list()
      est.mu<-list()
      est.xi<-list()
      sim.int<-list()
      pstm.pt<-rep(0,length(std.dose))
      pstm.pe<-rep(0,length(std.dose))
      ##################################
      
      dose.tried<-rep(0,length(std.dose))#######no skipping, indicate whether the dose level has been tried or not######
      fa=0.25 ######feasibility bound for toxicity####
      fb=0.25######feasibility bound for efficacy#####
      improve=0####whether to increase the feasibility bound for toxicity###
      ef.im=0###whether to increase the feasibility bound for efficacy####
      phi=0.5
      u<-rep(3,3) ####observation window for three cohorts####
      newdata1<-newdata
      for(i in 1:N) {
        ########## Update posterior sample by MCMC#######
        foo<-jags.model(textConnection(model),data=newdata1,inits=init,quiet=T)
        out <- coda.samples(model=foo,variable.names=c("gammat","gammae","rhot","rhoe","tau","lambda"), n.iter=2000)
        h.gammat[[i]]<-out[[1]][,"gammat"][1001:2000]
        h.gammae[[i]]<-out[[1]][,"gammae"][1001:2000]
        h.rhot[[i]]<-out[[1]][,"rhot"][1001:2000]
        h.rhoe[[i]]<-out[[1]][,"rhoe"][1001:2000]
        h.tau[[i]]<-out[[1]][,"tau"][1001:2000]
        if(fa<=0.45&ef.im==1) {fa=fa+0.05}
        if(fb<=0.45&improve==1)  {fb=fb+0.05}
        mtd<-quantile(h.gammat[[i]],prob=c(fb)) ###inference for MTD###
        med<-quantile(h.gammae[[i]],prob=c(fa)) ####inference for MED####
        if(mtd>1) {mtd=1}
        if(med>1) {med=1}
        adj.mtd<-std.dose[max(which(mtd-std.dose>=0))]###choose the one in pre-specified dose level most cloest to MTD#####
        adj.med<-std.dose[min(which(std.dose-med>=0))]
        sim.int[[i]]<-c(mtd,med) ### the interval for MTD and MED###
        sim.ad<-NULL ####accepted sets initialized####
        z<-NULL
        confint.l<-quantile(h.gammat[[i]]-h.gammae[[i]],prob=c(0.75))
        for (j in 1:length(std.dose)) {
          est.lpe[[j]]<-(1/(h.gammae[[i]]- Xmin))*(h.gammae[[i]]*log(h.rhoe[[i]]/(1-h.rhoe[[i]]))- Xmin*log(thetae/(1-thetae))+(log(thetae/(1-thetae))-log(h.rhoe[[i]]/(1-h.rhoe[[i]])))*std.dose[j])
          est.lpt[[j]]<-(1/(h.gammat[[i]]- Xmin))*(h.gammat[[i]]*log(h.rhot[[i]]/(1-h.rhot[[i]]))- Xmin*log(thetat/(1-thetat))+(log(thetat/(1-thetat))-log(h.rhot[[i]]/(1-h.rhot[[i]])))*std.dose[j])
          est.xi[[j]]<-mean(1/(1+exp(-est.lpt[[j]])))
          est.mu[[j]]<-mean(1/(1+exp(-est.lpe[[j]])))
          est.pt[[j]]<-est.xi[[j]]
          est.pe[[j]]<-est.mu[[j]]
          pstm.pt[j]<-mean(est.pt[[j]])
          pstm.pe[j]<-mean(est.pe[[j]])
          if(std.dose[j]>=adj.med&std.dose[j]<=adj.mtd){
            sim.ad<-c(std.dose[j],sim.ad)
            z<-c(mean(est.pe[[j]]-w*est.pt[[j]]),z) ###predictive utility####
          }
        }
        u[newdata$M==1]= u[newdata$M==1]+phi ###follow up time, for missing data increase 3 months###
        for (p.i in 1:length(newdata$M)){
          if(u[p.i]<a.T){
            for(c.k in 1:3){
              if (u[p.i]>=c.int[c.k] & u[p.i]<c.int[c.k+1]){newdata$s[p.i,c.k]=u[p.i]-c.int[c.k]}
              if (u[p.i]<c.int[c.k]) {newdata$s[p.i,c.k]=0}
              if (u[p.i]>=c.int[c.k+1]){newdata$s[p.i,c.k]=1}
            }
          }
          else {newdata$M[p.i]=0;newdata$C[p.i]=1}   
        }
        m.p<-exp(-lambda[1]*newdata$s[,1]-lambda[2]*newdata$s[,2]-lambda[3]*newdata$s[,3])###the probability for missing####
        missing.idx<-which(newdata$M==1&newdata$eff==1)
        if(length(missing.idx)>0){
          for(h in 1:length(missing.idx)){
            s.mis<-missing.idx[h]
            newdata$M[s.mis]<-rbinom(1,1,m.p[s.mis])####whether the data is still missing####
          }      
        }    
        res<-c(0,0,0)
        improve<-1
        ef.im<-1
        ###############stopping rule##############
        if (length(sim.ad)==0){
          if(confint.l[1]>=0) {
            sim.ad=adj.mtd ### if no recommendation when data is accumulating, then mtd#####
            z=1  
          }
          if(confint.l[1]<0) break ######if p(mtd < med)>0.75, then stop #####
        }
        #################################
        if(length(sim.ad)!=0){
          idx<-which(std.dose==sim.ad[which.max(z)])
          dose.tried[idx]=1
          nextdose<-std.dose[idx]
          
          #####generate new data based on true distribution####
          for (m in 1:3){
            set.seed(m+k+10*i+55*j)
            S = rtnorm(n = 1, mean = xi[idx], sd = sqrt(xi[idx]*(1-xi[idx])/N), lower = 0, upper = 1)
            eff = rtnorm(n = 1, mean = mu[idx] + tau*(S-xi[idx]), sd = sqrt(sqrt(std.dose[idx]^lambda)), lower = 0, upper = 1) 
            M = 1
            newdata$S=c(newdata$S,S)
            newdata$eff=c(newdata$eff,eff)
            newdata$X=c(newdata$X,nextdose)
            newdata$N=newdata$N+1 
            newdata$M=c(newdata$M,M)
            newdata$s=rbind(newdata$s,c(0,0,0))
            newdata$C=c(newdata$C,0)
            u=c(u,0)      
          }    
          newdata1<-newdata
          newdata1$eff[newdata1$M==1]<-NA ###make those missing data with response NA####   
        }
        else {break}
      }
      
      sim.list=list(newdata,sim=k,size=newdata$N,patientId=1:newdata$N)
      all.table<-rbind(data.frame(sim.list),all.table)
      final[k]=newdata$X[newdata$N] ###final dose recommendation####
      if (length(sim.ad)==0)  {
        final[k]=NA
      }
    }
    
    ####summarize result to write down in csv file###
    s.table<-all.table[all.table$patientId != c(37,38,39),]
    dr.table<-table(final,useNA="always")/sim
    df.table<-table(s.table[,4])/nrow(s.table)
    t.table<-nrow(subset(s.table, S >= thetat))/nrow(s.table)
    e.table<-nrow(subset(s.table, eff >= thetae))/nrow(s.table)
    p.table<-rep(0,length(std.dose))
    p.table[which(std.dose==as.numeric(names(df.table)[which.max(df.table)]))]=19
    u.table<-matrix(c(d.u,std.dose),ncol=2)
    colnames(u.table)<-c("utility","d")
    ut<-merge(s.table,u.table,by.x="X",by.y="d")
    eut<-sum(ut$utility)/nrow(ut)
    avgsample<-nrow(s.table)/sim
    r.rt<-cbind(t(round(real.table,2)),std.dose)
    cm<-data.frame(r.rt)
    c1<-data.frame(df.table);colnames(c1)<-c("dose","Tp")
    c2<-data.frame(dr.table);colnames(c2)<-c("dose","Rp")
    c3<-data.frame(t.table);colnames(c3)<-c("DLT")
    c4<-data.frame(e.table);colnames(c4)<-c("Efficacy")
    m1<-merge(cm,c1,by.x="std.dose",by.y="dose",all.x=T)
    m2<-merge(m1,c2,by.x="std.dose",by.y="dose",all.x=T)
    m3<-cbind(m2,c3,c4,eut,sname,design,avgsample)
    m3[,5:6]=m3[,5:6]*100
    colnames(m3)[7:8]=c("DLT","Efficacy")
    write.xlsx(m3,"EWOUC-NETS-DA-S3-new.xlsx")
    ##################
    return(m3)
  }
  
  
  ####simulate five scenarios####
  lat.s1<-simLate(std.dose[t.idx[1]],std.dose[e.idx[1]],2,"EWOUC-DA","Extremely Good")
  lat.s2<-simLate(std.dose[t.idx[2]],std.dose[e.idx[2]],2,"EWOUC-DA","Good")
  lat.s3<-simLate(std.dose[t.idx[3]],std.dose[e.idx[3]],2,"EWOUC-DA","Moderate")
  lat.s4<-simLate(std.dose[t.idx[4]],std.dose[e.idx[4]],2,"EWOUC-DA","Bad")
  lat.s5<-simLate(std.dose[t.idx[5]],std.dose[e.idx[5]],2,"EWOUC-DA","Extremely Bad")