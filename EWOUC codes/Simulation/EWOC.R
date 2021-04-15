#######EWOC design######

####scenario 1#####
require(rjags)
require(msm) ###package to calculate time####
logit<-function(p){return(log(p/(1-p)))}

doselevel<-seq(0.2,1,0.8/4)
std.dose<-c(doselevel)

model="model{
for (i in 1:N){
# Likelihood
   tox[i]~dbern(p[i])
   logit(p[i])<- (1/(gamma - Xmin))*(gamma*logit(rho0)- Xmin*logit(theta)+(logit(theta)-logit(rho0))*X[i])
  } 
# end of for loop
# Priors
gamma ~ dunif(Xmin,1.2)
rho0 ~ dunif(0,theta)
} # end of BUGS code"


emuOU<-function(doselevel,gammat,gammae,N,sim,cohort,w,sname,a.T,thetae,ph,design,arrival){
  arrival=7
  phi=ph
  ##### sencario setting######
  scenario=sname     ##scenario name###
  cohort=3           ##cohort size###
  thetat<-0.333      ###TTL#####
 
  
  
  ############model###############
  Xmin=0.2;Xmax=1;rhot=0.03;rhoe=0.08;
  b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
  b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
  b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
  b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
  s1.pt<-1/(1+exp(-(b0t+b1t*std.dose)))  ####toxicity based on the logistic model####
  s1.pe<-1/(1+exp(-(b0e+b1e*std.dose)))  ####efficacy based on the logistic model #####
  d.u<-s1.pe-w*s1.pt       ####true utility###
  real.table=rbind(s1.pt,s1.pe,d.u)  ####combine toxicity efficacy and utility###
  d.table<-NULL
  sim.pe<-s1.pe
  sim.pt<-s1.pt
  ##########################################################
  #####################joint probability###################
  p11<-sim.pe*sim.pt*(1+(1-sim.pe)*(1-sim.pt)*(exp(phi)-1)/(1+exp(phi)))
  p10<-sim.pt-p11
  p01<-sim.pe-p11
  p00<-1-p11-p10-p01
  ##########################################################
  all.table<-NULL  ####storing all data####
  s.table<-NULL    ###stoing treated patient data####
  sim.table<-NULL   
  sim.list<-NULL
  time<-rep(0,sim)  #####total time for every simulation####
  error=0
  
  set.seed(1234)
  r=round(runif(sim)*1000)
  final<-rep(0,sim)  ####final dose recommendation####
  ################scenario setting end ##########################
  
  ##############simulation begin##############
  for (k in 1:sim){
    fa=0.25
    fb=0.25
    init=list(rho0=rhot,gamma=gammat,.RNG.name="base::Wichmann-Hill",.RNG.seed=r[k])
    tox.int=NULL
    eff.int=NULL
    newdata=list(eff=c(0,0,0),tox=c(0,0,0),X=c(0.2,0.2,0.2),Xmin=0.2,theta=0.333,thetae=thetae,N=3)
    h.gammat<-list()
    h.rhot<-list()  
     dose.tried<-rep(0,length(doselevel))
    improve=0
    ef.im=0
    error2=0
    addT=0
    sim.seed=sample(seq(1,14,1),replace=T,size=1)
    dose.tried[1]=1
    for(i in 1:N) {
      time[k]<-time[k]+addT
      foo<-jags.model(textConnection(model),data=newdata,inits=init,quiet=T)
      out <- coda.samples(model=foo,variable.names=c("gamma","rho0"), n.iter=2000)
      h.gammat[[i]]<-out[[1]][,"gamma"][1001:2000]
      h.rhot[[i]]<-out[[1]][,"rho0"][1001:2000]  
      if(fb<=0.45&improve==1)  {fb=fb+0.05}
      mtd<-quantile(h.gammat[[i]],prob=c(fb))      
      if(mtd>1) mtd=1
      adj.mtd<-std.dose[max(which(mtd-std.dose>=0))]  
      idx<-which(std.dose==adj.mtd)
      res<-c(0,0,0)
      improve<-1
      ef.im<-1
      addT<-rpexp(1,rate=c(3/{a.T*(3-1+0.5)},3/(a.T*(3-2+0.5)),3/{a.T*(3-3+0.5)}),t=seq(0,a.T-a.T/3,a.T/3))
      if(addT>a.T){addT=a.T}
        dose.tried[idx]=1
        nextdose<-std.dose[idx]
        for (m in 1:3){
          res[m]<-sample(c(1,2,3,4),size=1,replace=T,prob=c(p11[idx],p10[idx],p01[idx],p00[idx]))
          if (res[m]==1) {tox=1;eff=1}
          else if (res[m]==2) {tox=1;eff=0}
          else if (res[m]==3) {tox=0;eff=1}
          else if (res[m]==4) {tox=0;eff=0}
          recuittime<-rexp(1,arrival)
          needmore<-recuittime-a.T
          if(needmore<0){needmore=0}              
          newdata$tox=c(newdata$tox,tox)
          newdata$eff=c(newdata$eff,eff)
          newdata$X=c(newdata$X,adj.mtd)
          newdata$N=newdata$N+1
        }
    }    
    
    sim.list=list(newdata,sim=k,size=newdata$N,patientId=1:newdata$N)
    all.table<-rbind(data.frame(sim.list),all.table)
    final[k]=newdata$X[newdata$N]
  }
  
  s.table<-all.table[all.table$patientId != c(37,38,39),]

  dr.table<-table(final,useNA="always")/sim
  
  df.table<-table(s.table[,3])/nrow(s.table)
  
  t.table<-table(s.table[,2])/nrow(s.table)
  
  e.table<-table(s.table[,1])/nrow(s.table)
  
  
  n.table<-all.table[all.table$size <=36,]
  
  aver.samplesize<-nrow(s.table)/sim

  u.table<-matrix(c(d.u,std.dose),ncol=2)
  colnames(u.table)<-c("utility","d")
  ut<-merge(s.table,u.table,by.x="X",by.y="d")
  eut<-sum(ut$utility)/nrow(ut)
  r.rt<-cbind(t(round(real.table,2)),std.dose)
  cm<-data.frame(r.rt)
  
  c1<-data.frame(df.table);colnames(c1)<-c("dose","Tp")
  c2<-data.frame(dr.table);colnames(c2)<-c("dose","Rp")
  c3<-data.frame(t.table);colnames(c3)<-c("t","Dp")
  c4<-data.frame(e.table);colnames(c4)<-c("e","Ep")
  expt<-mean(time)
  
  m1<-merge(cm,c1,by.x="std.dose",by.y="dose",all.x=T)
  m2<-merge(m1,c2,by.x="std.dose",by.y="dose",all.x=T)
  m3<-cbind(m2,c3[c3$t==1,2],c4[c4$e==1,2],eut,expt,sname,design, aver.samplesize)
  m3[,5:6]=m3[,5:6]*100
  
  colnames(m3)[7:8]=c("DLT","Efficacy")
  
  write.table(m3,"EWOUC2.csv",na="0",dec=".",qmethod="double",sep=",",col.names=F,row.names=T,quote=T,append="T")
  
  
  return(list(m3))
}



###
ewoc.s1<-emuOU(doselevel,std.dose[5],std.dose[2],12,1000,3,3,"Extremely Good",3,0.3,0,"EWOUC-comp",0.25)
ewoc.s2<-emuOU(doselevel,std.dose[4],std.dose[2],12,1000,3,3,"Good",3,0.3,0,"EWOUC-comp",7)
ewoc.s3<-emuOU(doselevel,std.dose[5],std.dose[3],12,1000,3,2,"Moderate",3,0.3,0,"EWOUC-comp",7)
ewoc.s4<-emuOU(doselevel,std.dose[2],std.dose[4],12,1000,3,3,"Bad",3,0.3,0,"EWOUC-comp",7)
ewoc.s5<-emuOU(doselevel,std.dose[2],std.dose[5],12,1000,3,3,"Extremely Bad",3,0.3,0,"EWOUC-comp",7)







