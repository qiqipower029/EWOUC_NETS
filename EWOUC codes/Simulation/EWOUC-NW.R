
########EWOUC Not Wait#######
require(rjags)

doselevel<-seq(0.2,1,0.2)
std.dose<-doselevel
t.idx<-c(5,4,5,2,2)
e.idx<-c(2,2,3,4,5)
logit<-function(p){return(log(p/(1-p)))} ####logit function###

#######jags model############
model="model
{
for (i in 1:N)
{
# Logistic regression model for extensification
tox[i] ~ dbern(P1[i])
eff[i]~ dbern(P2[i])
logit(P1[i])<-(1/(gammat - Xmin))*(gammat*logit(rhot)- Xmin*logit(thetat)+(logit(thetat)-logit(rhot))*X[i])
logit(P2[i])<-(1/(gammae - Xmin))*(gammae*logit(rhoe)- Xmin*logit(thetae)+(logit(thetae)-logit(rhoe))*X[i])
P11[i] <- P1[i]*P2[i]*(1+(1-P1[i])*(1-P2[i])*(exp(phi)-1)/(1+exp(phi)))
P10[i] <- P1[i]-P11[i]
P01[i] <- P2[i]-P11[i]
P00[i] <- 1-(P10[i]+P01[i]+P11[i])
}

phi~dnorm(0,1)
rhot~dunif(0,0.333)
rhoe~dunif(0,0.5)
gammat~dunif(Xmin,1.2)
gammae~dunif(Xmin,1.2)
}"
############################

#####simulation function########
simLate<-function(gammat,gammae,w,design,sname){
  N=12 ####12 cohorts, 3 patients with each cohort####
 
  phi=0
  lambda=c(0.4,0.667,2) ###true lambda###
  thetat<-0.333 #####target tolerated level###
  thetae<-0.3  ######efficacy lower bound####
 
  #####model setting###
  Xmin=0.2;Xmax=1;rhot=0.03;rhoe=0.08;
  b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
  b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
  b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
  b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
  s1.pt<-1/(1+exp(-(b0t+b1t*std.dose)))
  s1.pe<-1/(1+exp(-(b0e+b1e*std.dose)))
  d.u<-s1.pe-w*s1.pt
###############################################

  real.table=rbind(s1.pt,s1.pe,d.u) #########combind table for true toxicity, efficacy and utility########
  d.table<-NULL
  sim.pe<-s1.pe
  sim.pt<-s1.pt
  
  all.table<-NULL
  p11<-sim.pe*sim.pt*(1+(1-sim.pe)*(1-sim.pt)*(exp(phi)-1)/(1+exp(phi)))
  p10<-sim.pt-p11
  p01<-sim.pe-p11
  p00<-1-p11-p10-p01
  s.table<-NULL
  
  sim=1000 ########1000 replicates#######
  set.seed(1234)
  r=round(runif(sim)*1000)#######seed#######
  a.T<-3 ####assessment period###
  K<-3 
  c.int<-seq(0,a.T,a.T/K)
  res=rep(0,3)
  final=rep(0,sim)
  for (k in 1:sim){
    init=list(rhot=rhot,rhoe=rhoe,gammae=gammat,gammat=gammae, phi=0.1, lambda=c(0.2,0.4,2),.RNG.name="base::Wichmann-Hill",.RNG.seed=r[k])
    tox.int=NULL
    eff.int=NULL
   
    #####initial data###
    newdata=list(eff=c(0,0,0),M=c(0,0,0),tox=c(0,0,0),X=c(0.2,0.2,0.2),Xmin=0.2,thetat=0.333,thetae=0.3,N=3,C=c(1,1,1),s=matrix(c(1,1,1,1,1,1,1,1,1),ncol=3,nrow=3))
    h.gammat<-list()
    h.gammae<-list()
    h.rhot<-list()
    h.rhoe<-list()
    h.phi<-list()
    est.lpe<-list()
    est.lpt<-list()
    est.pe<-list() 
    est.pt<-list()
    sim.int<-list()
    pstm.pt<-rep(0,length(doselevel))
    pstm.pe<-rep(0,length(doselevel))
    dose.tried<-rep(0,length(doselevel))
    fa=0.25 ####feasibility bound###
    fb=0.25 
    improve=0 #####whether to improve feasibility bounds####
    ef.im=0
    tau=0.5 #####patient inter arrival####
    u<-rep(3,3) 
    newdata1<-newdata
    for(i in 1:N) {
      foo<-jags.model(textConnection(model),data=newdata1,inits=init,quiet=T)
      out <- coda.samples(model=foo,variable.names=c("gammat","gammae","rhot","rhoe","phi","lambda"), n.iter=2000)
      ######posterior samples######
      h.gammat[[i]]<-out[[1]][,"gammat"][1001:2000]
      h.gammae[[i]]<-out[[1]][,"gammae"][1001:2000]
      h.rhot[[i]]<-out[[1]][,"rhot"][1001:2000]
      h.rhoe[[i]]<-out[[1]][,"rhoe"][1001:2000]
      h.phi[[i]]<-out[[1]][,"phi"][1001:2000]
     #########increase feasibility bound until 0.5####
      if(fa<=0.45&ef.im==1) {fa=fa+0.05}
      if(fb<=0.45&improve==1)  {fb=fb+0.05}
     ################################################
      mtd<-quantile(h.gammat[[i]],prob=c(fb)) ####infer mtd by posterior quantile###
      med<-quantile(h.gammae[[i]],prob=c(fa)) ####infer med by posterior quantile###
      if(mtd>1) {mtd=1}
      if(med>1) {med=1}
     
      adj.mtd<-std.dose[max(which(mtd-std.dose>=0))]####choose the dose level most close to the prespecified one###########
      adj.med<-std.dose[min(which(std.dose-med>=0))]
      sim.int[[i]]<-c(mtd,med)
      sim.ad<-NULL
      z<-NULL
      confint.l<-quantile(h.gammat[[i]]-h.gammae[[i]],prob=c(0.75))
     ####calculate posterior predictive utility#####
      for (j in 1:length(std.dose)) {
        est.lpe[[j]]<-(1/(h.gammae[[i]]- Xmin))*(h.gammae[[i]]*log(h.rhoe[[i]]/(1-h.rhoe[[i]]))- Xmin*log(thetae/(1-thetae))+(log(thetae/(1-thetae))-log(h.rhoe[[i]]/(1-h.rhoe[[i]])))*std.dose[j])
        est.lpt[[j]]<-(1/(h.gammat[[i]]- Xmin))*(h.gammat[[i]]*log(h.rhot[[i]]/(1-h.rhot[[i]]))- Xmin*log(thetat/(1-thetat))+(log(thetat/(1-thetat))-log(h.rhot[[i]]/(1-h.rhot[[i]])))*std.dose[j])
        est.pt[[j]]<-1/(1+exp(-est.lpt[[j]]))
        est.pe[[j]]<-1/(1+exp(-est.lpe[[j]]))
        pstm.pt[j]<-mean(est.pt[[j]])
        pstm.pe[j]<-mean(est.pe[[j]])
        if(std.dose[j]>=adj.med&std.dose[j]<=adj.mtd){
          sim.ad<-c(std.dose[j],sim.ad)
          z<-c(mean(est.pe[[j]]-w*est.pt[[j]]),z) 
        }
      }
     ##############################################
      u[newdata$M==1]= u[newdata$M==1]+tau####################follow up time increase for those with missing reponses#####
     ####calculate time for piece wise exponential distribution for jags model#####
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
     ###############################################################################
     
     ####calculate missing probability for response#######
      m.p<-exp(-lambda[1]*newdata$s[,1]-lambda[2]*newdata$s[,2]-lambda[3]*newdata$s[,3])
      missing.idx<-which(newdata$M==1&newdata$eff==1)
      if(length(missing.idx)>0){
        for(h in 1:length(missing.idx)){
          s.mis<-missing.idx[h]
          newdata$M[s.mis]<-rbinom(1,1,m.p[s.mis]) ###whether data is still missing####
        }
      }
      #########################################################
      res<-c(0,0,0)
      improve<-1
      ef.im<-1
      ###############stopping rule##############
      if (length(sim.ad)==0){
        if(confint.l[1]>=0) {
          sim.ad=adj.mtd
          z=1  
        }
        if(confint.l[1]<0) break
      }
      #################################
      if(length(sim.ad)!=0){
        idx<-which(std.dose==sim.ad[which.max(z)])
        dose.tried[idx]=1
        nextdose<-std.dose[idx]
        for (m in 1:3){
          res[m]<-sample(c(1,2,3,4),size=1,replace=T,prob=c(p11[idx],p10[idx],p01[idx],p00[idx]))
          if (res[m]==1) {tox=1;eff=1;M=1;}
          else if (res[m]==2) {tox=1;eff=0;M=1;}
          else if (res[m]==3) {tox=0;eff=1;M=1;}
          else if (res[m]==4) {tox=0;eff=0;M=1;}
         
          newdata$tox=c(newdata$tox,tox)
          newdata$eff=c(newdata$eff,eff)
          newdata$X=c(newdata$X,nextdose)
          newdata$N=newdata$N+1
          newdata$M=c(newdata$M,M)
          newdata$s=rbind(newdata$s,c(0,0,0))
          newdata$C=c(newdata$C,0)
          u=c(u,0)      
        }
        newdata1<-newdata #######update data by true data####
        newdata1$eff[newdata1$M==1]<-0####treat missing outcome as non effect####
      }
      else {break}
    }
    
    sim.list=list(newdata,sim=k,size=newdata$N,patientId=1:newdata$N) ####accumulate data###
    all.table<-rbind(data.frame(sim.list),all.table)
    final[k]=newdata$X[newdata$N]####final dose recommendation###
    if (length(sim.ad)==0)  {
      final[k]=NA ####if early terminaated, the result is NA#####
    }
  }
  s.table<-all.table[all.table$patientId != c(37,38,39),]
  dr.table<-table(final,useNA="always")/sim
  df.table<-table(s.table[,4])/nrow(s.table)
  t.table<-table(s.table[,3])/nrow(s.table)
  e.table<-table(s.table[,1])/nrow(s.table)
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
  c3<-data.frame(t.table);colnames(c3)<-c("t","Dp")
  c4<-data.frame(e.table);colnames(c4)<-c("e","Ep")
  
  m1<-merge(cm,c1,by.x="std.dose",by.y="dose",all.x=T)
  m2<-merge(m1,c2,by.x="std.dose",by.y="dose",all.x=T)
  m3<-cbind(m2,c3[c3$t==1,2],c4[c4$e==1,2],eut,sname,design,avgsample)
  m3[,5:6]=m3[,5:6]*100
  colnames(m3)[7:8]=c("DLT","Efficacy")
  write.table(m3,"EWOUC2.csv",na="0",dec=".",qmethod="double",sep=",",col.names=F,row.names=T,quote=T,append="T")
  return(m3)
}

lat.s1<-simLate(std.dose[t.idx[1]],std.dose[e.idx[1]],3,"EWOUC-NW2","Extremely Good")
lat.s2<-simLate(std.dose[t.idx[2]],std.dose[e.idx[2]],3,"EWOUC-NW2","Good")
lat.s3<-simLate(std.dose[t.idx[3]],std.dose[e.idx[3]],2,"EWOUC-NW2","Moderate")
lat.s4<-simLate(std.dose[t.idx[4]],std.dose[e.idx[4]],3,"EWOUC-NW2","Bad")
lat.s5<-simLate(std.dose[t.idx[5]],std.dose[e.idx[5]],3,"EWOUC-NW2","Extremely Bad")

