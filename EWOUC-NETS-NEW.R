#######EWOUC-Comp######

logit<-function(p){return(log(p/(1-p)))} ####logit function####
std.dose<-seq(0.2,1,0.8/4)  ######dose level for simulation####
doselevel<-seq(0.2,1,0.8/4)

  #######################
  emuOU<-function(doselevel,gammat,gammae,N,sim,cohort,w,sname,a.T,thetae,ph,design,arrival){
    
    phi=ph ######correlation between toxicity and efficacy######  
    ##### sencario setting######
    scenario=sname     ##scenario name###
    cohort=3           ##cohort size###
    thetat<-0.333      ###target tolerated level###
    upe<-rep(0,sim)    ####Utility###
    res<-c(0,0,0)     ####initial result####    
    ######model####
    Xmin=0.2;Xmax=1;rhot=0.03;rhoe=0.08  #####initial value
    b0t=1/(gammat-Xmin)*(gammat*logit(rhot)-Xmin*logit(thetat))
    b1t=1/(gammat-Xmin)*(logit(thetat)-logit(rhot))
    b0e=1/(gammae-Xmin)*(gammae*logit(rhoe)-Xmin*logit(thetae))
    b1e=1/(gammae-Xmin)*(logit(thetae)-logit(rhoe))
    s1.pt<-1/(1+exp(-(b0t+b1t*std.dose)))  ####true toxicity based on the logistic model####
    s1.pe<-1/(1+exp(-(b0e+b1e*std.dose)))  ####true efficacy based on the logistic model #####
    d.u<-s1.pe-w*s1.pt       ####true utility###
    real.table=rbind(s1.pt,s1.pe,d.u)  ####combine true  toxicity efficacy and utility###
    
    d.table<-NULL
    sim.pe<-s1.pe
    sim.pt<-s1.pt
    #####################joint probability for toxicity and efficacy###################
    p11<-sim.pe*sim.pt*(1+(1-sim.pe)*(1-sim.pt)*(exp(phi)-1)/(1+exp(phi)))
    p10<-sim.pt-p11
    p01<-sim.pe-p11
    p00<-1-p11-p10-p01
    ####################################
    
    all.table<-NULL  ####storing all data####
    s.table<-NULL    ####storing treated patient data####
    sim.table<-NULL 
    sim.list<-NULL
    time<-rep(0,sim)  #####total time for every simulation####
    error=0            
    pstm.du=rep(0,sim)  ####posterior utility###
    
    set.seed(1234)
    r=round(runif(sim)*1000)
    final<-rep(0,sim)  ####final dose recommendation initial value####
    ##############################################
    
    ##############simulation begin##############
    for (k in 1:sim){
      fa=0.25   #####initial feasibility bound#####
      fb=0.25   
      init=list(rhot=rhot,rhoe=rhoe,gammae=gammae,gammat=gammat, phi=0,.RNG.name="base::Wichmann-Hill",.RNG.seed=r[k])
      tox.int=NULL   ####initial value for toxicity#####
      eff.int=NULL   ####initial value for efficacy#####
      
      newdata=list(eff=c(0,0,0),tox=c(0,0,0),X=c(0.2,0.2,0.2),Xmin=0.2,thetat=0.333,thetae=thetae,N=3)
      
      ############list to store MCMC posterior samples######
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
      ##########################################
      
      improve=0   ##########whether to improve toxicity bound########
      ef.im=0     ##########whether to improve efficacy bound######
      error2=0
      addT=0
      sim.seed=sample(seq(1,14,1),replace=T,size=1)
      dose.tried[1]=1  ######whether the dose has been tried#########
      for(i in 1:N) {
        
        time[k]<-time[k]+addT #####for every dose level add enrollment time######
        
        #######run jags to generate MCMC chain#########
        foo<-jags.model(textConnection(model),data=newdata,inits=init,quiet=T)
        out <- coda.samples(model=foo,variable.names=c("gammat","gammae","rhot","rhoe","phi"), n.iter=2000)
        ################################################
        
        ##########store posterior samples############
        h.gammat[[i]]<-out[[1]][,"gammat"][1001:2000]
        h.gammae[[i]]<-out[[1]][,"gammae"][1001:2000]
        h.rhot[[i]]<-out[[1]][,"rhot"][1001:2000]
        h.rhoe[[i]]<-out[[1]][,"rhoe"][1001:2000]
        h.phi[[i]]<-out[[1]][,"phi"][1001:2000]
        ##############################################
        
        #####improve feasibility bound###########
        if(fa<=0.45&ef.im==1) {fa=fa+0.05}
        if(fb<=0.45&improve==1)  {fb=fb+0.05}
        ########################################
        
        ##########choose MTD and MED#################
        mtd<-quantile(h.gammat[[i]],prob=c(fb)) ###infer mtd by its posterior quantile###
        med<-quantile(h.gammae[[i]],prob=c(fa)) ###infer med by its posterior quantile###
        if(mtd>1) mtd=1 
        if(med>1) med=1
        adj.mtd<-std.dose[max(which(mtd-std.dose>=0))]
        adj.med<-std.dose[min(which(std.dose-med>=0))]
        if(length(which(std.dose-med>0))==0) adj.med=max(std.dose)
        ##############################################
        
        sim.int[[i]]<-c(mtd,med)
        sim.ad<-NULL #####accept dose level set initialized###
        z<-NULL
        confint.l<-quantile(h.gammat[[i]]-h.gammae[[i]],prob=c(0.75)) #####posterior 75% quantile differece for mtd and med#####
        for (j in 1:length(std.dose)) {
          ##########predict posterior efficacy and toxicity############
          est.lpe[[j]]<-(1/(h.gammae[[i]]- Xmin))*(h.gammae[[i]]*log(h.rhoe[[i]]/(1-h.rhoe[[i]]))- Xmin*log(thetae/(1-thetae))+(log(thetae/(1-thetae))-log(h.rhoe[[i]]/(1-h.rhoe[[i]])))*std.dose[j])
          est.lpt[[j]]<-(1/(h.gammat[[i]]- Xmin))*(h.gammat[[i]]*log(h.rhot[[i]]/(1-h.rhot[[i]]))- Xmin*log(thetat/(1-thetat))+(log(thetat/(1-thetat))-log(h.rhot[[i]]/(1-h.rhot[[i]])))*std.dose[j])
          est.pt[[j]]<-1/(1+exp(-est.lpt[[j]]))
          est.pe[[j]]<-1/(1+exp(-est.lpe[[j]]))
          pstm.pt[j]<-mean(est.pt[[j]])
          pstm.pe[j]<-mean(est.pe[[j]])
          ###############################################################
          
          ##############choose the dose into accepted dose set#######
          if(std.dose[j]>=adj.med&std.dose[j]<=adj.mtd){        
            sim.ad<-c(std.dose[j],sim.ad)
            z<-c(mean(est.pe[[j]]-w*est.pt[[j]]),z)  
          }
          #######################################################
        }
        
        res<-c(0,0,0)
        improve<-1
        ef.im<-1
        ####################generate time to efficacy############
        addT<-rpexp(1,rate=c(3/{a.T*(3-1+0.5)},3/(a.T*(3-2+0.5)),3/{a.T*(3-3+0.5)}),t=seq(0,a.T-a.T/3,a.T/3))
        if(addT>a.T){addT=a.T}
        ########################################################
        
        ############stopping rule#############
        if (length(sim.ad)==0){
          if(confint.l[1]>=0) {
            sim.ad=adj.mtd ###if no dose in accepted dose set, treat patient at mtd####
            z=1  
          }
          if(confint.l[1]<0) break ####if p(mtd>med)>75%, stop the trial####
        }
        ######################################
        if(length(sim.ad)!=0){
          idx<-which(std.dose==sim.ad[which.max(z)])
          dose.tried[idx]=1
          nextdose<-std.dose[idx] ######choose the best utility dose level######
          for (m in 1:3){
            #####responses for next cohort#######
            res[m]<-sample(c(1,2,3,4),size=1,replace=T,prob=c(p11[idx],p10[idx],p01[idx],p00[idx]))
            if (res[m]==1) {tox=1;eff=1}
            else if (res[m]==2) {tox=1;eff=0}
            else if (res[m]==3) {tox=0;eff=1}
            else if (res[m]==4) {tox=0;eff=0}
            ###################################
            
            recuittime<-rexp(1,arrival)
            needmore<-recuittime-a.T
            if(needmore<0){needmore=0}       
            
            ########toxicity and efficacy for cumulated data#######
            newdata$tox=c(newdata$tox,tox)
            newdata$eff=c(newdata$eff,eff)
            newdata$X=c(newdata$X,nextdose)
            newdata$N=newdata$N+1
            ####################################
          }
        }
        
      }
      sim.list=list(newdata,sim=k,size=newdata$N,patientId=1:newdata$N)
      all.table<-rbind(data.frame(sim.list),all.table)###accumulate result####
      final[k]=newdata$X[newdata$N]####final dose recommendation####
      if (length(sim.ad)==0)  {
        final[k]=NA
      }
      pstm.du[k]<-mean(est.pe[[idx]]-est.pt[[idx]])
      upe[k]<-sum(est.pe[[idx]]>thetae)/length(est.pe[[idx]])
    }
    ####manipulate data to write down final result####
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
  ewouc.s1<-emuOU(doselevel,std.dose[5],std.dose[2],12,1000,3,3,"Extremely Good",3,0.3,0,"EWOUC-comp",0.25)
  ewouc.s2<-emuOU(doselevel,std.dose[4],std.dose[2],12,1000,3,3,"Good",3,0.3,0,"EWOUC-comp",7)
  ewouc.s3<-emuOU(doselevel,std.dose[5],std.dose[3],12,1000,3,2,"Moderate",3,0.3,0,"EWOUC-comp",7)
  ewouc.s4<-emuOU(doselevel,std.dose[2],std.dose[4],12,1000,3,3,"Bad",3,0.3,0,"EWOUC-comp",7)
  ewouc.s5<-emuOU(doselevel,std.dose[2],std.dose[5],12,1000,3,3,"Extremely Bad",3,0.3,0,"EWOUC-comp",7)
  
  
  
  