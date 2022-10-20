library(survival)
library(partlyconditional)


fp.func3=function(sk,t,theta){
  AUC.st.mat=c()
  for (s.curr in sk)
  { 
    tsvec=cbind(1, sqrt(s.curr), sqrt(t), sqrt(s.curr)^2, sqrt(t)^2, sqrt(s.curr*t), sqrt(t)^3, sqrt(s.curr)^3, (sqrt(t)^2)*sqrt(s.curr), sqrt(t)*(sqrt(s.curr)^2))
    teta=as.vector(as.matrix(tsvec)%*%theta)
    teta[teta>=100]=100
    teta[teta<=-100]=-100
    AUC.st=exp(teta)/(1+exp(teta))
    AUC.st.mat=rbind(AUC.st.mat, AUC.st)
  }
  return(AUC.st.mat)}  

logL=function(theta,sk,d.times,n1.h.mat,n2.h.mat,eta.sk.mat,sk.le.dh,Ewt.die){
  AUC.seq.mat=fp.func3(sk,d.times,theta) 
  logL.seq= Ewt.die*eta.sk.mat*sk.le.dh*(n1.h.mat*log(AUC.seq.mat)+n2.h.mat*log(1-AUC.seq.mat))
  logL=-sum(logL.seq)}

#PC.Cox.wt is the PC.Cox function with weight
PC.Cox.wt=function (id, stime, status, measurement.time, predictors, data, 
                    additional.formula.pars = NULL, da.short, input_weight) 
{
  call <- match.call()
  stopifnot(is.character(stime))
  stopifnot(is.character(status))
  stopifnot(is.data.frame(data))
  tmpnames <- c(id, stime, status, measurement.time, predictors)
  if (!all(is.element(tmpnames, names(data)))) 
    stop(paste("'", tmpnames[which(!is.element(tmpnames, 
                                               names(data)))], "' cannot be found in data.frame provided", 
               sep = ""))
  mycomplete <- complete.cases(data[, tmpnames])
  if (nrow(data) != sum(mycomplete)) {
    warning(paste(nrow(data) - sum(mycomplete), "observation(s) were removed due to missing data \n New sample size is now:", 
                  sum(mycomplete)))
    data <- data[mycomplete, ]
  }
  wonky.times <- data[[tmpnames[2]]] < data[[tmpnames[4]]]
  if (any(wonky.times)) {
    cat(paste0("... removing ", sum(wonky.times), " observations where outcome time", 
               tmpnames[2], " is less than measurement time ", 
               tmpnames[4], ".\n"))
  }
  data <- data[!wonky.times, ]
  t.star <- data[[stime]] - data[[measurement.time]]
  my.formula <- as.formula(paste0("Surv( t.star, ", status, 
                                  ") ~", paste(predictors, collapse = " + "), 
                                  "+", "cluster(", id, ")", ifelse(is.null(additional.formula.pars), 
                                                                   "", paste("+ ", additional.formula.pars))))
  
  da.short$wt=input_weight
  da.short.wt=da.short[, c("obs_id", "wt")]
  data.wt <- merge(data, da.short.wt,by="obs_id")
  
  fit <- coxph(my.formula, data = data.wt, x = TRUE, weights = wt)
  predictors.unique <- predictors[!is.element(predictors, c(id, stime, status, measurement.time))]
  out <- list(model.fit = fit, call = call, variable.names = c(id, stime, status, measurement.time, predictors.unique))
  class(out) <- "PC_cox"
  out
}


#Backtofull is used to re-organize the datasets
# the input datasets require da.short and da.long
# da.short should include variable obs_id, Y(observed survival time), delta, followed by baseline variables 
# da.long should include variables  obs_id, vtime, Y, delta followed by longitudinal variables. 
# the function will calcu;ate the number of baseline variables, and rename the baseline variables as X1, X2,.....
# the function will calculate the number of the longitudinal variables, and rename the longitudinal variables as Y1, Y2, ....
Backtofull=function(da.short, da.long, sk)
{ 
  s.num=length(sk)
  n=dim(da.short)[1]
  n.Yt=dim(da.long)[2]-4
  n.X=dim(da.short)[2]-3
  Yt.names=paste0("Y", seq(1:n.Yt))
  X.names =paste0("X", seq(1:n.X))
  colnames(da.short)=c(c("obs_id", "Y", "delta"), X.names)
  colnames(da.long) =c(c("obs_id", "vtime", "Y", "delta"), Yt.names)
  da.long.full=cbind(rep(sk, n),rep(seq(1,n), each=s.num), rep(1, n*s.num) )
  colnames(da.long.full)=c("vtime", "obs_id", "eta")
  da.long1=merge(da.long, da.long.full ,by=c("obs_id", "vtime" ), all = T) 
  da.long1$eta[is.na(da.long1[,c(Yt.names[1])])]=0
  da.long2=merge(da.long1[,c(c("obs_id", "vtime", "eta"),Yt.names)], da.short[,c(c("obs_id", "Y" , "delta"), X.names)], by=c("obs_id"))
  da.long2=da.long2[, c(c("obs_id", "Y", "delta", "vtime",  "eta"), X.names, Yt.names )]
  da.short1=da.long2[!duplicated(da.long2$obs_id), c("obs_id", "delta", "Y")]
  return(list(da.short=da.short,da.long=da.long2, X.names=X.names, Yt.names=Yt.names, n.X=n.X, n.Yt=n.Yt))}

main1.sub.func=function(da.short, da.long, sk, tseq.eval, par0,  N.CV, resample, nsap, d.per, seed=12345){
  s.num=length(sk)
  t.num=length(tseq.eval)
  n=dim(da.short)[1]
  
  data0=Backtofull(da.short=da.short, da.long=da.long, sk=sk)
  da.CV=data0$da.short
  marker.CV=data0$da.long
  X.names=data0$X.names
  Yt.names=data0$Yt.names
  n.X=data0$n.X
  n.Yt=data0$n.Yt
  
  if (is.null(d.per)){ d.per=quantile(da.CV$Y[da.CV$delta==1],0.9)} 
  if (is.null(d.per)==F){d.per=d.per}
  
  AUC.rep=c()   #save the repeated CV AUC(s,t) values
  Vseq.list=list() #save the binary variable which splits the dataset into training and testing sets.
  par.rep=c()       #save all repeated CV thetas (3 repeat * 2 fold-> 6 rows)
  for (rep.CV in 1:N.CV){
    Vseq=rbinom(n,1,0.4)+1  # +1 to make the group indicator become 1&2 instead of 0&1
    AUC.fd=c()      #save the AUC(s,t) for the two fold CV.
    par.fd=c()      #save all CV thetas (2 fold-> 2 rows)
    for (fd in 1:2) {
      da.Tr=da.CV[Vseq==fd,] 
      marker.Tr=merge(da.Tr[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime", "eta"), X.names, Yt.names)], by=c("obs_id"))
      marker.Tr$logvt=log(marker.Tr$vtime+1)
      pc.coef.Tr <-  PC.Cox(
        id = "obs_id",
        stime = "Y",
        status = "delta",
        measurement.time = "vtime",  ##survival time and measurement time must be on the same scale.
        predictors =c(c("logvt"),X.names, Yt.names),
        data = marker.Tr[marker.Tr$eta==1,])$model.fit$coefficients[-1]  
      da.Te=da.CV[Vseq==3-fd,]
      marker.Te=merge(da.Te[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime", "eta"), X.names, Yt.names)], by=c("obs_id"))
      pc.coef.Tr.mat=t(matrix(rep(pc.coef.Tr, dim(marker.Te)[1] ), ncol =  dim(marker.Te)[1], nrow=n.X+n.Yt))
      marker.Te$Zt=apply(pc.coef.Tr.mat*marker.Te[,c(X.names,Yt.names)], 1, sum)
      
      ########################################Compute EST of AUC(s,t)######################################
      #We use the da.Te as the survival data
      #We use marker.Te as the biomarker data
      
      n.Te=dim(da.Te)[1]
      d.times=sort(unique(da.Te$Y[da.Te$delta==1 & da.Te$Y<=d.per]))
      md=length(d.times)
      Ymat=matrix(rep(da.Te$Y,md),n.Te)
      dt.mat=matrix(rep(d.times,rep(n.Te,md)),n.Te)
      dt.mat.tb=dt.mat+0.00001; 
      n.h=apply((Ymat>dt.mat.tb),2,sum)
      death.mat=(Ymat==dt.mat)*da.Te$delta #currently do not allow tied failure times;
      marker.Te$eta[marker.Te$eta==0]=0
      RS.data=marker.Te
      
      RS.seq.list=RS.mat.list=RS.d.list=RS.dmat.list=n1.h.list=n2.h.list=list()
      for (i in 1:s.num)
      {
        RS.seq.list[[i]]=RS.data[RS.data$vtime==sk[i],]
        RS.mat.list[[i]]<-matrix(rep(RS.seq.list[[i]]$Zt,md),n.Te)
        RS.d.list[[i]]=apply(RS.mat.list[[i]]*death.mat,2,sum,na.rm=TRUE)
        RS.d.list[[i]][RS.d.list[[i]]==0]<-NA
        RS.dmat.list[[i]]=matrix(rep(RS.d.list[[i]],rep(n.Te,md)),n.Te)
        n1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)
        n2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb),2,sum,na.rm=TRUE)
      }
      n1.h.mat=n2.h.mat=c()
      eta.sk.mat=c()
      for (i in 1:s.num)
      {
        n1.h.mat=rbind(n1.h.mat,n1.h.list[[i]])
        n2.h.mat=rbind(n2.h.mat,n2.h.list[[i]])
        eta.sk.mat=rbind(eta.sk.mat, RS.d.list[[i]])
      }  
      eta.sk.mat=(eta.sk.mat/eta.sk.mat)
      eta.sk.mat[is.na(eta.sk.mat)]=0
      sk.mat=matrix(rep(sk,each=md), ncol=md, byrow=TRUE)
      d.times.mat=matrix(rep(d.times,rep(s.num,md)),s.num)
      sk.le.dh=(sk.mat<=d.times.mat)
      
      fit.eta.Te=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=n1.h.mat, n2.h.mat=n2.h.mat, sk=sk, 
                       eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=1)
      AUC.Te=fp.func3(sk,tseq.eval,fit.eta.Te$par)
      AUC.Te.vec=as.vector(t(AUC.Te))
      par.Te=c(fit.eta.Te$convergence, fit.eta.Te$par)
      AUC.fd=rbind(AUC.fd, AUC.Te.vec)
      par.fd=rbind(par.fd, par.Te)
    }
    AUC.fd.avg=colMeans(AUC.fd)
    AUC.rep=rbind(AUC.rep, AUC.fd.avg)
    par.rep=rbind(par.rep, par.fd)
    Vseq.list[[rep.CV]]=Vseq
  }
  AUC.rep.CV=colMeans(AUC.rep)
  s.mat=matrix(rep(sk,each=t.num), ncol=t.num, byrow=TRUE)
  t.mat=matrix(rep(tseq.eval,rep(s.num,t.num)),s.num)
  slet.mat=1*(s.mat<=t.mat)
  slet.mat[slet.mat==0]=NA
  AUC.rep.CV=t(matrix(AUC.rep.CV, nrow=length(tseq.eval), ncol=s.num))*slet.mat #estimated AUC(s,t)
  theta=par.rep[,-1] #all repeated-CV thetas
  conv=par.rep[,1]
  ##########################################################################################################
  ##########################################################################################################
  #calculate the SE
  resap.AUC.mat=c() #save resampled AUC
  set.seed(seed)
  for( sap in 1:nsap)
  { 
    if (resample==1){ EXPwt=rexp(n, 1)} 
    if (resample==0){ EXPwt=rep(1, n) }
    AUC.rep.sap=c() #save repeated CV AUC
    for (rep.CV in 1:N.CV) {
      Vseq=Vseq.list[[rep.CV]]
      AUC.fd.sap=c()
      for (fd in 1:2) {
        EXPwt.Tr=EXPwt[Vseq==fd]
        EXPwt.Te=EXPwt[Vseq==3-fd]
        da.Tr=da.CV[Vseq==fd,] 
        marker.Tr=merge(da.Tr[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime", "eta"), X.names, Yt.names)], by=c("obs_id"))
        marker.Tr$logvt=log(marker.Tr$vtime+1)
        resap.coef.Tr=PC.Cox.wt(
          id = "obs_id",
          stime = "Y",
          status = "delta",
          measurement.time = "vtime",  ##survival time and measurement time must be on the same scale.
          predictors =c(c("logvt"), X.names, Yt.names),
          data = marker.Tr[marker.Tr$eta==1,],
          da.short=da.Tr, input_weight = EXPwt.Tr )$model.fit$coefficients[-1]
        da.Te=da.CV[Vseq==3-fd,] 
        marker.Te=merge(da.Te[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime", "eta"), X.names, Yt.names)], by=c("obs_id"))
        resap.coef.Tr.mat=t(matrix(rep(resap.coef.Tr, dim(marker.Te)[1]), nrow=n.X+n.Yt, ncol=dim(marker.Te)[1]))
        marker.Te$Zt=apply((resap.coef.Tr.mat)*marker.Te[,c(X.names, Yt.names)], 1,sum)
        
        ########################################Calculate resapAUC     ######################################
        #We use the da.Te as the survival data
        #We use marker.Te as the biomarker data
        
        n.Te=dim(da.Te)[1]
        d.times=sort(unique(da.Te$Y[da.Te$delta==1 & da.Te$Y<=d.per]))
        md=length(d.times)
        Ymat=matrix(rep(da.Te$Y,md),n.Te)
        dt.mat=matrix(rep(d.times,rep(n.Te,md)),n.Te)
        dt.mat.tb=dt.mat+0.00001; #tie breaking;
        n.h=apply((Ymat>dt.mat.tb),2,sum)
        death.mat=(Ymat==dt.mat)*da.Te$delta
        marker.Te$eta[marker.Te$eta==0]=0
        RS.data=marker.Te
        
        # Calculate the perturbed n1 and n2
        RS.seq.list=RS.mat.list=RS.d.list=RS.dmat.list=list() 
        Wn1.h.list=Wn2.h.list=list() 
        n1n2.wt=matrix(rep(EXPwt.Te,md),n.Te)
        for (i in 1:s.num)
        {
          RS.seq.list[[i]]=RS.data[RS.data$vtime==sk[i],]
          RS.mat.list[[i]]<-matrix(rep(RS.seq.list[[i]]$Zt,md),n.Te)
          RS.d.list[[i]]=apply(RS.mat.list[[i]]*death.mat,2,sum,na.rm=TRUE)
          RS.d.list[[i]][RS.d.list[[i]]==0]<-NA 
          RS.dmat.list[[i]]=matrix(rep(RS.d.list[[i]],rep(n.Te,md)),n.Te)
          Wn1.h.list[[i]]=apply((RS.dmat.list[[i]]>RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
          Wn2.h.list[[i]]=apply((RS.dmat.list[[i]]<=RS.mat.list[[i]])*(Ymat>dt.mat.tb)*n1n2.wt,2,sum,na.rm=TRUE)
        }
        Wn1.h.mat=Wn2.h.mat=c()  
        eta.sk.mat=c()
        for (i in 1:s.num)
        {
          Wn1.h.mat=rbind(Wn1.h.mat, Wn1.h.list[[i]])
          Wn2.h.mat=rbind(Wn2.h.mat, Wn2.h.list[[i]])
          eta.sk.mat=rbind(eta.sk.mat, RS.d.list[[i]])
        }  
        eta.sk.mat=(eta.sk.mat/eta.sk.mat)
        eta.sk.mat[is.na(eta.sk.mat)]=0
        sk.mat=matrix(rep(sk,each=md), ncol=md, byrow=TRUE)
        d.times.mat=matrix(rep(d.times,rep(s.num,md)),s.num)
        sk.le.dh=(sk.mat<=d.times.mat)
        
        # Extract the weights(Ewt.die) for subjects who experienced the event of interest
        da0<-cbind(da.Te$Y,da.Te$delta,EXPwt.Te)
        da0.sort=da0[order(da0[,1], decreasing = F),]
        Ewt.dseq= da0.sort[ da0.sort[,2]==1 &  da0.sort[,1]<=d.per,3]
        Ewt.die=matrix(rep(Ewt.dseq,rep(s.num,md)),s.num)
        
        #Computed the perturbed estimates of theta and AUC(s,t)
        resap.fit=optim(par=par0, fn=logL, d.times=d.times, n1.h.mat=Wn1.h.mat, n2.h.mat=Wn2.h.mat, sk=sk, 
                        eta.sk.mat=eta.sk.mat, sk.le.dh=sk.le.dh, control = list(maxit = 10000,reltol=1e-8),Ewt.die=Ewt.die)
        s.mat=matrix(rep(sk,each=t.num), ncol=t.num, byrow=TRUE)
        t.mat=matrix(rep(tseq.eval,rep(s.num,t.num)),s.num)
        slet.mat=1*(s.mat<=t.mat)
        slet.mat[slet.mat==0]=NA
        
        resap.AUC0=fp.func3(sk,tseq.eval,resap.fit$par)*slet.mat
        resap.AUC.vec= c(as.vector(t(resap.AUC0)),resap.fit$convergence)
        AUC.fd.sap=rbind(AUC.fd.sap, resap.AUC.vec)
      }
      AUC.fd.sap.avg=colMeans(AUC.fd.sap)
      AUC.rep.sap=rbind(AUC.rep.sap,AUC.fd.sap.avg)
    }
    AUC.rep.sap.avg=colMeans(AUC.rep.sap)
    resap.AUC.mat=rbind(resap.AUC.mat, AUC.rep.sap.avg)
    print(sap)
  }#resap.AUC.mat is the matrix combined 250 resampling AUC matrix
  
  resap.AUC.mat=as.data.frame(resap.AUC.mat)
  names(resap.AUC.mat)[ncol(resap.AUC.mat)]<-"cov1"
  resap.AUC.mat$minAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=min, na.rm=T)
  resap.AUC.mat$maxAUC=apply(resap.AUC.mat[,1:(t.num*s.num)], 1, FUN=max, na.rm=T)
  resap.AUC.mat$cov2=ifelse(resap.AUC.mat$minAUC<0.00001,1,0)
  resap.AUC.mat$cov3=ifelse(resap.AUC.mat$maxAUC>0.99999,1,0)
  resap.AUC.mat.cov=as.matrix(resap.AUC.mat[resap.AUC.mat$cov1==0&resap.AUC.mat$cov2==0&resap.AUC.mat$cov3==0,])
  RESAP.sd=apply(resap.AUC.mat.cov[,1:(t.num*s.num)], 2, FUN=sd)
  RESAP.sd.mat=t(matrix(RESAP.sd, t.num, s.num))
  return(list(AUC=AUC.rep.CV,conv=conv, SE=RESAP.sd.mat))
  }


da.sim.long <- read.csv("C:/Drive D/Thesis2/code for paper2/Code for regular visit/da.sim.long_reg.txt")
da.sim.short <- read.csv("C:/Drive D/Thesis2/code for paper2/Code for regular visit/da.sim.short_reg.txt")

AUC.st=main1.sub.func(da.short=da.sim.short, da.long=da.sim.long, sk=seq(0,2,0.25),
                  tseq.eval=seq(0.1,5,0.05), par0=c(0.5,0,0,0,0,0,0,0,0,0),  N.CV=3, resample=1, nsap=5, d.per=3.62)




