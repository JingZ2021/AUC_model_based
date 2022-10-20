#prepare the dataset
Backtofull.irreg=function(da.short, da.long)
{ 
  n=dim(da.short)[1]
  n.Yt=dim(da.long)[2]-4
  n.X=dim(da.short)[2]-3
  Yt.names=paste0("Y", seq(1:n.Yt))
  X.names =paste0("X", seq(1:n.X))
  colnames(da.short)=c(c("obs_id", "Y", "delta"), X.names)
  colnames(da.long) =c(c("obs_id", "vtime", "Y", "delta"), Yt.names)
  da.long1=merge(da.long[,c(c("obs_id", "vtime"),Yt.names)], da.short[,c(c("obs_id", "Y" , "delta"), X.names)], by=c("obs_id"))
  da.long1=da.long1[, c(c("obs_id", "Y", "delta", "vtime"), X.names, Yt.names )]
  da.short1=da.long1[!duplicated(da.long1$obs_id), c("obs_id", "delta", "Y")]
  return(list(da.short=da.short1,da.long=da.long1, X.names=X.names, Yt.names=Yt.names, n.X=n.X, n.Yt=n.Yt))}



#### get the list contains the S.j, X.j, and Z.j
Get_J_list<-function(da.long, death.long.list,n,c)
{
  risk.long.list=list()
  for (i in 1:length(death.long.list))
  { 
    
    s.i.curr=death.long.list[[i]]$S.i
    obs_id.i.curr=death.long.list[[i]]$obs_id.i
    names.ext=c("obs_id","vtime","Y","delta","Zt", "EXPwt")
    k.h=function(x,h){#biquadratic kernel
      k.h=15/16*(1-(x/h)^2)^2*(abs(x/h)<=1) }
    da.calc=da.long[da.long$obs_id!=obs_id.i.curr,]
    da.calc=da.calc[,names.ext]
    da.calc$ttime=abs(da.calc$vtime-s.i.curr)
    a.calc=by(da.calc,da.calc$obs_id, function(x) x[which.min(x$ttime), ] )
    da.ext=data.frame(t(matrix(unlist(a.calc),dim(da.calc)[2])))
    names(da.ext)=c("obs_id.j","S.j","X.j","delta","Z.j","EXPwt.j","ttime")
    hn=c*(n^(-1/3))
    da.ext$k.s=k.h(da.ext$S.j-s.i.curr, hn)
    risk.long.list[[i]]=da.ext
  }
  return(risk.long.list)
}

#### function 4 get the LogL
logL=function(theta,N12.list){
  AUC.list=fp.func3(N12.list,theta) ######d.times is t, #sk is s  
  AUC.list$ln_L= AUC.list$n1*log(AUC.list$AUC.st)+AUC.list$n2*log(1-AUC.list$AUC.st)
  logL=-sum(AUC.list$ln_L,na.rm=T)}

logL.resap=function(theta,N12.list){
  AUC.list=fp.func3.resap(N12.list,theta) ######d.times is t, #sk is s  
  AUC.list$ln_L= AUC.list$EXPwt.i*(AUC.list$n1*log(AUC.list$AUC.st)+AUC.list$n2*log(1-AUC.list$AUC.st))
  logL=-sum(AUC.list$ln_L,na.rm=T)}

fp.func3=function(ST.list, theta){
  ST.mat<-matrix(unlist(ST.list), ncol = 4, byrow = TRUE)
  theta.mat=t(matrix(rep(theta,length(ST.mat[,1])),ncol=length(ST.mat[,1])))
  ST.theta.mat=cbind(ST.mat, theta.mat)
  ST.theta.dat=as.data.frame(ST.theta.mat)
  names(ST.theta.dat)=c("n1", "n2", "d.i", "S.i", "theta1","theta2", "theta3"
                        ,"theta4", "theta5", "theta6", "theta7", "theta8", "theta9", "theta10")
  n1=ST.theta.dat$n1
  n2=ST.theta.dat$n2
  d.i=ST.theta.dat$d.i
  S.i=ST.theta.dat$S.i
  theta1=ST.theta.dat$theta1
  theta2=ST.theta.dat$theta2
  theta3=ST.theta.dat$theta3
  theta4=ST.theta.dat$theta4
  theta5=ST.theta.dat$theta5
  theta6=ST.theta.dat$theta6
  theta7=ST.theta.dat$theta7
  theta8=ST.theta.dat$theta8
  theta9=ST.theta.dat$theta9
  theta10=ST.theta.dat$theta10
  ST.theta.dat$teta=theta1*1+theta2*sqrt(S.i)+theta3*sqrt(d.i)+theta4*sqrt((S.i)^2)+theta5*sqrt((d.i)^2)+theta6*sqrt((S.i*d.i))+theta7*sqrt(((d.i)^3))+theta8*sqrt(((S.i)^3))+theta9*sqrt((S.i*((d.i)^2)))+theta10*sqrt(((S.i)^2)*d.i)
  
  ST.theta.dat$teta[ST.theta.dat$teta>=100]=100
  ST.theta.dat$teta[ST.theta.dat$teta<=-100]=-100
  ST.theta.dat$AUC.st=exp(ST.theta.dat$teta)/(1+exp(ST.theta.dat$teta))
  return(ST.theta.dat)} 

fp.func3.resap=function(ST.list, theta){
  ST.mat<-matrix(unlist(ST.list), ncol = 5, byrow = TRUE)
  theta.mat=t(matrix(rep(theta,length(ST.mat[,1])),ncol=length(ST.mat[,1])))
  ST.theta.mat=cbind(ST.mat, theta.mat)
  ST.theta.dat=as.data.frame(ST.theta.mat)
  names(ST.theta.dat)=c("n1", "n2", "d.i", "S.i", "EXPwt.i","theta1","theta2", "theta3"
                        ,"theta4", "theta5", "theta6", "theta7", "theta8", "theta9", "theta10")
  
  n1=ST.theta.dat$n1
  n2=ST.theta.dat$n2
  d.i=ST.theta.dat$d.i
  S.i=ST.theta.dat$S.i
  theta1=ST.theta.dat$theta1
  theta2=ST.theta.dat$theta2
  theta3=ST.theta.dat$theta3
  theta4=ST.theta.dat$theta4
  theta5=ST.theta.dat$theta5
  theta6=ST.theta.dat$theta6
  theta7=ST.theta.dat$theta7
  theta8=ST.theta.dat$theta8
  theta9=ST.theta.dat$theta9
  theta10=ST.theta.dat$theta10
  ST.theta.dat$teta=theta1*1+theta2*sqrt(S.i)+theta3*sqrt(d.i)+theta4*sqrt((S.i)^2)+theta5*sqrt((d.i)^2)+theta6*sqrt((S.i*d.i))+theta7*sqrt(((d.i)^3))+theta8*sqrt(((S.i)^3))+theta9*sqrt((S.i*((d.i)^2)))+theta10*sqrt(((S.i)^2)*d.i)
  
  ST.theta.dat$teta[ST.theta.dat$teta>=100]=100
  ST.theta.dat$teta[ST.theta.dat$teta<=-100]=-100
  ST.theta.dat$AUC.st=exp(ST.theta.dat$teta)/(1+exp(ST.theta.dat$teta))
  return(ST.theta.dat)}

fp.calc=function(s,t,theta){
  AUC.st.mat=c()
  for (s.curr in s)
  { 
    tsvec=cbind(1, sqrt(s.curr), sqrt(t), sqrt(s.curr^2), sqrt(t^2), sqrt(s.curr*t), sqrt(t^3), sqrt(s.curr^3), sqrt((t^2)*s.curr), sqrt(t*(s.curr^2)))
    teta=as.vector(as.matrix(tsvec)%*%theta)
    teta[teta>=100]=100
    teta[teta<=-100]=-100
    AUC.st=exp(teta)/(1+exp(teta))
    AUC.st.mat=rbind(AUC.st.mat, AUC.st)
  }
  return(AUC.st.mat)}


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


main1.sub.func.ireg=function(da.short, da.long, sk, tseq.eval, par0, N.CV, resample, nsap, d.per, c, seed=12345){
  s.num=length(sk)
  t.num=length(tseq.eval)
  n=dim(da.short)[1]
  
  data0=Backtofull.irreg(da.short=da.short, da.long=da.long)
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
      marker.Tr=merge(da.Tr[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime"), X.names, Yt.names)], by=c("obs_id"))
      marker.Tr$logvt=log(marker.Tr$vtime+1)
      pc.coef.Tr <-  PC.Cox(
        id = "obs_id",
        stime = "Y",
        status = "delta",
        measurement.time = "vtime",  ##survival time and measurement time must be on the same scale.
        predictors =c(c("logvt"),X.names, Yt.names),
        data = marker.Tr)$model.fit$coefficients[-1] 
      da.Te=da.CV[Vseq==3-fd,]
      marker.Te=merge(da.Te[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime"), X.names, Yt.names)], by=c("obs_id"))
      pc.coef.Tr.mat=t(matrix(rep(pc.coef.Tr, dim(marker.Te)[1] ), ncol =  dim(marker.Te)[1], nrow=n.X+n.Yt))
      marker.Te$Zt=apply(pc.coef.Tr.mat*marker.Te[,c(X.names,Yt.names)], 1, sum)
      
      ########################################Compute EST of AUC(s,t)######################################
      #We use the da.Te as the survival data
      #We use marker.Te as the biomarker data
      #step1
      ptm <- proc.time()
      #d.long is the list contains dh, Zi, Si, for each death people on their biomarker visit time
      n.Te=dim(da.Te)[1]
      da.Te$EXPwt=1
      marker.Te$EXPwt=1
      d.long=marker.Te[marker.Te$delta==1&marker.Te$Y<=d.per,c("obs_id", "Y", "Zt", "vtime", "EXPwt")]
      d.long.sort=d.long[order(d.long$Y,d.long$vtime),]
      md=dim(d.long.sort)[1]
      row.names(d.long.sort)=as.character(c(1:md))
      d.long.list=as.list(as.data.frame(t(as.matrix(d.long.sort))))
      d.long.list=lapply(d.long.list, function(x){x=data.frame(t(x))})
      d.long.colname=c("obs_id.i", "d.i", "Z.i", "S.i", "EXPwt.i")
      d.long.list=lapply(d.long.list, setNames, d.long.colname)
      proc.time()-ptm #0.25 second
      
      #step2 #c=0.5
      ptm <- proc.time()
      risk.long.list=Get_J_list(da.long = marker.Te, death.long.list = d.long.list, n=n.Te,c=c)
      CALC.list<-mapply(cbind,d.long.list,risk.long.list, SIMPLIFY = F)
      proc.time()-ptm #108.72 second
      
      #step3
      ptm <- proc.time()
      CALC.list=lapply(CALC.list, within, n1<-k.s*(Z.i>Z.j)*(X.j>d.i))
      CALC.list=lapply(CALC.list, within, n2<-k.s*(Z.i<=Z.j)*(X.j>d.i))
      CALC.list=lapply(CALC.list, within, n1.W<-k.s*(Z.i>Z.j)*(X.j>d.i)*EXPwt.j)
      CALC.list=lapply(CALC.list, within, n2.W<-k.s*(Z.i<=Z.j)*(X.j>d.i)*EXPwt.j)
      CALC.sum.list=lapply(CALC.list, colSums, na.rm=T)
      CALC.mean.list=lapply(CALC.list, colMeans, na.rm=T)
      n1n2.list=lapply(CALC.sum.list, `[`, c("n1", "n2"))
      n1n2.list=lapply(n1n2.list, function(x){x=data.frame(t(x))})
      DiSi.list=lapply(CALC.mean.list, `[`, c("d.i", "S.i"))
      DiSi.list=lapply(DiSi.list, function(x){x=data.frame(t(x))})
      N12.list=mapply(cbind,n1n2.list,DiSi.list, SIMPLIFY = F)
      
      n1n2W.list=lapply(CALC.sum.list, `[`, c("n1.W", "n2.W"))
      n1n2W.list=lapply(n1n2W.list, function(x){x=data.frame(t(x))})
      DiSi.W.list=lapply(CALC.mean.list, `[`, c("d.i", "S.i", "EXPwt.i"))
      DiSi.W.list=lapply(DiSi.W.list, function(x){x=data.frame(t(x))})
      N12.W.list=mapply(cbind,n1n2W.list,DiSi.W.list, SIMPLIFY = F)
      proc.time()-ptm #2.96 second
      fit.eta.Te=optim(par=par0, fn=logL, N12.list=N12.list, control = list(maxit = 10000,reltol=1e-8))
      AUC.Te=fp.calc(s=sk,t=tseq.eval,theta=fit.eta.Te$par)
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
    #resap.AUC.curr.CV.mat=c()
    for (rep.CV in 1:N.CV) {
      Vseq=Vseq.list[[rep.CV]]
      AUC.fd.sap=c()
      for (fd in 1:2){
        EXPwt.Tr=EXPwt[Vseq==fd]
        EXPwt.Te=EXPwt[Vseq==3-fd]
        da.Tr=da.CV[Vseq==fd,]
        marker.Tr=merge(da.Tr[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime"), X.names, Yt.names)], by=c("obs_id"))
        marker.Tr$logvt=log(marker.Tr$vtime+1)
        resap.coef.Tr=PC.Cox.wt(
          id = "obs_id",
          stime = "Y",
          status = "delta",
          measurement.time = "vtime",  ##survival time and measurement time must be on the same scale.
          predictors =c(c("logvt"), X.names, Yt.names),
          data = marker.Tr,
          da.short=da.Tr, input_weight = EXPwt.Tr )$model.fit$coefficients[-1]
        da.Te=da.CV[Vseq==3-fd,] 
        marker.Te=merge(da.Te[,c("obs_id", "Y", "delta")], marker.CV[, c(c("obs_id", "vtime"), X.names, Yt.names)], by=c("obs_id"))
        resap.coef.Tr.mat=t(matrix(rep(resap.coef.Tr, dim(marker.Te)[1]), nrow=n.X+n.Yt, ncol=dim(marker.Te)[1]))
        marker.Te$Zt=apply((resap.coef.Tr.mat)*marker.Te[,c(X.names, Yt.names)], 1,sum)
        
        ########################################Calculate resapAUC     ######################################
        #We use the da.Te as the survival data
        #We use marker.Te as the biomarker data
        n.Te=dim(da.Te)[1]  
        da.Te$EXPwt=EXPwt.Te
        EXPwt.Te.vec=da.Te[,c("obs_id", "EXPwt")]
        marker.Te<- merge(marker.Te,EXPwt.Te.vec,by=c("obs_id"))
        
        #step1
        ptm <- proc.time()
        #death.long is the list contains dh, Zi, Si, for each death people on their biomarker visit time
        d.long=marker.Te[marker.Te$delta==1&marker.Te$Y<=d.per,c("obs_id", "Y", "Zt", "vtime", "EXPwt")]
        d.long.sort=d.long[order(d.long$Y,d.long$vtime),]
        md=dim(d.long.sort)[1]
        row.names(d.long.sort)=as.character(c(1:md))
        d.long.list=as.list(as.data.frame(t(as.matrix(d.long.sort))))
        d.long.list=lapply(d.long.list, function(x){x=data.frame(t(x))})
        d.long.colname=c("obs_id.i", "d.i", "Z.i", "S.i", "EXPwt.i")
        d.long.list=lapply(d.long.list, setNames, d.long.colname)
        proc.time()-ptm #0.25 second
        
        #step2 #c=0.5
        ptm <- proc.time()
        risk.long.list=Get_J_list(da.long = marker.Te, death.long.list = d.long.list, n=n.Te,c=c)
        CALC.list<-mapply(cbind,d.long.list,risk.long.list, SIMPLIFY = F)
        proc.time()-ptm #108.72 second
        
        #step3
        ptm <- proc.time()
        CALC.list=lapply(CALC.list, within, n1<-k.s*(Z.i>Z.j)*(X.j>d.i))
        CALC.list=lapply(CALC.list, within, n2<-k.s*(Z.i<=Z.j)*(X.j>d.i))
        CALC.list=lapply(CALC.list, within, n1.W<-k.s*(Z.i>Z.j)*(X.j>d.i)*EXPwt.j)
        CALC.list=lapply(CALC.list, within, n2.W<-k.s*(Z.i<=Z.j)*(X.j>d.i)*EXPwt.j)
        CALC.sum.list=lapply(CALC.list, colSums, na.rm=T)
        CALC.mean.list=lapply(CALC.list, colMeans, na.rm=T)
        n1n2.list=lapply(CALC.sum.list, `[`, c("n1", "n2"))
        n1n2.list=lapply(n1n2.list, function(x){x=data.frame(t(x))})
        DiSi.list=lapply(CALC.mean.list, `[`, c("d.i", "S.i"))
        DiSi.list=lapply(DiSi.list, function(x){x=data.frame(t(x))})
        N12.list=mapply(cbind,n1n2.list,DiSi.list, SIMPLIFY = F)
        
        n1n2W.list=lapply(CALC.sum.list, `[`, c("n1.W", "n2.W"))
        n1n2W.list=lapply(n1n2W.list, function(x){x=data.frame(t(x))})
        DiSi.W.list=lapply(CALC.mean.list, `[`, c("d.i", "S.i", "EXPwt.i"))
        DiSi.W.list=lapply(DiSi.W.list, function(x){x=data.frame(t(x))})
        N12.W.list=mapply(cbind,n1n2W.list,DiSi.W.list, SIMPLIFY = F)
        proc.time()-ptm #2.96 second
        
        resap.fit=optim(par=par0, fn=logL.resap, N12.list=N12.W.list, control = list(maxit = 10000,reltol=1e-8))
        s.mat=matrix(rep(sk,each=t.num), ncol=t.num, byrow=TRUE)
        t.mat=matrix(rep(tseq.eval,rep(s.num,t.num)),s.num)
        slet.mat=1*(s.mat<=t.mat)
        slet.mat[slet.mat==0]=NA
        resap.AUC0=fp.calc(sk,tseq.eval,resap.fit$par)*slet.mat
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


da.sim.short <- read.csv("C:/Drive D/Thesis2/code for paper2/Code for irregular visit/da.sim.short_irreg.txt")
da.sim.long <- read.csv("C:/Drive D/Thesis2/code for paper2/Code for irregular visit/da.sim.long_irreg.txt")

AUC.st=main1.sub.func.ireg(da.short=da.sim.short,
                           da.long=da.sim.long,
                           sk=seq(0,2,0.25),
                           tseq.eval=seq(0.1,5,0.05),
                           par0=c(0.5,0,0,0,0,0,0,0,0,0),
                           N.CV=3,
                           resample=1,
                           nsap=4,
                           d.per=3.03,
                           c=1)
