LASER <-
function(nsample=length(z), X,z, X.target, m=c(6,8), centering='LP', coef.smooth='BIC', parallel=FALSE){
  X<-as.matrix(X)
  X.target<-matrix(X.target,1,ncol(X))
  n<-length(z)
  zm.target<-0
  zmean<-rep(0,length(z))
  if(is.null(centering)){
    y<-z
  }else if(centering=='LP'){
    Tx<-eLP.poly(X,m[1])
    reg.dat=as.data.frame(cbind(z,Tx))
    if(ncol(Tx)>50){big.flag=TRUE}else{big.flag=FALSE}
    fit1 <- leaps::regsubsets(z~., data = reg.dat,intercept=TRUE,really.big=big.flag)
    id<-which.min(summary(fit1)$bic)
    coefi <- coef(fit1, id = id)
    if(length(names(coefi))<1){
      zmean<-0
      zm.target<-0
    }else{
      zmean=cbind(rep(1,n),matrix(Tx[,names(coefi)[-1]],n,length(coefi)-1))%*%as.matrix(coefi)
      y<-z-zmean
      ## predicting z mean at x.target:
      mi<-rep(1,ncol(X)+1)
      Txapprox=matrix(0,1,ncol(Tx))
      colnames(Txapprox)<-colnames(Tx)
      for(i in 1:ncol(X)){
        mi[i+1]<-min(m[1],length(unique(X[,i]))-1)
        Txi<-Predict.LP.poly(X[,i],as.matrix(Tx[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]),X.target[,i])
        Txapprox[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]=Txi
      }
      zm.target<-cbind(1,matrix(Txapprox[,names(coefi)[-1]],nrow=1))%*%as.matrix(coefi)
    }
  }else if(centering=='lm'){
    lmfit<-lm(z~X)
    zmean<-fitted(lmfit)
    y<-z-zmean
    zm.target<-matrix(c(1,X.target),nrow=1)%*%as.matrix(lmfit$coefficients)
  }else if(centering=='spline' & ncol(X)==1){
    x1<-as.numeric(X)
    splinfit<-smooth.spline(x1,z,df=8)
    zmean<-fitted(splinfit)
    y<-z-zmean
    xnew<-data.frame(x1=as.numeric(X.target))
    zm.target<-predict(splinfit,xnew)$y
  }
  zm.target<-as.numeric(zm.target)
  
  Lcoef<-LPcden(X,y,m,X.test=X.target,method=coef.smooth)
  if(sum(abs(Lcoef))==0){
    y.sample<-y
  }else{
    if(parallel==TRUE){
      numCores<-round(detectCores()/2)
      cl<-makeCluster(numCores)
    }else{
      cl<-NULL
    }
    y.sample<-g2l.sampler(nsample,LP.par=t(Lcoef),Y=y,clusters=cl)
    if(parallel==TRUE){
      stopCluster(cl)
    }
  }
  out<-list(data=y.sample+as.numeric(zm.target)[1],LPcoef=Lcoef)
  return(out)
}
