LASER <-function(nsample=length(z), X,z, X.target, m=c(6,8), centering='LP', coef.smooth='BIC', parallel=FALSE){
  X<-as.matrix(X)
  X.target<-matrix(X.target,1,ncol(X))
  n<-length(z)
  if(is.null(centering)){
    y<-z
    zmean<-rep(0,length(z))
  }else if(centering=='LP'){
    Tx<-eLP.poly(X,m[1])
    reg.dat=as.data.frame(cbind(z,Tx))
    if(ncol(Tx)>50){big.flag=TRUE}else{big.flag=FALSE}
    fit1 <- leaps::regsubsets(z~., data = reg.dat,intercept=TRUE,really.big=big.flag)
    id<-which.min(summary(fit1)$bic)
    coefi <- coef(fit1, id = id)
    zmean=cbind(rep(1,n),matrix(Tx[,names(coefi)[-1]],n,length(coefi)-1))%*%as.matrix(coefi)
    y<-z-zmean
  }else if(centering=='lm'){
	lmfit<-lm(z~X)
	zmean<-fitted(lmfit)
	y<-z-zmean
  }else if(centering=='spline' & ncol(X)==1){
	splinfit<-smooth.spline(X,z,df=8)
	zmean<-fitted(splinfit)
	y<-z-zmean
  }
  
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
  
  zmean.ind<-which(apply(X,1,function(x) all(x==X.target)))
  out<-list(data=y.sample+zmean[zmean.ind][1],
            LPcoef=Lcoef)
  return(out)
}