LPcden <-
function(X,y,X.test,m=c(6,8),method='BIC'){ 
  X.test<-as.matrix(X.test)
  X<-as.matrix(X)
  multivar=0
  if(ncol(X)>1){multivar=1}
  
  TX<-eLP.poly(X,m[1])
  ## need adjust to multivar X:
  #colnames(TX)<-paste('TX',1:m[1],sep='')
  Ty<-eLP.poly(y,m[2])
  colnames(Ty)<-paste('Ty',1:m[2],sep='')
  LPregression<-function(j,TX,Ty,X,X.test){
    reg.dat<-as.data.frame(cbind( Ty[,j],TX))
    colnames(reg.dat)[1]<-'Tyj'
    
    if(ncol(TX)<=1){
      frmla<-'Tyj~.-1'
    }else{
      if(ncol(TX)>50){big.flag=TRUE}else{big.flag=FALSE}
      fit1 <- leaps::regsubsets(Tyj~., data = reg.dat,intercept=FALSE,really.big=big.flag)
      
      id<-which.min(summary(fit1)$bic)
      coefi <- coef(fit1, id = id)
      
      frmla<-paste0('Tyj~',names(coefi)[1])
      if(length(coefi)>1){
        for(i in 2:length(coefi)){
          frmla<-paste0(frmla,'+',names(coefi)[i])
        }
      }
      frmla<-paste0(frmla,'-1')
    }
    mi<-rep(1,ncol(X)+1)
    Txapprox=matrix(0,nrow(X.test),ncol(TX))
    for(i in 1:ncol(X)){
      mi[i+1]<-min(m[1],length(unique(X[,i]))-1)
      Txi<-Predict.LP.poly(X[,i],as.matrix(TX[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]),X.test[,i])
      Txapprox[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]=as.matrix(Txi)
    }
    newdat=as.data.frame(Txapprox)
    colnames(newdat)<-colnames(TX)  
    
    lp.lm <- lm( as.formula(frmla),data=reg.dat)
    lp.pred<-predict(lp.lm,newdata=newdat,se.fit=TRUE)
    
    out<-rbind(lp.pred$fit,lp.pred$se.fit)
    
    return(out)
  }
  
  LP.coef<-matrix(0,nrow(X.test),m[2])
  colnames(LP.coef)<-paste('LP[',1:m[2],']',sep='')
  
  for(i in 1:nrow(X.test)){
    X.test0<-matrix(X.test[i,],nrow=1)
    LP.coef0<-sapply(1:m[2],LPregression,TX,Ty,X,X.test0)
    LP.coef[i,]<-LP.smooth(LP.coef0,n=length(y),method=method)
  }
  return(LP.coef)
}
