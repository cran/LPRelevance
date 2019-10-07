LPregression <-
function(j,Tx,Ty,X,X.test,m){
  reg.dat<-as.data.frame(cbind( Ty[,j],Tx))
  colnames(reg.dat)[1]<-'Tyj'
  
  if(ncol(Tx)<=1){
    frmla<-'Tyj~.-1'
  }else{
    if(ncol(Tx)>50){big.flag=TRUE}else{big.flag=FALSE}
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
  #mi<-rep(1,ncol(X)+1)
  #Txapprox=matrix(0,nrow(X.test),ncol(Tx))
  #for(i in 1:ncol(X)){
  #  mi[i+1]<-min(m[1],length(unique(X[,i]))-1)
  #  Txi<-Predict.LP.poly(X[,i],as.matrix(Tx[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]),X.test[,i])
  #  Txapprox[,sum(mi[1:i]):(sum(mi[1:(i+1)])-1)]=as.matrix(Txi)
  #}
  #newdat=as.data.frame(Txapprox)
  #colnames(newdat)<-colnames(Tx) 
  
  newdat<-eLP.poly.predict(X,Tx,X.test,m[1])
  
  lp.lm <- lm( as.formula(frmla),data=reg.dat)
  lp.pred<-predict(lp.lm,newdata=newdat,se.fit=TRUE)
  
  out<-rbind(lp.pred$fit,lp.pred$se.fit)
  
  return(out)
}
