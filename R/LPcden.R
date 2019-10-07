LPcden <-
function(X,y,X.test,m=c(6,8),method='BIC'){ 
  X<-as.matrix(X)
  X.test<-matrix(X.test,ncol=ncol(X))
  multivar=0
  if(ncol(X)>1){multivar=1}
  
  Tx<-eLP.poly(X,m[1])
  ## need adjust to multivar X:
  #colnames(Tx)<-paste('Tx',1:m[1],sep='')
  Ty<-eLP.poly(y,m[2])
  colnames(Ty)<-paste('Ty',1:m[2],sep='')
  
  LP.coef<-matrix(0,nrow(X.test),m[2])
  colnames(LP.coef)<-paste('LP[',1:m[2],']',sep='')
  for(i in 1:nrow(X.test)){
      X.test0<-matrix(X.test[i,],nrow=1)
      LP.coef0<-sapply(1:m[2],LPregression,Tx,Ty,X,X.test0,m)
      if(!is.null(method)){
        LP.coef[i,]<-LP.smooth(LP.coef0,n=length(y),method=method)
      }else{
        LP.coef[i,]<-LP.coef0[1,]
      }
    }
  return(LP.coef)
}
