g2l.infer.boot <-
function(X,z,m=c(6,8),X.test=NULL,alpha=.1,B=100,nsample=length(z), fdr.curve.approx=c("direct", "indirect"),
                          ngrid=2000,centering='LP', locfdr.df=10,fdr.th.fixed=NULL, coef.smooth='BIC',
                          rel.null=c('custom','th'),parallel=FALSE){
  X<-as.matrix(X)
  
  n<-length(z)
  X.axe<-unique(X)
  macro.flag<-FALSE
  if(is.null(X.test)){X.test<-X.axe; macro.flag<-TRUE}
  
  #centralize
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
  }else if(centering=='spline'& ncol(X)==1){
    splinfit<-smooth.spline(X,z,df=8)
    zmean<-fitted(splinfit)
    y<-z-zmean
  }
  
  bre<-200;df<-locfdr.df
  if(rel.null=='th'){nulltype<-0}else{nulltype<-1}
  
  w.full<-locfdr::locfdr(y,bre=bre,df=df,nulltype=nulltype,plot=0)
  fdr.pool<-approxfun(y,w.full$fdr,method='linear',rule=2)
  p0.pool<-w.full$fp0[2*nulltype+1,]
  if(p0.pool[3]>1){p0.pool[3]<-1}
  
  Ty<-eLP.poly(y,m[2])
  
  lpfdr.th<-list()
  lpfdr.fp0<-data.frame(mu=rep(0,nrow(X.test)),sd=rep(0,nrow(X.test)),pi0=rep(0,nrow(X.test)))
  
  lpfdr.th$ti<-lpfdr.th$u<-lpfdr.th$l<-rep(NA,nrow(X.test))
  lpfdr.grid<-lpfdr<-matrix(0,nrow(X.test),ngrid)
  if(macro.flag==FALSE){
    lpfdr<-array(0,c(nrow(X.test),B,ngrid))
  }
  
  
  simu_tags<-rep(0,length(z))
  z.fdr.tab<-data.frame(X=X,z=z,prob_null=0)
  
  if(parallel==TRUE){
    numCores<-round(detectCores()/2)
    cl<-makeCluster(numCores)
  }else{
    cl<-NULL
  }
  
  pb<-txtProgressBar(min = 0, max = nrow(X.test)*B, initial = 0, style = 3)
  setTxtProgressBar(pb,0)
  
  for(i in 1:nrow(X.test)){
    
    iX.test<-X.test[i,]
    ind0<-which.min(apply(X.axe,1,function(x) sum(abs(x-iX.test))))
    x0<-matrix(X.axe[ind0,],nrow=length(ind0))
    zmean.ind<-which(apply(X,1,function(x) all(x==x0)))
    
    Lcoef<-LPcden(X,y,m,X.test=x0,method=coef.smooth)
    
    if(sum(abs(Lcoef))==0){
      y.sample<-y
      ygrid<-seq(min(y.sample),max(y.sample),length.out=ngrid)
      lpfdr0<-fdr.pool
      parms0<-p0.pool
      lpfdrmat0<-matrix(1,B,1)%*%matrix(lpfdr0(ygrid),1,ngrid)
    }else{
      lpfdrmat0<-matrix(0,B,ngrid)
      for(iter in 1:B){  
        tryCatch({
          y.sample<-g2l.sampler(nsample,LP.par=t(Lcoef),Y=y,clusters=cl)
          ygrid<-seq(min(y.sample,y[zmean.ind]),max(y.sample,y[zmean.ind]),length.out=ngrid)
          w0 <- locfdr::locfdr(y.sample,bre=bre,df=df,nulltype=nulltype,plot=0)
          parms0<-w0$fp0[2*nulltype+1,]
          if(fdr.curve.approx=='direct'){
            lpfdr_iter<-approxfun(y.sample,w0$fdr,method='linear',rule=2)
          }else if(fdr.curve.approx=='indirect'){
            if(parms0[3]>1){parms0[3]<-1}
            fdr.y0<-fdr.pool(y)*(parms0[3]/p0.pool[3])*
              (dnorm(y,mean=parms0[1],sd=parms0[2])/
                 dnorm(y,mean=p0.pool[1],sd=p0.pool[2]))/
              (1+Lcoef%*%t(Ty))
            fdr.y0[which(fdr.y0>=1)]<-1
            lpfdr_iter<-approxfun(y,fdr.y0,method='linear',rule=2)
          }
          lpfdrmat0[iter,]<-lpfdr_iter(ygrid)
          setTxtProgressBar(pb,B*(i-1)+iter)
        },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      }
      lpfdr0<-approxfun(ygrid,colSums(lpfdrmat0)/B,method='linear',rule=2)
    }
    zgrid<-ygrid+zmean[zmean.ind[1]]
    x0fdr<-lpfdr0(ygrid)
    if(macro.flag==TRUE){
      lpfdr[i,]<-x0fdr
    }else{
      lpfdr[i,,]<-lpfdrmat0
    }
    
    zfdr<-lpfdr0(y[zmean.ind])
    zfdr[zfdr>1]<-1
    zfdr[zfdr<0]<-0
    
    lpfdr.grid[i,]<-zgrid
    
    z.fdr.tab$prob_null[zmean.ind]<-zfdr
    lpfdr.fp0$mu[i]<-parms0[1]+zmean[zmean.ind[1]]
    lpfdr.fp0$sd[i]<-parms0[2]
    lpfdr.fp0$pi0[i]<-parms0[3]
    
    #threshold finding: (control alpha level for each x)
    if(is.null(fdr.th.fixed)){
      th_out<-fdr.thresh(y.sample,lpfdr0,alpha=alpha)
      th.fdr<-min(th_out$th.fdr,2*alpha)
    }else{
      th.fdr<-fdr.th.fixed
    }
    lpfdr.th$ti[i]<-th.fdr
    
    #finding upper limit:
    l0<--Inf;u0<-Inf
    zgrid.revind<-order(zgrid,decreasing =TRUE)
    for(ind in 1:length(zgrid)){
      if(is.na(th.fdr)){
        u0<-max(zgrid)
        break
      }else if(x0fdr[zgrid.revind[ind]]>th.fdr){
        if(ind>1){
          u0<-zgrid[zgrid.revind[ind-1]]
        }
        break
      }
    }
    #finding lower limit:
    for(ind in 1:length(zgrid)){
      if(is.na(th.fdr)){
        l0<-min(zgrid)
        break
      }else if(x0fdr[ind]>th.fdr){
        if(ind>1){
          l0<-zgrid[ind-1]
        }
        break
      }
    }
    if(!is.infinite(u0)){lpfdr.th$u[i]<-u0}
    if(!is.infinite(l0)){lpfdr.th$l[i]<-l0}
    
    simu_tags[zmean.ind]=(z[zmean.ind]<l0)+(z[zmean.ind]>u0)
  }
  
  
  if(parallel==TRUE){stopCluster(cl)}
  
  
  out<-list(data=list(X=X,z=z,X.test=X.test,z.means=zmean,rej=simu_tags),
            fdr.z=z.fdr.tab,
            lp.m=m,
            lpfdr=lpfdr,
            z.grid=lpfdr.grid,
            lpfdr.fp0=lpfdr.fp0,
            lpfdr.th=lpfdr.th)
  
  
  return(out)
  
}
