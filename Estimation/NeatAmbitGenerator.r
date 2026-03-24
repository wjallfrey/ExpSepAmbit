setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")
library(stats,lattice,MASS,parallel,ggplot2,reshape2,lattice,gridExtra,gganimate,pracma,lamW,ggpubr)

B<-matrix(c(1,0,0,1),2,2)

ambitfield3<-function(Lt,Lx,deltat,c,B,parameters,tol,depth){ #depth not actually included anywhere
  if(missing(tol)){
    tol<-10^(-3)
  }
  if(missing(depth)){
    depth=1
  }
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  g11=function(t,x,s,xi) g(t,x,s,xi,k=k11,mu=mu11,lambda=lambda11)
  g12=function(t,x,s,xi) 0
  g21=function(t,x,s,xi) g(t,x,s,xi,k=k21,mu=mu21,lambda=lambda21)
  g22=function(t,x,s,xi) g(t,x,s,xi,k=k22,mu=mu22,lambda=lambda22)
  
  deltax=c*deltat
  Ttop=Lt+Lx/2*(deltat/deltax)
  Xtop=0
  
  #Making depth of intergal - Should also extract lambdas and have max(1/2mu,1/(mu+lambda))
  intdep<-max(5,ceiling(1/(min(2*mu11,2*mu21,2*mu22,mu11+lambda11,mu21+lambda21,mu22+lambda22))*log(1/tol)))
  #Making random field
  fieldgridsize<-Lt/deltat+Lx/deltax+2*intdep/deltat
  Wvec<-mvrnorm(fieldgridsize^2,c(0,0),B*(2*deltat*deltax))
  W1=matrix(Wvec[,1],fieldgridsize,fieldgridsize)
  W2=matrix(Wvec[,2],fieldgridsize,fieldgridsize)
  ijtos<-function(i,j){
    return(Ttop-(i+j-1)*deltat)
  }
  ijtoxi<-function(i,j){
    return(Xtop+deltax*(j-i))
  }
  xistoi<-function(xi,s){
    return(round(((xi-Xtop-c*(Ttop-s))/(-deltax)+1)*0.5,2))
  }
  xistoj<-function(xi,s){
    return(round(((xi-Xtop+c*(Ttop-s))/(deltax)+1)*0.5,2))
  }
  integraly<-function(t,x,g,W){ #need (T-t)/deltat and (X-x)/deltax to be same parity (even)
    istar<-xistoi(x,t-deltat)
    jstar<-xistoj(x,t-deltat)
    #check istar and jstar are integers
    mu<--log(g(1,0,0,0)/g(0,0,0,0))
    lambda<--log(g(0,1,0,0)/g(0,0,0,0))
    if(istar%%1==0&&jstar%%1==0){
      int<-0
      for(i in (istar):(istar+intdep/deltax)){
        for(j in (jstar):(jstar+intdep/deltax)){
          if((i-istar)==(j-jstar)){
            int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j))*W[i,j]*(1-c*lambda*deltat/3+deltat^2*(mu^2+c^2*lambda^2)/12)
          }
          else{
            int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j))*W[i,j]*(1+deltat^2*(mu^2+c^2*lambda^2)/12)
          }
        }
      }
      return(as.numeric(int))
    }
    else{print("x,t value not on grid")}
  }
  #How to get Y(t+2*deltat,x) from Y(t,x)
  y2dtfromytx<-function(t,x,ytx,g,tdec,W){
    mu<--log(g(1,0,0,0)/g(0,0,0,0))
    lambda<--log(g(0,1,0,0)/g(0,0,0,0))
    if(t+2*deltat<=Lt){
      int<-ytx*exp(-2*tdec*deltat)
      iapex<-xistoi(x,t+deltat)
      japex<-xistoj(x,t+deltat)
      if(iapex%%1==0&&japex%%1==0){
        for(i in (iapex+1):(iapex+intdep/deltat)){
          int<-int+g(t+2*deltat,x,ijtos(i,japex),ijtoxi(i,japex))*W[i,japex]*(1+deltat^2*(mu^2+c^2*lambda^2)/12)
        }
        for(j in (japex+1):(japex+intdep/deltat)){
          int<-int+g(t+2*deltat,x,ijtos(iapex,j),ijtoxi(iapex,j))*W[iapex,j]*(1+deltat^2*(mu^2+c^2*lambda^2)/12)
        }
        int<-int+g(t+2*deltat,x,ijtos(iapex,japex),ijtoxi(iapex,japex))*W[iapex,japex]*(1-c*lambda*deltat/3+deltat^2*(mu^2+c^2*lambda^2)/12)
        return(as.numeric(int))
      }
      else{print("error")}
    }
    else{print("increment goes out of t range")}
  }
  mtox<-function(m){
    return(round(-Lx/2+(m-1)*2*deltax,2))
  }
  ntot<-function(n){
    return(round((n-1)*2*deltat,2))
  }
  #Compute field
  y11matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  y12matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  y21matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  y22matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  for(m in (1):(Lx/(2*deltax)+1)){
    y11matrix[1,m]<-as.numeric(integraly(0,mtox(m),g11,W1))
    #y12matrix[1,m]<-as.numeric(integraly(0,mtox(m),g12,W2)) ALL ZEROS
    y21matrix[1,m]<-as.numeric(integraly(0,mtox(m),g21,W1))
    y22matrix[1,m]<-as.numeric(integraly(0,mtox(m),g22,W2))
  }
  print("done first layer")
  for(n in 2:(Lt/(2*deltat)+1)){
    for(m in (1):(Lx/(2*deltax)+1)){
      y11matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y11matrix[n-1,m]),g11,mu11,W1))
      #y12matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y12matrix[n-1,m]),g12,0,W2)) ALL ZEROS
      y21matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y21matrix[n-1,m]),g21,mu21,W1))
      y22matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y22matrix[n-1,m]),g22,mu22,W2))
    }
  }
  print("done whole grid")
  y1matrix=y11matrix+y12matrix
  y2matrix=y21matrix+y22matrix
  
  return(list(y1=y1matrix,y2=y2matrix,Lx=Lx,Lt=Lt,c=c,deltat=deltat,parameters=parameters))
}

set.seed(123)
ambitgenerator1<-function(Lt,Lx,deltat,c,B,parameters,alpha,tol,W1,W2,getmoments){
  if((Lt/(2*deltat))%%1!=0&&(Lx/(2*c*deltat))%%1!=0){
    return("Incompatible parameters: need Nx,Nt to be integers")
  }
  g<-function(t,x,s,xi,k,mu,lambda){
    return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
  }
  if(missing(alpha)&&missing(tol)){
    print("tolerance set to 0.001")
    tol<-0.001
  }
  if(missing(getmoments)){getmoments=F}
  
  truncationerror<-function(alpha,mu,lambda,c){
      if(c*lambda==mu){
        return(exp(-alpha*deltat*mu)*(1+2*mu*(alpha*deltat/2)))
      }
      else{
        return(1-exp(-alpha*deltat*mu)*(1+2*mu*(1-exp(-alpha*deltat/2*(c*lambda-mu)))/(c*lambda-mu)))
      }
    }
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  
  deltax=c*deltat

  if(missing(alpha)){
    if((truncationerror(100,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),1)-tol)<0){
      alpha<-ceiling(uniroot(function(alphavar){truncationerror(alphavar,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),1)-tol},c(0,100))$root)
      print(paste("alpha set to",alpha))
    }
    else{print("tolerance too low")}
  }
  xisfieldsize<-Lx/deltax+(Lt)/deltat+alpha
  if(missing(W1)||missing(W2)){Wvec<-mvrnorm(xisfieldsize^2,c(0,0),B*(0.5*deltat*deltax))
    W1<-matrix(Wvec[,1],xisfieldsize,xisfieldsize)
    W2<-matrix(Wvec[,2],xisfieldsize,xisfieldsize)
  }
  #levelplot(W1)
  ijtos<-function(i,j){
    Lx/(2*c)+Lt-(i+j)*deltat/2+deltat/2
  }
  ijtoxi<-function(i,j){
    (j-i)*deltax/2
  }
  xistoi<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat-xi/deltax+1/2,3)
  }
  xistoj<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat+xi/deltax+1/2,3)
  }  
  
  integraldirect<-function(t,x,k,mu,lambda,W){
    int<-0
    iapex<-xistoi(x,t-deltat/2)
    japex<-xistoj(x,t-deltat/2)
    for(i in iapex:(iapex+alpha-1)){
      for(j in japex:(japex+alpha-1)){
        int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)*W[i,j]
      }
    }
    return(int)
  }
  
  
  
  y11matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y12matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y21matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y22matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  
  for(m in 1:(Lt/(deltat)+1)){
    for(n in 1:(Lx/(deltax)+1)){
      y11matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k11,mu11,lambda11,W1)
      y21matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k21,mu21,lambda21,W1)
      y22matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k22,mu22,lambda22,W2)
    }
  }
  
  
  
  
  y1matrix=y11matrix+y12matrix
  y2matrix=y21matrix+y22matrix
  
  #diagnostics
  diagnostics<-NULL
  if(getmoments){
    varY1<-var(as.vector(y1matrix))
    varY2<-var(as.vector(y2matrix))
    covY1Y2<-cov(as.vector(y1matrix),as.vector(y2matrix))
    diagnostics<-data.frame("mean Y1"=numeric(),"mean Y2"=numeric(),"var Y1"=numeric(),"var Y2"=numeric(),"cov Y1Y2"=numeric())
    diagnostics[1,]<-c(0,0,c*k11^2/(2*mu11*(mu11+c*lambda11)),c*k21^2/(2*mu21*(mu21+c*lambda21))+(c*k22^2/(2*mu22*(mu22+c*lambda22))),(2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21))))
    diagnostics[2,]<-c(mean(y1matrix),mean(y2matrix),varY1,varY2,covY1Y2)
    rownames(diagnostics)<-c("Expected","Computed")
  }
  return(list(y1=y1matrix,y2=y2matrix,Lx=Lx,Lt=Lt,c=c,deltat=deltat,parameters=parameters,moments=diagnostics))
}

ambitgenerator2<-function(Lt,Lx,deltat,c,B,parameters,alpha,tol,W1,W2,getmoments){
  if((Lt/(2*deltat))%%1!=0&&(Lx/(2*c*deltat))%%1!=0){
    return("Incompatible parameters: need Nx,Nt to be integers")
  }
  g<-function(t,x,s,xi,k,mu,lambda){
    return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
  }
  if(missing(alpha)&&missing(tol)){
    print("tolerance set to 0.001")
    tol<-0.001
  }
  if(missing(getmoments)){getmoments=F}
  
  truncationerror<-function(alpha,mu,lambda,c){
    if(c*lambda==mu){
      return(exp(-alpha*deltat*mu)*(1+2*mu*(alpha*deltat/2)))
    }
    else{
      return(1-exp(-alpha*deltat*mu)*(1+2*mu*(1-exp(-alpha*deltat/2*(c*lambda-mu)))/(c*lambda-mu)))
    }
  }
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  
  deltax=c*deltat
  
  if(missing(alpha)){
    if((truncationerror(100,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),c)>(1-tol))){
      alpha<-ceiling(uniroot(function(alphavar){truncationerror(alphavar,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),1)-1+tol},c(0,100))$root)
      print(paste("alpha set to",alpha))
    }
    else{print("tolerance too low")}
  }
  if(alpha%%1!=0||alpha<1){
    return("Need alpha integer greater than 1")
  }
  xisfieldsize<-Lx/deltax+(Lt)/deltat+alpha
  if(missing(W1)||missing(W2)){Wvec<-mvrnorm(xisfieldsize^2,c(0,0),B*(0.5*deltat*deltax))
  W1<-matrix(Wvec[,1],xisfieldsize,xisfieldsize)
  W2<-matrix(Wvec[,2],xisfieldsize,xisfieldsize)
  }
  #levelplot(W1)
  ijtos<-function(i,j){
    Lx/(2*c)+Lt-(i+j)*deltat/2+deltat/2
  }
  ijtoxi<-function(i,j){
    (j-i)*deltax/2
  }
  xistoi<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat-xi/deltax+1/2,3)
  }
  xistoj<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat+xi/deltax+1/2,3)
  }  
  
  integraldirect<-function(t,x,k,mu,lambda,W){
    int<-0
    variance<-0
    iapex<-xistoi(x,t-deltat/2)
    japex<-xistoj(x,t-deltat/2)
    for(i in iapex:(iapex+alpha-1)){
      for(j in japex:(japex+alpha-1)){
        int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)*W[i,j]
        variance<-variance+g(t,x,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)^2*(0.5*deltat*deltax)*(kappa(mu,lambda)+((i-iapex)==(j-japex))*kappaKD(mu,lambda))
      }
    }
    print(variance)
    return(int)
  }
  ydtfromytx<-function(t,x,ytx,k,mu,lambda,W){
    int<-ytx*exp(-mu*deltat)
    iapex<-xistoi(x,t-deltat/2)
    japex<-xistoj(x,t-deltat/2)
    for(i in (iapex+1):(iapex+alpha-1)){
      int<-int+g(t,x,ijtos(i,japex),ijtoxi(i,japex),k,mu,lambda)*W[i,japex]
    }
    for(j in (japex+1):(japex+alpha-1)){
      int<-int+g(t,x,ijtos(iapex,j),ijtoxi(iapex,j),k,mu,lambda)*W[iapex,j]
    }
    int<-int+g(t,x,ijtos(iapex,japex),ijtoxi(iapex,japex),k,mu,lambda)*W[iapex,japex]
    return(int)
  }
  
  y11matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y12matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y21matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y22matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  
  for(n in 1:(Lx/(deltax)+1)){
    y11matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k11,mu11,lambda11,W1)
    y21matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k21,mu21,lambda21,W1)
    y22matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k22,mu22,lambda22,W2)
  }
  for(m in 2:(Lt/(deltat)+1)){
    for(n in 1:(Lx/(deltax)+1)){
      y11matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y11matrix[m-1,n],k11,mu11,lambda11,W1)
      y21matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y21matrix[m-1,n],k21,mu21,lambda21,W1)
      y22matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y22matrix[m-1,n],k22,mu22,lambda22,W2)
    }
  }
  
  y1matrix=y11matrix+y12matrix
  y2matrix=y21matrix+y22matrix
  
  #diagnostics
  diagnostics<-NULL
  if(getmoments){
    varY1<-var(as.vector(y1matrix))
    varY2<-var(as.vector(y2matrix))
    covY1Y2<-cov(as.vector(y1matrix),as.vector(y2matrix))
    diagnostics<-data.frame("mean Y1"=numeric(),"mean Y2"=numeric(),"var Y1"=numeric(),"var Y2"=numeric(),"cov Y1Y2"=numeric())
    diagnostics[1,]<-c(0,0,c*k11^2/(2*mu11*(mu11+c*lambda11)),c*k21^2/(2*mu21*(mu21+c*lambda21))+(c*k22^2/(2*mu22*(mu22+c*lambda22))),(2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21))))
    diagnostics[2,]<-c(mean(y1matrix),mean(y2matrix),varY1,varY2,covY1Y2)
    rownames(diagnostics)<-c("Expected","Computed")
  }
  return(list(y1=y1matrix,y2=y2matrix,Lx=Lx,Lt=Lt,c=c,deltat=deltat,parameters=parameters,moments=diagnostics))
}

ambitgenerator3<-function(Lt,Lx,deltat,c,B,parameters,alpha,tol,W1,W2,getmoments){
  if((Lt/(2*deltat))%%1!=0&&(Lx/(2*c*deltat))%%1!=0){
    return("Incompatible parameters: need Nx,Nt to be integers")
  }
  g<-function(t,x,s,xi,k,mu,lambda){
    return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
  }
  if(missing(alpha)&&missing(tol)){
    print("tolerance set to 0.001")
    tol<-0.001
  }
  if(missing(getmoments)){getmoments=F}
  
  truncationerror<-function(alpha,mu,lambda,c){
    if(c*lambda==mu){
      return(exp(-alpha*deltat*mu)*(1+2*mu*(alpha*deltat/2)))
    }
    else{
      return(1-exp(-alpha*deltat*mu)*(1+2*mu*(1-exp(-alpha*deltat/2*(c*lambda-mu)))/(c*lambda-mu)))
    }
  }
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  
  deltax=c*deltat
  
  if(missing(alpha)){
    if((truncationerror(100,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),c)>(1-tol))){
      alpha<-ceiling(uniroot(function(alphavar){truncationerror(alphavar,min(mu11,mu21,mu22),min(lambda11,lambda21,lambda22),1)-1+tol},c(0,100))$root)
      print(paste("alpha set to",alpha))
    }
    else{print("tolerance too low")}
  }
  if(alpha%%1!=0||alpha<1){
    return("Need alpha integer greater than 1")
  }
  xisfieldsize<-Lx/deltax+(Lt)/deltat+alpha
  if(missing(W1)||missing(W2)){Wvec<-mvrnorm(xisfieldsize^2,c(0,0),B*(0.5*deltat*deltax))
  W1<-matrix(Wvec[,1],xisfieldsize,xisfieldsize)
  W2<-matrix(Wvec[,2],xisfieldsize,xisfieldsize)
  }
  #levelplot(W1)
  ijtos<-function(i,j){
    Lx/(2*c)+Lt-(i+j)*deltat/2+deltat/2
  }
  ijtoxi<-function(i,j){
    (j-i)*deltax/2
  }
  xistoi<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat-xi/deltax+1/2,3)
  }
  xistoj<-function(xi,s){
    round((Lx/(2*c)+Lt-s)/deltat+xi/deltax+1/2,3)
  }
  kappa<-function(mu,lambda){
    if(mu!=c*lambda){
      1/(2*deltat^2*(mu+c*lambda)*(mu-c*lambda))*(cosh(2*mu*deltat)-cosh(2*c*lambda*deltat))
    }
    else{
      1/(2*deltat^2*(mu+c*lambda))*(2*deltat*sinh(2*deltat*mu))
    }
  }
  kappaKD<-function(mu,lambda){
    if(mu!=c*lambda){
      1/(2*deltat^2*(mu+c*lambda)*(mu-c*lambda))*(-c*lambda/mu*sinh(2*mu*deltat)+sinh(2*c*lambda*deltat))
    }
    else{
      1/(2*deltat^2*(mu+c*lambda))*(-2*deltat*cosh(2*deltat*mu)+sinh(2*deltat*mu)/mu)
    }
  }
  integraldirect<-function(t,x,k,mu,lambda,W){
    int<-0
    iapex<-xistoi(x,t-deltat/2)
    japex<-xistoj(x,t-deltat/2)
    for(i in iapex:(iapex+alpha-1)){
      for(j in japex:(japex+alpha-1)){
        int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)*W[i,j]*(kappa(mu,lambda)+((i-iapex)==(j-japex))*kappaKD(mu,lambda))^0.5
      }
    }
    return(int)
  }
  ydtfromytx<-function(t,x,ytx,k,mu,lambda,W){
    int<-ytx*exp(-mu*deltat)
    iapex<-xistoi(x,t-deltat/2)
    japex<-xistoj(x,t-deltat/2)
    for(i in (iapex+1):(iapex+alpha-1)){
      int<-int+g(t,x,ijtos(i,japex),ijtoxi(i,japex),k,mu,lambda)*W[i,japex]*kappa(mu,lambda)^0.5
    }
    for(j in (japex+1):(japex+alpha-1)){
      int<-int+g(t,x,ijtos(iapex,j),ijtoxi(iapex,j),k,mu,lambda)*W[iapex,j]*kappa(mu,lambda)^0.5
    }
    int<-int+g(t,x,ijtos(iapex,japex),ijtoxi(iapex,japex),k,mu,lambda)*W[iapex,japex]*(kappa(mu,lambda)+kappaKD(mu,lambda))^0.5
    return(int)
  }
  
  y11matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y12matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y21matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  y22matrix=matrix(0,Lt/(deltat)+1,Lx/(deltax)+1)
  
  for(n in 1:(Lx/(deltax)+1)){
    y11matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k11,mu11,lambda11,W1)
    y21matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k21,mu21,lambda21,W1)
    y22matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k22,mu22,lambda22,W2)
  }
  for(m in 2:(Lt/(deltat)+1)){
    for(n in 1:(Lx/(deltax)+1)){
      y11matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y11matrix[m-1,n],k11,mu11,lambda11,W1)
      y21matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y21matrix[m-1,n],k21,mu21,lambda21,W1)
      y22matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y22matrix[m-1,n],k22,mu22,lambda22,W2)
    }
  }
  
  y1matrix=y11matrix+y12matrix
  y2matrix=y21matrix+y22matrix
  
  #diagnostics
  diagnostics<-NULL
  if(getmoments){
    varY1<-var(as.vector(y1matrix))
    varY2<-var(as.vector(y2matrix))
    covY1Y2<-cov(as.vector(y1matrix),as.vector(y2matrix))
    diagnostics<-data.frame("mean Y1"=numeric(),"mean Y2"=numeric(),"var Y1"=numeric(),"var Y2"=numeric(),"cov Y1Y2"=numeric())
    diagnostics[1,]<-c(0,0,c*k11^2/(2*mu11*(mu11+c*lambda11)),c*k21^2/(2*mu21*(mu21+c*lambda21))+(c*k22^2/(2*mu22*(mu22+c*lambda22))),(2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21))))
    diagnostics[2,]<-c(mean(y1matrix),mean(y2matrix),varY1,varY2,covY1Y2)
    rownames(diagnostics)<-c("Expected","Computed")
  }
  return(list(y1=y1matrix,y2=y2matrix,Lx=Lx,Lt=Lt,c=c,deltat=deltat,parameters=parameters,moments=diagnostics))
}

parasall04<-list(k11=0.4,k21=0.4,k22=0.4,mu11=0.4,mu21=0.4,mu22=0.4,lambda11=0.4,lambda21=0.4,lambda22=0.4)
parameters<-parasall04
set.seed(123)
Y<-ambitgenerator(100,200,1,1,B=diag(2),parasall04,tol=0.01,getmoments=T)
plotfield(Y)

Ygen2<-ambitgenerator2(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),parasall04,alpha=24,W1=W1Test,W2=W2Test,getmoments=T)
Ygen3<-ambitgenerator3(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),parasall04,alpha=24,W1=W1Test,W2=W2Test,getmoments=T)

paraset2<-list(k11=0.5,k21=0.3,k22=0.4,mu11=0.2,mu21=0.8,mu22=0.4,lambda11=0.3,lambda21=0.9,lambda22=0.5)
sd<-1
set.seed(sd)
Y2<-ambitgenerator2(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),paraset2,alpha=100,getmoments=T)
set.seed(sd)
Y3<-ambitgenerator3(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),paraset2,alpha=100,getmoments=T)
print(rbind(Y2[[8]][1,],Y2[[8]][2,],Y3[[8]][2,]))

momentsa2<-matrix(0,100,5)
momentsa3<-matrix(0,100,5)
for(i in 1:100){
  print(i)
  WvecTest<-mvrnorm(1000^2,c(0,0),diag(2)*(0.5*deltat*deltax))
  W1Test<-matrix(WvecTest[,1],1000,1000)
  W2Test<-matrix(WvecTest[,2],1000,1000)
  Y2<-ambitgenerator2(Lt=100,Lx=200,deltat=0.5,c=1,B=diag(2),parasall04,W1=W1Test,W2=W2Test,alpha=10,getmoments=T)
  Y3<-ambitgenerator3(Lt=100,Lx=200,deltat=0.5,c=1,B=diag(2),parasall04,W1=W1Test,W2=W2Test,alpha=10,getmoments=T)
  momentsa2[i,]<-as.numeric(Y2[[8]][2,])
  momentsa3[i,]<-as.numeric(Y3[[8]][2,])
}
boxplot(momentsa2[,1],momentsa2[,2],momentsa3[,1],momentsa3[,2])
boxplot(momentsa2[,3],momentsa3[,3])
boxplot(momentsa2[,4],momentsa3[,4])
boxplot(momentsa2[,5],momentsa3[,5])


Y<-ambitgenerator3(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),parasall04,alpha=24,getmoments=T)
as.numeric(Y[[8]][2,])

moments<-matrix(0,100,5)
for(i in 1:100){
  print(i)
  WvecTest<-mvrnorm(1000^2,c(0,0),diag(2)*(0.5*deltat*deltax))
  W1Test<-matrix(WvecTest[,1],1000,1000)
  W2Test<-matrix(WvecTest[,2],1000,1000)
  Y3<-ambitgenerator3(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),parasall04,W1=W1Test,W2=W2Test,alpha=24,getmoments=T)
  moments[i,]<-as.numeric(Y3[[8]][2,])
}



set.seed(123)
deltat<-1
c<-1
deltax<-c*deltat
WvecTest<-mvrnorm(1000^2,c(0,0),diag(2)*(0.5*deltat*deltax))
W1Test<-matrix(WvecTest[,1],1000,1000)
W2Test<-matrix(WvecTest[,2],1000,1000)
YTestseed123<-ambitgenerator(100,200,1,1,B=diag(2),parasall04,tol=0.01,W1=W1Test,W2=W2Test,getmoments=T)


W1constant<-matrix(1,1000,1000)
W2constant<-matrix(1,1000,1000)

Lx<-200
Lt<-100
deltat<-1
c<-1
deltax<-c*deltat
alpha<-5
B=diag(2)

Y<-ambitgenerator(Lt=100,Lx=200,deltat=1,c=1,B=diag(2),parasall04,tol=0.01,W1=W1constant*0.5,W2=W1constant*0.5)

W2constant<-matrix(0.5*deltat*deltax,xisfieldsize,xisfieldsize)

#Y<-ambitgenerator(100,200,1,1,B=diag(2),parasall04,alpha,W1constant,W2constant)

xistoi(0,0-deltat/2)
xistoj(0,0-deltat/2)
alpha<-10
grid=matrix(0,alpha,alpha)
xisfieldsize<-Lx/deltax+(Lt)/deltat+alpha
W1constant<-matrix(0.5*deltat*deltax,xisfieldsize,xisfieldsize)
for(i in 0:(alpha-1)){
  for(j in 0:(alpha-1)){
    grid[i+1,j+1]<-g(0,0,ijtos(201+i,201+j),ijtoxi(201+i,201+j),0.4,0.4,0.4)*W1constant[i+1,j+1]
  }
}
sum(grid)
2.5*f(alpha)

Y[[1]][1]




plotfield(Y)
