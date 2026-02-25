setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")
library(lattice)
library(MASS)

ambitfield<-function(Lt,Lx,deltat,c,B,tol,g11,g12,g21,g22){
  deltax=c*deltat
  Ttop=Lt+Lx/2*(deltat/deltax)
  Xtop=0
  #Extracting mu from functions g
  mu11<- -log(g11(1,0,0,0)/g11(0,0,0,0))
  mu21<- -log(g21(1,0,0,0)/g21(0,0,0,0))
  mu22<- -log(g22(1,0,0,0)/g22(0,0,0,0))
  lambda11<- -log(g11(0,0,0,1)/g11(0,0,0,0))
  lambda21<- -log(g21(0,0,0,1)/g21(0,0,0,0))
  lambda22<- -log(g22(0,0,0,1)/g22(0,0,0,0))
  #Making depth of intergal - Should also extract lambdas and have max(1/2mu,1/(mu+lambda))
  intdep<-max(5,ceiling(1/(min(2*mu11,2*mu21,2*mu22,mu11+lambda11,mu21+lambda21,mu22+lambda22))*log(1/tol)))
  #Making random field
  fieldgridsize<-Lt/deltat+Lx/deltax+2*intdep/deltat
  Wvec<-mvrnorm(fieldgridsize^2,c(0,0),B*(2*deltat*deltax)) #check resolution multiplier correct
  W1=matrix(Wvec[,1],fieldgridsize,fieldgridsize)
  W2=matrix(Wvec[,2],fieldgridsize,fieldgridsize)
  print(c(intdep,fieldgridsize))
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
    if(istar%%1==0&&jstar%%1==0){
      int<-0
      for(i in (istar):(istar+intdep/deltax)){
        for(j in (jstar):(jstar+intdep/deltax)){
          #For direct variance calculation: int<-int+(1*exp(-1*abs(t-ijtos(i,j))-1*abs(x-ijtoxi(i,j))))^2*(2*deltat*deltax)
          int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j))*W[i,j]
        }
      }
      return(as.numeric(int))
    }
    else{print("x,t value not on grid")}
  }
  #How to get Y(t+2*deltat,x) from Y(t,x)
  y2dtfromytx<-function(t,x,ytx,g,tdec,W){
    if(t+2*deltat<=Lt){
      int<-ytx*exp(-2*tdec*deltat)
      iapex<-xistoi(x,t+deltat)
      japex<-xistoj(x,t+deltat)
      for(i in (iapex+1):(iapex+intdep/deltat)){
        int<-int+g(t+2*deltat,x,ijtos(i,japex),ijtoxi(i,japex))*W[i,japex]
      }
      for(j in (japex+1):(japex+intdep/deltat)){
        int<-int+g(t+2*deltat,x,ijtos(iapex,j),ijtoxi(iapex,j))*W[iapex,j]
      }
      int<-int+g(t+2*deltat,x,ijtos(iapex,japex),ijtoxi(iapex,japex))*W[iapex,japex]
      return(as.numeric(int))
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
  
  return(list(y1matrix,y2matrix,Lx,Lt,c,deltat))
}

g<-function(t,x,s,xi,k,mu,lambda){
  return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
}

Y<-ambitfield(Lt=100,Lx=200,deltat=0.2,c=1,B=matrix(c(1,0,0,1),2,2),tol=10^(-6),
              g11=function(t,x,s,xi) g(t,x,s,xi,k=1,mu=0.3,lambda=1),
              g12=function(t,x,s,xi) 0,
              g21=function(t,x,s,xi) g(t,x,s,xi,k=1,mu=0.6,lambda=1),
              g22=function(t,x,s,xi) g(t,x,s,xi,k=1,mu=0.4,lambda=1))

saveRDS(Y,file=paste("Yfieldtests",1,".rds",sep=""))
