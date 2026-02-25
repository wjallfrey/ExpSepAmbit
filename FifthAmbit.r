setwd("/Users/willallfrey/Documents/R")
library(MASS)

k11=1
lambda11=1
mu11=0.1
c=1
B=matrix(c(1,0,0,1),2,2)

g<-function(t,x,s,xi){
  return(k11*exp(-mu11*abs(t-s)-lambda11*abs(x-xi)))
}



#Area where we want to find Y(t,x):
Lt=100 #t in [0,Lt]
Lx=200 #x in [-Lx/2,Lx/2]
deltat=0.5
deltax=c*deltat

Ttop=Lt+Lx/2*(deltat/deltax)
Xtop=0

#
integraldepthcoef<-20
#Random field
fieldgridsize<-Lt/deltat+Lx/deltax+2*integraldepthcoef/deltat
Wvec<-mvrnorm(fieldgridsize^2,c(0,0),B*(2*deltat*deltax)) #check resolution multiplier correct
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

#integraldepthcoef<-max(5,ceiling(1/mu11*log(10^6)))
#Evaluating integral directly
integraly<-function(t,x){ #need (T-t)/deltat and (X-x)/deltax to be same parity (even)
  integraldepthcoef<-max(5,ceiling(1/(2*mu11)*log(10^6)),ceiling(1/(lambda11+mu11)*log(10^6)))
  istar<-xistoi(x,t-deltat)
  jstar<-xistoj(x,t-deltat)
  #check istar and jstar are integers
  if(istar%%1==0&&jstar%%1==0){
    int<-0
    for(i in (istar):(istar+integraldepthcoef/deltax)){
      for(j in (jstar):(jstar+integraldepthcoef/deltax)){
        #For direct variance calculation: int<-int+(1*exp(-1*abs(t-ijtos(i,j))-1*abs(x-ijtoxi(i,j))))^2*(2*deltat*deltax)
        int<-int+g(t,x,ijtos(i,j),ijtoxi(i,j))*W1[i,j]
      }
    }
    print(g(t,x,ijtos(istar+integraldepthcoef/deltax+1,jstar),ijtoxi(istar+integraldepthcoef/deltax+1,jstar)))
    print(g(t,x,ijtos(istar,jstar+integraldepthcoef/deltax+1),ijtoxi(istar,jstar+integraldepthcoef/deltax+1)))
    print(g(t,x,ijtos(istar+integraldepthcoef/deltax+1,jstar+integraldepthcoef/deltax+1),ijtoxi(istar+integraldepthcoef/deltax+1,jstar+integraldepthcoef/deltax+1)))
    return(as.numeric(int))
  }
  else{return("x,t value not on grid")}
}
#How to get Y(t+2*deltat,x) from Y(t,x)
y2dtfromytx<-function(t,x,ytx){
  if(t+2*deltat<=Lt){
    int<-ytx*exp(-2*mu11*deltat)
    iapex<-xistoi(x,t+deltat)
    japex<-xistoj(x,t+deltat)
    for(i in (iapex+1):(iapex+integraldepthcoef/deltax)){
      int<-int+g(t+2*deltat,x,ijtos(i,japex),ijtoxi(i,japex))*W1[i,japex]
    }
    for(j in (japex+1):(japex+integraldepthcoef/deltax)){
      int<-int+g(t+2*deltat,x,ijtos(iapex,j),ijtoxi(iapex,j))*W1[iapex,j]
    }
    int<-int+g(t+2*deltat,x,ijtos(iapex,japex),ijtoxi(iapex,japex))*W1[iapex,japex]
    return(as.numeric(int))
  }
  else{return("increment goes out of t range")}
}

mtox<-function(m){
  return(round(-Lx/2+(m-1)*2*deltax,2))
}
ntot<-function(n){
  return(round((n-1)*2*deltat,2))
}

generatefield<-function(){
  ymatrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  for(m in (1):(Lx/(2*deltax)+1)){
    ymatrix[1,m]<-as.numeric(intergaly(0,mtox(m)))
  }
  for(n in 2:(Lt/(2*deltat)+1)){
    for(m in (1):(Lx/(2*deltax)+1)){
      ymatrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(ymatrix[n-1,m])))
    }
  }
  return(ymatrix)
}

f<-generatefield()
heatmap(f,Rowv = NA,Colv=NA)



ambitfield<-function(Lt,Lx,deltat,c,B,intdep,g11,g12,g21,g22){
  deltax=c*deltat
  Ttop=Lt+Lx/2*(deltat/deltax)
  Xtop=0
  #Making random field
  fieldgridsize<-Lt/deltat+Lx/deltax+2*integraldepthcoef/deltat
  Wvec<-mvrnorm(fieldgridsize^2,c(0,0),B*(2*deltat*deltax)) #check resolution multiplier correct
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
  integraly<-function(t,x,g1,g2){ #need (T-t)/deltat and (X-x)/deltax to be same parity (even)
    istar<-xistoi(x,t-deltat)
    jstar<-xistoj(x,t-deltat)
    #check istar and jstar are integers
    if(istar%%1==0&&jstar%%1==0){
      int<-0
      for(i in (istar):(istar+intdep/deltax)){
        for(j in (jstar):(jstar+intdep/deltax)){
          #For direct variance calculation: int<-int+(1*exp(-1*abs(t-ijtos(i,j))-1*abs(x-ijtoxi(i,j))))^2*(2*deltat*deltax)
          int<-int+g1(t,x,ijtos(i,j),ijtoxi(i,j))*W1[i,j]+g2(t,x,ijtos(i,j),ijtoxi(i,j))*W2[i,j]
        }
      }
      return(as.numeric(int))
    }
    else{print("x,t value not on grid")}
  }
  #How to get Y(t+2*deltat,x) from Y(t,x)
  y2dtfromytx<-function(t,x,ytx,g1,g2){
    if(t+2*deltat<=Lt){
      int<-ytx*exp(-2*mu11*deltat) ####WRONG
      iapex<-xistoi(x,t+deltat)
      japex<-xistoj(x,t+deltat)
      for(i in (iapex+1):(iapex+intdep/deltax)){
        int<-int+g1(t+2*deltat,x,ijtos(i,japex),ijtoxi(i,japex))*W1[i,japex]+g2(t+2*deltat,x,ijtos(i,japex),ijtoxi(i,japex))*W2[i,japex]
      }
      for(j in (japex+1):(japex+intdep/deltax)){
        int<-int+g1(t+2*deltat,x,ijtos(iapex,j),ijtoxi(iapex,j))*W1[iapex,j]+g2(t+2*deltat,x,ijtos(iapex,j),ijtoxi(iapex,j))*W2[iapex,j]
      }
      int<-int+g1(t+2*deltat,x,ijtos(iapex,japex),ijtoxi(iapex,japex))*W1[iapex,japex]+g2(t+2*deltat,x,ijtos(iapex,japex),ijtoxi(iapex,japex))*W2[iapex,japex]
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
  y1matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  y2matrix=matrix(0,Lt/(2*deltat)+1,Lx/(2*deltax)+1)
  for(m in (1):(Lx/(2*deltax)+1)){
    y1matrix[1,m]<-as.numeric(integraly(0,mtox(m),g11,g12))
    y2matrix[1,m]<-as.numeric(integraly(0,mtox(m),g21,g22))
  }
  for(n in 2:(Lt/(2*deltat)+1)){
    for(m in (1):(Lx/(2*deltax)+1)){
      y1matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y1matrix[n-1,m]),g11,g12))
      y2matrix[n,m]<-as.numeric(y2dtfromytx(ntot(n-1),mtox(m),(y2matrix[n-1,m]),g21,g22))
    }
  }
  return(list(y1matrix,y2matrix))
}

k11=1
k21=1
k22=1
mu11=1
mu21=1
mu22=1
lambda11=1
lambda21=1
lambda22=1


Y<-ambitfield(Lt=100,Lx=200,deltat=1,c=1,B=matrix(c(1,0,0,1),2,2),intdep=5,
           g11=function(t,x,s,xi) g(t,x,s,xi,k11,mu11,lambda11),
           g12=function(t,x,s,xi) 0,
           g21=function(t,x,s,xi) g(t,x,s,xi,k21,mu21,lambda21),
           g22=function(t,x,s,xi) g(t,x,s,xi,k22,mu22,lambda22))
y1<-Y[[1]]
y2<-Y[[2]]

heatmap(y1,Rowv = NA,Colv=NA)
heatmap(y2,Rowv = NA,Colv=NA)
