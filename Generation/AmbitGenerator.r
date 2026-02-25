setwd("/Users/willallfrey/Documents/R")
library(lattice)
library(ggplot2)
library(plotly)

ambitfield<-function(g11,g21,g22,B,c,Lt,Lx,buffert,bufferx,deltat){
  g12<-function(t,x,s,xi){return(0)}
  Lxi=Lx+2*bufferx
  Ls=Lt+buffert
  deltax=c*deltat
  Ns=Ls/deltat
  Nxi=Lxi/deltax
  Nt=Lt/deltat
  Nx=Lx/deltax
  #Construction of the Gaussian random field:
  Wvec<-mvrnorm((Nxi+1)*(Ns+1),c(0,0),B*(deltat*deltax)) #check resolution multiplier correct
  W1=matrix(Wvec[,1],Nxi+1,Ns+1)
  W2=matrix(Wvec[,2],Nxi+1,Ns+1)
  #Functions to change from index space to xi-s space:
  jtos<-function(j){
    if(j>=1&&j<=(Ns+1)){
      return(-buffert+(j-1)*Ls/Ns)
    }
    else{
      print("j out of range")
    }
  }
  itoxi<-function(i){
    if(1<=i&&i<=(Nxi+1)){
      return(-Lxi/2+(i-1)*Lxi/Nxi)
    }
    else{
      print("i out of range")
    }
  }
  #Full integral method:
  integraly<-function(t,x,g1,g2){
    int<-0
    for(j in 1:(Ns/Ls*(t+buffert)+1)){
      for(i in (max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)):(min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j)))))){
        int<-int+g1(t,x,jtos(j),itoxi(i))*W1[i,j]+g2(t,x,jtos(j),itoxi(i))*W2[i,j]
      }
    }
    return(int)
  }
  #Function to find Y(t+deltat,x) from Y(t,x) here written ytx:
  ytplusfromyt<-function(t,x,ytx,g1,g2){
    if(t+deltat<=Lt){
      output=exp(-mu11*deltat)*ytx #Unique to this exponential form of g
      for(j in 1:(Ns/Ls*(t+buffert)+1)){
        iplus<-Nxi/Lxi*(x+Lxi/2+c*(t+deltat-jtos(j)))+1
        iminus<-Nxi/Lxi*(x+Lxi/2-c*(t+deltat-jtos(j)))+1
        if(iminus>=1){
          output<-output+g1(t+deltat,x,jtos(j),itoxi(iminus))*W1[iminus,j]+g2(t+deltat,x,jtos(j),itoxi(iminus))*W2[iminus,j]
        }
        if(iplus<=(Nxi+1)){
          output<-output+g1(t+deltat,x,jtos(j),itoxi(iplus))*W1[iplus,j]+g2(t+deltat,x,jtos(j),itoxi(iplus))*W2[iplus,j]
        }
      }
      apexj<-(Ns/Ls*(t+deltat+buffert)+1)
      apexi<-Nxi/Lxi*(x+Lxi/2)+1
      output<-output+g1(t+deltat,x,jtos(apexj),itoxi(apexi))*W1[apexi,apexj]+g2(t+deltat,x,jtos(apexj),itoxi(apexi))*W2[apexi,apexj]
      return(output)
    }
    else{
      print("increment goes out of range")
    }
  }
  #index space to x-t space functions:
  ntox<-function(n){
    if(1<=n&&n<=Nx+1){
      return(-Lx/2+(n-1)*Lx/Nx)
    }
    else{
      print("n out of range")
    }
  }
  mtot<-function(m){
    if(1<=m&&m<=Nt+1){
      return((m-1)*Lt/Nt)
    }
    else{
      print("m out of range")
    }
  }
  #Construction of output matrix and its calculation:
  y1matrix=matrix(0,Nt+1,Nx+1)
  y2matrix=matrix(0,Nt+1,Nx+1)
  for(n in 1:(Nx+1)){
    y1matrix[1,n]=integraly(0,ntox(n),g11,g12)
    y2matrix[1,n]=integraly(0,ntox(n),g21,g22)
  }
  for(m in 2:(Nt+1)){
    for(n in 1:(Nx+1)){
      y1matrix[m,n]=ytplusfromyt(mtot(m-1),ntox(n),y1matrix[m-1,n],g11,g12)
      y2matrix[m,n]=ytplusfromyt(mtot(m-1),ntox(n),y2matrix[m-1,n],g21,g22)
    }
  }
  return(list(y1matrix,y2matrix))
}


g<-function(t,x,s,xi,k,mu,lambda){
  return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
}

k11=1
lambda11=1
mu11=1
k21=1
lambda21=1
mu21=1
k22=1
lambda22=1
mu22=1
c=1
B=matrix(c(1,0,0,1),2,2)
deltat=0.1
Lt=100 #t in [0,Lt]
Lx=200 #x in [-Lx/2,Lx/2]
bufferx=15 #Should change depending on lambda
buffert=15 #Should change depending on mu 

#deltax=c*deltat
#Nt=Lt/deltat
#Nx=Lx/deltax


fieldlist<-ambitfield(function(t,x,s,xi) g(t,x,s,xi,k11,mu11,lambda11),
                      function(t,x,s,xi) g(t,x,s,xi,k21,mu21,lambda21),
                      function(t,x,s,xi) g(t,x,s,xi,k22,mu22,lambda22),
                      B,c,Lt=100,Lx=200,buffert,bufferx,deltat=1)

y1<-fieldlist[[1]]
y2<-fieldlist[[2]]
#highresy1<-y1
#highresy2<-y2
heatmap(highresy1,Rowv = NA,Colv=NA)
heatmap(highresy2,Rowv = NA,Colv=NA)
var(as.vector(y1))

deltat=1
deltax=c*deltat
Nt=Lt/deltat
Nx=Lx/deltax
Xdomain<-((-Nx/2):(Nx/2))*(Lx/Nx)
Tdomain<-((0:Nt))*(Lt/Nt)


grid<-expand.grid(x = Xdomain,t = Tdomain)
grid$y1vec<-as.vector(y1)
grid$y2vec<-as.vector(y2)
levelplot(y1vec~x*t,grid)

#theoretical variance
Lxi=Lx+2*bufferx
Ls=Lt+buffert
Ns=Ls/deltat
Nxi=Lxi/deltax

varfun<-function(t,x){
  int<-0
  for(j in 1:(Ns/Ls*(t+buffert))){
    for(i in (max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)+1):(min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j))))-1)){
      int<-int+(1*exp(-1*abs(t-jtos(j))-1*abs(x-itoxi(i))))^2*deltat*deltax*1
    }
    iplus<-min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j))))
    iminus<-max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)
    print((1*exp(-1*abs(t-jtos(j))-1*abs(x-itoxi(iplus))))^2*deltat*deltax)
    int<-int+(1*exp(-1*abs(t-jtos(j))-1*abs(x-itoxi(iplus))))^2*deltat*deltax*0.316
    int<-int+(1*exp(-1*abs(t-jtos(j))-1*abs(x-itoxi(iminus))))^2*deltat*deltax*0.316
  }
  japex<-Ns/Ls*(t+buffert)+1
  iapex<-Nxi/Lxi*(x+Lxi/2)+1
  print((1*exp(-1*abs(t-jtos(japex))-1*abs(x-itoxi(iapex))))^2)
  int<-int+(1*exp(-1*abs(t-jtos(japex))-1*abs(x-itoxi(iapex))))^2*deltat*deltax*0.1
  return(int)
}
#varfun(0,0,function(t,x,s,xi) g(t,x,s,xi,k11,mu11,lambda11),function(t,x,s,xi) g(t,x,s,xi,0,mu21,lambda21))

varfun(0,0)
varfun2(0,0)

varfun2<-function(t,x){
  int<-0
  for(j in 1:(Ns/Ls*(t+buffert))+1){
    for(i in (max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)):(min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j)))))){
      int<-int+(1*exp(-1*abs(t-jtos(j))-1*abs(x-itoxi(i))))^2*deltat*deltax
    }
  }
  return(int)
}

#other stuff
integralySPLIT<-function(t,x,g1,g2){
  int<-0
  for(j in 1:(Ns/Ls*(t+buffert)+1)){
    for(i in (max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)):(min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j)))))){
      int<-int+g1(t,x,jtos(j),itoxi(i))*W1[i,j]+g2(t,x,jtos(j),itoxi(i))*W2[i,j]
    }
  }
  return(int)
}

y1matrix=matrix(0,Nt+1,Nx+1)
y2matrix=matrix(0,Nt+1,Nx+1)
for(n in 1:(Nx+1)){
  y1matrix[1,n]=integraly(0,ntox(n),function(t,x,s,xi) g(t,x,s,xi,k11,mu11,lambda11),function(t,x,s,xi){return(0)})
  y2matrix[1,n]=integraly(0,ntox(n),function(t,x,s,xi) g(t,x,s,xi,k21,mu21,lambda21),function(t,x,s,xi) g(t,x,s,xi,k22,mu22,lambda22))
}
