library(stats)
library(ggplot2)
library(reshape2)
library(lattice)
library(gridExtra)
library(gganimate)
library(MASS)
library(parallel)

library(lattice,MASS,parallel,ggplot2,reshape2,gridExtra,gganimate,pracma,lamW,ggpubr)

paraset2<-list(k11=0.5,k21=0.3,k22=0.4,mu11=0.2,mu21=0.8,mu22=0.4,lambda11=0.3,lambda21=0.9,lambda22=0.5)
parameters<-paraset2
deltat=0.5
c=1
deltax=c*deltat
tol=0.0001
Lx<-200
Lt<-100
B<-diag(2)
set.seed(123)

k11<-parameters$k11
k21<-parameters$k21
k22<-parameters$k22
mu11<-parameters$mu11
mu21<-parameters$mu21
mu22<-parameters$mu22
lambda11<-parameters$lambda11
lambda21<-parameters$lambda21
lambda22<-parameters$lambda22

g<-function(t,x,s,xi,k,mu,lambda){
  return(k*exp(-mu*abs(t-s)-lambda*abs(x-xi)))
}


truncationerror<-function(alpha,mu,lambda,c){
  if(c*lambda==mu){
    return(exp(-alpha*deltat*mu)*(1+2*mu*(alpha*deltat/2)))
  }
  else{
    return(1-exp(-alpha*deltat*mu)*(1+2*mu*(1-exp(-alpha*deltat/2*(c*lambda-mu)))/(c*lambda-mu)))
  }
}
alpha<-ceiling(uniroot(function(alphavar){truncationerror(alphavar,2*min(mu11,mu21,mu22),2*min(lambda11,lambda21,lambda22),c)-1+tol},c(0,100))$root)

xisfieldsize<-Lx/deltax+(Lt)/deltat+alpha

Wvec<-mvrnorm(xisfieldsize^2,c(0,0),B*(0.5*deltat*deltax))
W1<-matrix(Wvec[,1],xisfieldsize,xisfieldsize)
W2<-matrix(Wvec[,2],xisfieldsize,xisfieldsize)

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
    2/(deltat^2*(mu+c*lambda)*(mu-c*lambda))*(cosh(mu*deltat)-cosh(c*lambda*deltat))
  }
  else{
    2/(deltat^2*(mu+c*lambda))*(deltat*sinh(deltat*mu))
  }
}
kappaKD<-function(mu,lambda){
  if(mu!=c*lambda){
    2/(deltat^2*(mu+c*lambda)*(mu-c*lambda))*(-c*lambda/mu*sinh(mu*deltat)+sinh(c*lambda*deltat))
  }
  else{
    2/(deltat^2*(mu+c*lambda))*(-1*deltat*cosh(deltat*mu)+sinh(deltat*mu)/mu)
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
vardirect<-function(k,mu,lambda,alpha){
  variance<-0
  variancenokappa<-0
  iapex<-401
  japex<-401
  for(i in iapex:(iapex+alpha-1)){
    for(j in japex:(japex+alpha-1)){
      variance<-variance+g(0,0,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)^2*(kappa(mu,lambda)+((i-iapex)==(j-japex))*kappaKD(mu,lambda))*(0.5*deltat*deltax)
      variancenokappa<-variancenokappa+g(0,0,ijtos(i,j),ijtoxi(i,j),k,mu,lambda)^2*(0.5*deltat*deltax)
    }
  }
  return(c(variance,variancenokappa))
}
vardirect(k11,mu11,lambda11,1)
g(0,0,ijtos(401,401),ijtoxi(401,401),k11,mu11,lambda11)^2*(0.5*deltat*deltax)*(kappa(mu11,lambda11)+kappaKD(mu11,lambda11))

#ydtfromytx<-function(t,x,ytx,k,mu,lambda,W){
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

#for(n in 1:(Lx/(deltax)+1)){
#  y11matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k11,mu11,lambda11,W1)
#  y21matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k21,mu21,lambda21,W1)
#  y22matrix[1,n]<-integraldirect(0,-Lx/2+(n-1)*deltax,k22,mu22,lambda22,W2)
#}
#for(m in 2:(Lt/(deltat)+1)){
#  for(n in 1:(Lx/(deltax)+1)){
#    y11matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y11matrix[m-1,n],k11,mu11,lambda11,W1)
#    y21matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y21matrix[m-1,n],k21,mu21,lambda21,W1)
#    y22matrix[m,n]<-ydtfromytx((m-1)*deltat,-Lx/2+(n-1)*deltax,y22matrix[m-1,n],k22,mu22,lambda22,W2)
#  }
#}
for(m in 1:(Lt/(deltat)+1)){
  for(n in 1:(Lx/(deltax)+1)){
    y11matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k11,mu11,lambda11,W1)
    y21matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k21,mu21,lambda21,W1)
    y22matrix[m,n]<-integraldirect((m-1)*deltat,-Lx/2+(n-1)*deltax,k22,mu22,lambda22,W2)
  }
}

y1matrix=y11matrix+y12matrix
y2matrix=y21matrix+y22matrix
