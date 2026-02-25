library(MASS)

Lxi=400
Ls=200
Nxi=400
Ns=200

k11=1
lambda11=1
mu11=1
c=1


xivec=-Lxi/2+Lxi/Nxi*(0:Nxi)
svec=-Ls/2+Ls/Ns*(0:Ns)

B=Lxi*Ls/(Nxi*Ns)*matrix(c(1,0,0,1),2,2)

Wvec<-mvrnorm((Nxi+1)*(Ns+1),c(0,0),B)
W1=Wvec[,1]
W2=Wvec[,2]
W=matrix(rnorm((Nxi+1)*(Ns+1),0,Lxi*Ls/(Nxi*Ns)),Nxi+1,Ns+1)

g<-function(t,x,s,xi){
  return(k11*exp(-mu11*abs(t-s)-lambda11*abs(x-xi)))
}
sxiInAmbit<-function(t,x,s,xi){
  if(s<=t&&abs(x-xi)<=c*(t-s)){
    return(1)
  }
  else{
    return(0)
  }
}
y<-function(t,x){
  output=0
  for(i in 1:(Nxi+1)){
    for(j in 1:(Ns+1)){
      xi<-xivec[i]
      s<-svec[j]
      output=output+g(t,x,s,xi)*sxiInAmbit(t,x,s,xi)*W[i,j]
    }
  }
  return(output)
}
#we want to work out y at the grid points in say 0<t<50, -10<x<10, take a while

tvalues<-0:50
xvalues<-(-40):40
ymatrix<-matrix(0,51,81)
for(i in 1:51){
  for(j in 1:81){
    ymatrix[i,j]=y(tvalues[i],xvalues[j])
  }
}

#We estimate C_11(T)/C_11(0) at each x value. #Theory suggest should be indep of x, so we average over all values of x
mu11hat<-function(ymatrix,T){
  C110<-rep(0,81)
  C11T<-rep(0,81)
  for(j in 1:81){
    C110[j]=cov(ymatrix[1:51,j],ymatrix[1:51,j])
    C11T[j]=cov(ymatrix[1:(51-T),j],ymatrix[(1+T):51,j])
  }
  return(1/T*log(mean(C110)/mean(C11T)))
}
#lambda11hat<-function(ymatrix,L,mu11){} Can be defined in terms of Lambert W function

heatmap(ymatrix,Rowv = NA,Colv=NA)
mu11hat(ymatrix,1)

#spatial
cov(ymatrix[3,1:21],ymatrix[3,1:21])
cov(ymatrix[1,1:20],ymatrix[1,2:21])
