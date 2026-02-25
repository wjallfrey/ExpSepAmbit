setwd("/Users/willallfrey/Documents/R")

k11=1
lambda11=1
mu11=1
c=1
B=matrix(c(1,0,0,1),2,2)

g<-function(t,x,s,xi){
  return(k11*exp(-mu11*abs(t-s)-lambda11*abs(x-xi)))
}

#Area where we want to find Y(t,x):
Lt=100 #t in [0,Lt]
Lx=200 #x in [-Lx/2,Lx/2]
bufferx=15 #Should change depending on lambda
buffert=15 #Should change depending on mu 
xres=1 #should perhaps change with c??
tres=1
Nx=Lx/xres #must be integer
Nt=Lt/tres #must be integer

#Area to be integrated over: Need c*Nxi*Ls=Lxi*Ns
Lxi=Lx+2*bufferx
Ls=Lt+buffert
xires=1
sres=1
Nxi=Lxi/xires
Ns=Ls/sres

#Construction of the gaussian random field
Wvec<-mvrnorm((Nxi+1)*(Ns+1),c(0,0),B*(xires*sres)) #check resolution multiplier correct
W1=matrix(Wvec[,1],Nxi+1,Ns+1)
W2=matrix(Wvec[,2],Nxi+1,Ns+1)

jtos<-function(j){
  if(j>=1&&j<=(Ns+1)){
    return(-buffert+(j-1)*Ls/Ns)
  }
  else{
    print("j out of range")
  }
}
stoj<-function(s){
  if(s<=Lt&&s>=-buffert){
    return(Ns/Ls*(s+buffert)+1)
  }
  else{
    print("s out of range")
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

xitoi<-function(xi){
  if(-Lxi/2<=xi&&xi<=Lxi/2){
    return(Nxi/Lxi*(xi+Lxi/2)+1)
  }
  else{
    print("xi out of range")
  }
}

integraly<-function(t,x){
  int<-0
  for(j in 1:(Ns/Ls*(t+buffert)+1)){
    for(i in (max(1,Nxi/Lxi*(x+Lxi/2-c(t-jtos(j)))+1)):(min(Nxi+1,1+Nxi/Lxi*(x+Lxi/2+c(t-jtos(j)))))){
      int<-int+g(t,x,jtos(j),itoxi(i))*W1[i,j]
    }
  }
  return(int)
}


deltat=Ls/Ns
ytplusfromyt<-function(t,x,ytx){
  if(t+deltat<=Lt){
    output=exp(-mu11*deltat)*ytx #Unique to this exponential form of g
    for(j in 1:(Ns/Ls*(t+buffert)+1)){
      iplus<-Nxi/Lxi*(x+Lxi/2+c*(t+deltat-jtos(j)))+1
      iminus<-Nxi/Lxi*(x+Lxi/2-c*(t+deltat-jtos(j)))+1
      if(iminus>=1){
        output<-output+g(t+deltat,x,jtos(j),itoxi(iminus))*W1[iminus,j]
      }
      if(iplus<=(Nxi+1)){
        output<-output+g(t+deltat,x,jtos(j),itoxi(iplus))*W1[iplus,j]
      }
    }
    apexj<-(Ns/Ls*(t+deltat+buffert)+1)
    apexi<-Nxi/Lxi*(x+Lxi/2)+1
    output<-output+g(t+deltat,x,jtos(apexj),itoxi(apexi))*W1[apexi,apexj]
    return(output)
  }
  else{
    print("increment goes out of range")
  }
}

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

ytplusfromyt(0,0,integraly(0,0))
integraly(1,0)



#Make the output matrix
ymatrix=matrix(0,Nt+1,Nx+1)

#use full integral method for first row
for(n in 1:(Nx+1)){
  ymatrix[1,n]=integraly(0,ntox(n))
}

#use iterative method for the rest
for(m in 2:(Nt+1)){
  for(n in 1:(Nx+1)){
    ymatrix[m,n]=ytplusfromyt(mtot(m-1),ntox(n),ymatrix[m-1,n])
  }
}

heatmap(ymatrix,Rowv = NA,Colv=NA)
