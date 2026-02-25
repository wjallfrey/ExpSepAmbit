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

#Area to be integrated over
Lxi=Lx+2*bufferx
Ls=Lt+buffert
Nxi=Lxi/xres
Ns=Ls/tres

#Construction of the gaussian random field
Wvec<-mvrnorm((Nxi+1)*(Ns+1),c(0,0),B*(xres*tres)) #check resolution multiplier correct
W1=matrix(Wvec[,1],Nxi+1,Ns+1)
W2=matrix(Wvec[,2],Nxi+1,Ns+1)


#Make the output matrix
ymatrix=matrix(0,Nt+1,Nx+1)

#define the inefficient y(t,x) function
y<-function(t,x){
  if(t<=Lt&&x<=Lx/2&&x>=-Lx/2&&t>=0){
    output=0
    for(j in 0:(Ns/Ls*(t+buffert))){
      for(i in max(0,(x*Nxi/Lxi+Nxi/2-c*Nxi/Lxi*(t+buffert-j*Ls/Ns))):min((x*Nxi/Lxi+Nxi/2+c*Nxi/Lxi*(t+buffert-j*Ls/Ns)),Nxi)){ 
        output<-output+g(t,x,-buffert+j*Ls/Ns,-Lxi/2+i*Lxi/Nxi)*W1[i+1,j+1]
      }
    }
    return(output)
  }
  else{
    print("invalid argument")
  }
}

#Slow method for comaparison
ymatrixfun<-function(){
  outputmatrix=matrix(0,Nt+1,Nx+1)
  for(i in 1:(Nt+1)){
    for(j in 1:(Nx+1)){
      tvalue<-(i-1)*Lt/Nt
      xvalue<-(-Lx/2+(j-1)*Lx/Nx)
      outputmatrix[i,j]<-y(tvalue,xvalue)
    }
  }
  return(outputmatrix)
}
ymatrixlong<-ymatrixfun()

#Use inefficient funciton for the t=0 row of the output matrix
for(j in 1:(Nx+1)){
  xvalue<-(-Lx/2+(j-1)*Lx/Nx)
  ymatrix[1,j]=y(0,xvalue)
}

deltat=1
output=exp(-mu11*deltat)*as.numeric(ymatrix[1,101])
for(j in 0:(Ns/Ls*(0+buffert))){
  iminus=0*Nxi/Lxi+Nxi/2-c*Nxi/Lxi*(0+deltat+buffert-j*Ls/Ns)
  iplus=0*Nxi/Lxi+Nxi/2+c*Nxi/Lxi*(0+deltat+buffert-j*Ls/Ns)
  if(iminus>=0){
    output=output+g(0,0,-buffert+j*Ls/Ns,0-c*(0+deltat+buffert-j*Ls/Ns))*W1[iminus+1,j+1]
  }
  if(iplus<=Nxi){
    output=output+g(0,0,-buffert+j*Ls/Ns,0+c*(0+deltat+buffert-j*Ls/Ns))*W1[iplus+1,j+1]
  }
  
}
j=(Ns/Ls*(t+buffert+deltat))
i=0*Nxi/Lxi+Nxi/2
output=output+g(0,0,-buffert+j*Ls/Ns,0+c*(0+deltat+buffert-j*Ls/Ns))*W1[i+1,j+1]


