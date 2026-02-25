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
Lt=100
Lx=200
bufferx=15 #Should change depending on lambda
buffert=15 #Should change depending on mu 
xres=0.1 #should perhaps change with c??
tres=0.1
Nx=Lx/xres #must be integer
Nt=Lt/tres #must be integer

#Area to be integrated over
Lxi=Lx+2*bufferx
Ls=Lt+buffert
Nxi=Lxi/xres
Ns=Ls/tres

#Construction of the gaussian random field.
Wvec<-mvrnorm((Nxi+1)*(Ns+1),c(0,0),B*(xres*tres)) #check resolution multiplier correct
W1=matrix(Wvec[,1],Nxi+1,Ns+1)
W2=matrix(Wvec[,2],Nxi+1,Ns+1)



#
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

integralofg2<-function(t,x){
  if(t<=Lt&&x<=Lx/2&&x>=-Lx/2&&t>=0){
    output=0
    for(j in 0:(Ns/Ls*(t+buffert))){
      for(i in max(0,(x*Nxi/Lxi+Nxi/2-c*Nxi/Lxi*(t+buffert-j*Ls/Ns))):min((x*Nxi/Lxi+Nxi/2+c*Nxi/Lxi*(t+buffert-j*Ls/Ns)),Nxi)){ 
        output<-output+g(t,x,-buffert+j*Ls/Ns,-Lxi/2+i*Lxi/Nxi)^2*xres*tres
        #output<-output+1*xres*tres
        #print(c(-buffert+j*Ls/Ns,-Lxi/2+i*Lxi/Nxi))
      }
    }
    #output<-output+g(t,x,0,0)^2*xres*tres
    return(output)
  }
  else{
    print("invalid argument")
  }
}


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

outputfield<-ymatrixfun()
outputfield2<-ymatrixfun()
outputfield3<-ymatrixfun()

diagnostics<-function(ymatrix){
  heatmap(ymatrix,Rowv = NA,Colv=NA)
  return(var(as.vector(ymatrix)))
}
diagnostics(outputfield2)

heatmap(outputfield2,Rowv = NA,Colv=NA)

###ESTIMATION SECTION

#Estimation of mu11
mu11hat<-function(ymatrix,T){
  lentime<-dim(ymatrix)[1]
  lenspace<-dim(ymatrix)[2]
  C110<-rep(0,lenspace)
  C11T<-rep(0,lenspace)
  for(j in 1:lenspace){
    C110[j]=cov(ymatrix[1:lentime,j],ymatrix[1:lentime,j])
    C11T[j]=cov(ymatrix[1:(lentime-T),j],ymatrix[(1+T):lentime,j])
  }
  return(1/T*log(mean(C110)/mean(C11T)))
}
#Estimation of lambda11
lambda11hat<-function(ymatrix,L,mu11){
  lentime<-dim(ymatrix)[1]
  lenspace<-dim(ymatrix)[2]
  C110<-rep(0,lentime)
  C11L<-rep(0,lentime)
  for(i in 1:lentime){
    C110[i]=cov(ymatrix[i,1:lenspace],ymatrix[i,1:lenspace])
    C11L[i]=cov(ymatrix[i,1:(lenspace-L)],ymatrix[i,(1+L):lenspace])
  }
  print(mean(C11L))
  print(mean(C110))
  lambda11fun<-function(l){
    mean(C110)/(mean(C11L)*mu11)*exp(-2*mu11*L/c-l*L)*(mu11*exp(mu11*L/c)+c*l*(exp(mu11*L/c)-1))-1
    #1/(0.22*mu11)*exp(-2*mu11*L/c-l*L)*(mu11*exp(mu11*L/c)+c*l*(exp(mu11*L/c)-1))-1
  }
  #print(lambda11fun(0))
  print(lambda11fun(1))
  #print(lambda11fun(2))
  rootsolve<-uniroot(lambda11fun,c(0,3))
  return(rootsolve$root)
}
lambda11hat(outputfield2,1,1)

mu11hat(outputfield,1)
spatialcovs=rep(0,dim(outputfield)[1])
for(i in 1:dim(outputfield)[1]){
  spatialcovs[i]<-cov(outputfield[i,1:200],outputfield[i,2:201])
}

