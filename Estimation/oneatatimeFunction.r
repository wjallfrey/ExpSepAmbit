library(ggplot2)
library(pracma)
library(lamW)
library(ggpubr)
library(reshape2)

gammaijTimeHAT<-function(i,j,Tlag,Y,mean){
  if(missing(mean)){
    mean=F
  }
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-colMeans(0.5*(yi[(Tlag/(2*deltat)+1):Nt,1:Nx]-yi[1:(Nt-Tlag/(2*deltat)),1:Nx])*(yj[(Tlag/(2*deltat)+1):Nt,1:Nx]-yj[1:(Nt-Tlag/(2*deltat)),1:Nx]))
  if(mean){
    return(mean(gam))
  }
  else{return(gam)}
}
gammaijSpaceHAT<-function(i,j,Llag,Y,mean){
  if(missing(mean)){
    mean=F
  }
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-rowMeans(0.5*(yi[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yi[1:Nt,1:(Nx-Llag/(2*deltat*c))])*(yj[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yj[1:Nt,1:(Nx-Llag/(2*deltat*c))]))
  if(mean){
    return(mean(gam))
  }
  else{return(gam)}
}
Cijhat0<-function(i,j,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  return(mean(yi*yj))
}
CijhatTime<-function(i,j,Tlag,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  return(colMeans(yi[1:(Nt-Tlag/(2*deltat)),1:(Nx)]*yj[(Tlag/(2*deltat)+1):Nt,(1):Nx]))
}
CijhatSpace<-function(i,j,Llag,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  return(rowMeans(yi[1:(Nt),1:(Nx-Llag/(2*deltat*c))]*yj[(1):Nt,(1+Llag/(2*deltat*c)):Nx]))
}

C11Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  mu11<-parameters$mu11
  lambda11<-parameters$lambda11
  return(c*k11^2/(2*mu11*(mu11+c*lambda11))*(exp(-mu11*Tlag)))
}
C12Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  return((2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21)))*(exp(-mu21*Tlag)))
}
C21Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  return((2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21)))*(exp(-mu11*Tlag)))
}
C22Time<-function(Tlag,parameters,c){
  k21<-parameters$k21
  k22<-parameters$k22
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return(c*k21^2/(2*mu21*(mu21+c*lambda21))*(exp(-mu21*Tlag))+(c*k22^2/(2*mu22*(mu22+c*lambda22)))*(exp(-mu22*Tlag)))
}
C11Space<-function(Llag,parameters,c){
  k11<-parameters$k11
  mu11<-parameters$mu11
  lambda11<-parameters$lambda11
  return(c*k11^2/(2*mu11*(mu11+c*lambda11))*(exp(-(mu11/c+lambda11)*Llag)*(1+c*lambda11/mu11*(1-exp(-mu11/c*Llag)))))
}
C12Space<-function(Llag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  return(c*k11*k21*exp(-(mu11+mu21)*(Llag/c))/(mu11+mu21)*(
    (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*lambda11+c*lambda21)+
      exp(-lambda11*Llag)*(exp((mu11+mu21+c*lambda11-c*lambda21)*Llag/(2*c))-1)/(mu11+mu21+c*lambda11-c*lambda21)+
      exp(-lambda21*Llag)*(exp((mu11+mu21-c*lambda11+c*lambda21)*Llag/(2*c))-1)/(mu11+mu21-c*lambda11+c*lambda21)
  )
  )
}
C22Space<-function(Llag,parameters,c){
  k21<-parameters$k21
  k22<-parameters$k22
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return(c*k21^2/(2*mu21*(mu21+c*lambda21))*(exp(-(mu21/c+lambda21)*Llag)*(1+c*lambda21/mu21*(1-exp(-mu21/c*Llag))))+
           c*k22^2/(2*mu22*(mu22+c*lambda22))*(exp(-(mu22/c+lambda22)*Llag)*(1+c*lambda22/mu22*(1-exp(-mu22/c*Llag)))))
}

normalisedgammaijTimeHAT<-function(i,j,Tlag,Y,mean){
  if(missing(mean)){
    mean=F
  }
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-colMeans(0.5*(yi[(Tlag/(2*deltat)+1):Nt,1:Nx]-yi[1:(Nt-Tlag/(2*deltat)),1:Nx])*(yj[(Tlag/(2*deltat)+1):Nt,1:Nx]-yj[1:(Nt-Tlag/(2*deltat)),1:Nx]))
  cov0<-mean((yi[1:Nt,1:Nx])*(yj[1:Nt,1:Nx]))
  if(mean){
    return(mean(gam/cov0))
  }
  else{return(gam/cov0)}
}
normalisedgammaijSpaceHAT<-function(i,j,Llag,Y,mean){
  if(missing(mean)){
    mean=F
  }
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-rowMeans(0.5*(yi[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yi[1:Nt,1:(Nx-Llag/(2*deltat*c))])*(yj[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yj[1:Nt,1:(Nx-Llag/(2*deltat*c))]))
  cov0<-mean((yi[1:Nt,1:Nx])*(yj[1:Nt,1:Nx]))
  if(mean){
    return(mean(gam/cov0))
  }
  else{return(gam/cov0)}
}

normalisedtheoreticalgamma11Time<-function(Tlag,parameters,c){1-C11Time(Tlag,parameters,c)/C11Time(0,parameters,c)}
normalisedtheoreticalgamma12Time<-function(Tlag,parameters,c){1-0.5*(C12Time(Tlag,parameters,c)+C21Time(Tlag,parameters,c))/C12Time(0,parameters,c)}
normalisedtheoreticalgamma22Time<-function(Tlag,parameters,c){1-C22Time(Tlag,parameters,c)/C22Time(0,parameters,c)}
normalisedtheoreticalgamma11Space<-function(Llag,parameters,c){1-C11Space(Llag,parameters,c)/C11Space(0,parameters,c)}
normalisedtheoreticalgamma12Space<-function(Llag,parameters,c){1-C12Space(Llag,parameters,c)/C12Space(0,parameters,c)}
normalisedtheoreticalgamma22Space<-function(Llag,parameters,c){1-C22Space(Llag,parameters,c)/C22Space(0,parameters,c)}

aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}

gammahatsijTime<-function(i,j,Y,timelags){
  Nx<-dim(Y[[1]])[2]
  output<-matrix(0,NumberofLags,Nx)
  NumberofLags<-length(timelags)
  for(k in 1:NumberofLags){
    Tlag<-timelags[k]
    output[k,]<-normalisedgammaijTimeHAT(i,j,Tlag,Y)
  }
  return(output)
}
gammahatsijSpace<-function(i,j,Y,spacelags){
  Nt<-dim(Y[[1]])[1]
  output<-matrix(0,NumberofLags,Nt)
  NumberofLags<-length(spacelags)
  for(k in 1:NumberofLags){
    Llag<-spacelags[k]
    output[k,]<-normalisedgammaijSpaceHAT(i,j,Llag,Y)
  }
  return(output)
}
ChatsijTime<-function(i,j,Y,timelags){
  timelags0<-c(0,timelags)
  Nx<-dim(Y[[1]])[2]
  output<-matrix(0,length(timelags0),Nx)
  for(k in 1:length(timelags0)){
    Tlag<-timelags0[k]
    output[k,]<-CijhatTime(i,j,Tlag,Y)
  }
  return(output)
}

E11Time<-function(mu11,estimatedparas,Y,timelags,gammahats){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  E11timemat<-matrix(0,NumberofLags,Nx)
  parameters<-estimatedparas
  parameters$mu11<-mu11
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E11timemat[i,]<-gammahats[i,]-normalisedtheoreticalgamma11Time(Tlag,parameters,c)
  }
  return(E11timemat)
}
loss11time<-function(mu11,estimatedparas,Y,timelags,S,gammahats){
  e11<-E11Time(mu11,estimatedparas,Y,timelags,gammahats)
  lossfunction<-t(rowSums(e11))%*%S%*%(rowSums(e11))
  return(lossfunction)
} #Should be scaled maybe if combining with others but for now we are not

E12Time<-function(mu21,estimatedparas,Y,timelags,gammahats){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  E12timemat<-matrix(0,NumberofLags,Nx)
  parameters<-estimatedparas
  parameters$mu21<-mu21
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E12timemat[i,]<-gammahats[i,]-normalisedtheoreticalgamma12Time(Tlag,parameters,c)
  }
  return(E12timemat)
}
loss12time<-function(mu21,estimatedparas,Y,timelags,S,gammahats){
  e12<-E12Time(mu21,estimatedparas,Y,timelags,gammahats)
  lossfunction<-t(rowSums(e12))%*%S%*%(rowSums(e12))
  return(lossfunction)
}

E11Space<-function(lambda11,estimatedparas,Y,spacelags,gammahats){
  NumberofLags<-length(spacelags)
  c<-Y[[5]]
  Nt<-dim(Y[[1]])[1]
  E11spaceemat<-matrix(0,NumberofLags,Nt)
  parameters<-estimatedparas
  parameters$lambda11<-lambda11
  for(i in 1:NumberofLags){
    Llag<-spacelags[i]
    E11spaceemat[i,]<-gammahats[i,]-normalisedtheoreticalgamma11Space(Llag,parameters,c)
  }
  return(E11spaceemat)
}
loss11space<-function(lambda11,estimatedparas,Y,spacelags,S,gammahats){
  e11<-E11Space(lambda11,estimatedparas,Y,spacelags,gammahats)
  lossfunction<-t(rowSums(e11))%*%S%*%(rowSums(e11))
  return(lossfunction)
}

E12Space<-function(lambda21,estimatedparas,Y,spacelags,gammahats){
  NumberofLags<-length(spacelags)
  c<-Y[[5]]
  Nt<-dim(Y[[1]])[1]
  E12spaceemat<-matrix(0,NumberofLags,Nt)
  parameters<-estimatedparas
  parameters$lambda21<-lambda21
  for(i in 1:NumberofLags){
    Llag<-spacelags[i]
    E12spaceemat[i,]<-gammahats[i,]-normalisedtheoreticalgamma12Space(Llag,parameters,c)
  }
  return(E12spaceemat)
}
loss12space<-function(lambda21,estimatedparas,Y,spacelags,S,gammahats){
  e12<-E12Space(lambda21,estimatedparas,Y,spacelags,gammahats)
  lossfunction<-t(rowSums(e12))%*%S%*%(rowSums(e12))
  return(lossfunction)
}


CovE11Time<-function(k11,estimatedparas,Y,timelags,Chats){
  timelags0<-c(0,timelags) #including the zero lag as we now have covaraince
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  CovE11timeemat<-matrix(0,length(timelags0),Nx)
  parameters<-estimatedparas
  parameters$k11<-k11
  for(i in 1:length(timelags0)){
    Tlag<-timelags0[i]
    CovE11timeemat[i,]<-Chats[i,]-C11Time(timelags0[i],parameters,c)
  }
  return(CovE11timeemat)
}
lossCov11time<-function(k11,estimatedparas,Y,timelags,S,Chats){
  e11<-CovE11Time(k11,estimatedparas,Y,timelags,Chats)
  lossfunction<-t(rowSums(e11))%*%S%*%rowSums(e11)
  return(lossfunction)
}

CovE12Time<-function(k21,estimatedparas,Y,timelags,Chats){
  timelags0<-c(0,timelags) #including the zero lag as we now have covaraince
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  CovE12timeemat<-matrix(0,length(timelags0),Nx)
  parameters<-estimatedparas
  parameters$k21<-k21
  for(i in 1:length(timelags0)){
    Tlag<-timelags0[i]
    CovE12timeemat[i,]<-Chats[i,]-C12Time(timelags0[i],parameters,c)
  }
  return(CovE12timeemat)
}
lossCov12time<-function(k21,estimatedparas,Y,timelags,S,Chats){
  e12<-CovE12Time(k21,estimatedparas,Y,timelags,Chats)
  lossfunction<-t(rowSums(e12))%*%S%*%rowSums(e12)
  return(lossfunction)
}

half22normalisedgamma22Time<-function(Tlag,parameters,c){
  parametersk22is0<-parameters
  parametersk22is0$k22<-0
  1-(C22Time(Tlag,parameters,c)-C22Time(Tlag,parametersk22is0,c))/(C22Time(0,parameters,c)-C22Time(0,parametersk22is0,c))
}
half22normalisedgamma22TimeHAT<-function(Tlag,parameters,c,Y){
  parametersk22is0<-parameters
  parametersk22is0$k22<-0
  1-(CijhatTime(2,2,Tlag,Y)-C22Time(Tlag,parametersk22is0,c))/(CijhatTime(2,2,0,Y)-C22Time(0,parametersk22is0,c))
}
E22timeHalf<-function(mu22,estimatedparas,Y,timelags){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  E22timemat<-matrix(0,NumberofLags,Nx)
  parameters<-estimatedparas
  parameters$mu22<-mu22
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E22timemat[i,]<-half22normalisedgamma22TimeHAT(Tlag,parameters,c,Y)-half22normalisedgamma22Time(Tlag,parameters,c)
  }
  return(E22timemat)
}
loss22timeHalf<-function(mu22,estimatedparas,Y,timelags,S){
  e22<-E22timeHalf(mu22,estimatedparas,Y,timelags)
  lossfunction<-t(rowSums(e22))%*%S%*%(rowSums(e22))
  return(lossfunction)
}
half22normalisedgamma22Space<-function(Llag,parameters,c){
  parametersk22is0<-parameters
  parametersk22is0$k22<-0
  1-(C22Space(Llag,parameters,c)-C22Space(Llag,parametersk22is0,c))/(C22Space(0,parameters,c)-C22Space(0,parametersk22is0,c))
}
half22normalisedgamma22SpaceHAT<-function(Llag,parameters,c,Y){
  parametersk22is0<-parameters
  parametersk22is0$k22<-0
  1-(CijhatSpace(2,2,Llag,Y)-C22Space(Llag,parametersk22is0,c))/(CijhatSpace(2,2,0,Y)-C22Space(0,parametersk22is0,c))
}
E22spaceHalf<-function(lambda22,estimatedparas,Y,spacelags){
  NumberofLags<-length(spacelags)
  c<-Y[[5]]
  Nt<-dim(Y[[1]])[1]
  E22spacemat<-matrix(0,NumberofLags,Nt)
  parameters<-estimatedparas
  parameters$lambda22<-lambda22
  for(i in 1:NumberofLags){
    Llag<-spacelags[i]
    E22spacemat[i,]<-half22normalisedgamma22SpaceHAT(Llag,parameters,c,Y)-half22normalisedgamma22Space(Llag,parameters,c)
  }
  return(E22spacemat)
}
loss22spaceHalf<-function(lambda22,estimatedparas,Y,spacelags,S){
  e22<-E22spaceHalf(lambda22,estimatedparas,Y,spacelags)
  lossfunction<-t(rowSums(e22))%*%S%*%(rowSums(e22))
  return(lossfunction)
}

CovE22Time<-function(k22,estimatedparas,Y,timelags){
  timelags0<-c(0,timelags) #including the zero lag as we now have covaraince
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  CovE22timeemat<-matrix(0,length(timelags0),Nx)
  parameters<-estimatedparas
  parameters$k22<-k22
  for(i in 1:length(timelags0)){
    Tlag<-timelags0[i]
    CovE22timeemat[i,]<-CijhatTime(2,2,timelags0[i],Y)-C22Time(timelags0[i],parameters,c)
  }
  return(CovE22timeemat)
}
lossCov22time<-function(k22,estimatedparas,Y,timelags,S){
  e22<-CovE22Time(k22,estimatedparas,Y,timelags)
  lossfunction<-t(rowSums(e22))%*%S%*%rowSums(e22)
  return(lossfunction)
}


getmu11<-function(Y,initialmu11,estimatedparas,timelags,tol,Weights){
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(timelags)
  S<-diag(NumberofLags)
  dif<-1
  gammahats<-gammahatsijTime(1,1,Y,timelags)
  while(dif>tol){
    mu11hat<-optim(initialmu11,function(mu11){loss11time(mu11,estimatedparas,Y,timelags,S,gammahats)},method="BFGS")$par
    if(!Weights){break}
    e11<-E11Time(mu11hat,estimatedparas,Y,timelags,gammahats)
    omega<-e11%*%t(e11)
    S<-solve(omega)
    dif<-abs(initialmu11-mu11hat)/mu11hat
    initialmu11<-mu11hat
    #print(c(mu11hat,dif))
  }
  return(mu11hat)
}
estimatedparas$mu11<-getmu11(Y,1,estimatedparas,timelags,0.01,T)
getmu21<-function(Y,initialmu21,estimatedparas,timelags,tol,Weights){
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(timelags)
  S<-diag(NumberofLags)
  dif<-1
  gammahats<-gammahatsijTime(1,2,Y,timelags)
  while(dif>tol){
    mu21hat<-optim(initialmu21,function(mu21){loss12time(mu21,estimatedparas,Y,timelags,S,gammahats)},method="BFGS")$par
    if(!Weights){break}
    e12<-E12Time(mu21hat,estimatedparas,Y,timelags,gammahats)
    omega<-e12%*%t(e12)
    S<-solve(omega)
    dif<-abs(initialmu21-mu21hat)/mu21hat
    initialmu21<-mu21hat
  }
  return(mu21hat)
}
estimatedparas$mu21<-getmu21(Y,1,estimatedparas,timelags,0.01,T)
getlambda11<-function(Y,initiallambda11,estimatedparas,spacelags,tol,Weights){
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(spacelags)
  S<-diag(NumberofLags)
  dif<-1
  gammahats<-gammahatsijSpace(1,1,Y,spacelags)
  while(dif>tol){
    lambda11hat<-optim(initiallambda11,function(lambda11){loss11space(lambda11,estimatedparas,Y,spacelags,S,gammahats)},method="BFGS")$par
    if(!Weights){break}
    e11<-E11Space(lambda11hat,estimatedparas,Y,spacelags,gammahats)
    omega<-e11%*%t(e11)
    S<-solve(omega)
    dif<-abs(initiallambda11-lambda11hat)/lambda11hat
    initiallambda11<-lambda11hat

  }
  return(lambda11hat)
}
estimatedparas$lambda11<-getlambda11(Y,1,estimatedparas,timelags,0.01,T)
getlambda21<-function(Y,initiallambda21,estimatedparas,spacelags,tol,Weights){
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(spacelags)
  S<-diag(NumberofLags)
  dif<-1
  gammahats<-gammahatsijSpace(1,2,Y,spacelags)
  while(dif>tol){
    lambda21hat<-optim(initiallambda21,function(lambda21){loss12space(lambda21,estimatedparas,Y,spacelags,S,gammahats)},method="BFGS")$par
    if(!Weights){break}
    e12<-E12Space(lambda21hat,estimatedparas,Y,spacelags,gammahats)
    omega<-e12%*%t(e12)
    S<-solve(omega)
    dif<-abs(initiallambda21-lambda21hat)/lambda21hat
    initiallambda21<-lambda21hat
  }
  return(lambda21hat)
}
estimatedparas$lambda21<-getlambda21(Y,1,estimatedparas,timelags,0.01,T)
getk11<-function(Y,initialk11,estimatedparas,timelags,tol,Weights){
  Chats<-ChatsijTime(1,1,Y,timelags)
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(timelags)
  S<-diag(NumberofLags+1)
  dif<-1
  while(dif>tol){
    k11hat<-optim(initialk11,function(k11){lossCov11time(k11,estimatedparas,Y,timelags,S,Chats)},method="BFGS")$par
    if(!Weights){break}
    e11<-CovE11Time(k11hat,estimatedparas,Y,timelags,Chats)
    omega<-e11%*%t(e11)
    S<-solve(omega)
    dif<-abs(initialk11-k11hat)/k11hat
    initialk11<-k11hat
  }
  return(k11hat)
}
estimatedparas$k11<-getk11(Y,1,estimatedparas,timelags,0.01,T)
getk21<-function(Y,initialk21,estimatedparas,timelags,tol,Weights){
  Chats<-ChatsijTime(1,2,Y,timelags)
  if(missing(Weights)){Weights<-T}
  NumberofLags<-length(timelags)
  S<-diag(NumberofLags+1)
  dif<-1
  while(dif>tol){
    k21hat<-optim(initialk21,function(k21){lossCov12time(k21,estimatedparas,Y,timelags,S,Chats)},method="BFGS")$par
    if(!Weights){break}
    e12<-CovE12Time(k21hat,estimatedparas,Y,timelags,Chats)
    omega<-e12%*%t(e12)
    S<-solve(omega)
    dif<-abs(initialk21-k21hat)/k21hat
    initialk21<-k21hat
  }
  return(k21hat)
}
estimatedparas$k21<-getk21(Y,1,estimatedparas,timelags,0.01,T)




getlambda21(Y,1,estimatedparas,spacelags,0.01,T)


initialmu21<-1
S<-diag(NumberofLags)
for(i in 1:5){
  mu21hat<-optim(initialmu21,function(mu21){loss12time(mu21,estimatedparas,Y,timelags,S)})$par
  e12<-E12Time(mu21hat,estimatedparas,Y,timelags)
  omega<-e12%*%t(e12)
  S<-solve(omega)
  initialmu21<-mu21hat
  print(mu21hat)
}
estimatedparas$mu21<-mu21hat


initiallambda11<-1
S<-diag(NumberofLags)
for(i in 1:5){
  lambda11hat<-optim(initiallambda11,function(lambda11){loss11space(lambda11,estimatedparas,Y,spacelags,S)})$par
  e11<-E11Space(lambda11hat,estimatedparas,Y,spacelags)
  omega<-e11%*%t(e11)
  S<-solve(omega)
  initiallambda11<-lambda11hat
  print(lambda11hat)
}
estimatedparas$lambda11<-lambda11hat


initiallambda21<-1
S<-diag(NumberofLags)
for(i in 1:5){
  lambda21hat<-optim(initiallambda21,function(lambda21){loss12space(lambda21,estimatedparas,Y,spacelags,S)})$par
  e12<-E12Space(lambda21hat,estimatedparas,Y,spacelags)
  omega<-e12%*%t(e12)
  S<-solve(omega)
  initiallambda21<-lambda21hat
  print(lambda21hat)
}
estimatedparas$lambda21<-lambda21hat


initialk11<-1
S<-diag(NumberofLags+1)
for(i in 1:5){
  k11hat<-optim(initialk11,function(k11){lossCov11time(k11,estimatedparas,Y,timelags,S)})$par
  e11<-CovE11Time(k11hat,estimatedparas,Y,timelags)
  omega<-e11%*%t(e11)
  S<-solve(omega)
  initialk11<-k11hat
  print(k11hat)
}
estimatedparas$k11<-k11hat


initialk21<-1
S<-diag(NumberofLags+1)
for(i in 1:5){
  k21hat<-optim(initialk21,function(k21){lossCov12time(k21,estimatedparas,Y,timelags,S)})$par
  e12<-CovE12Time(k21hat,estimatedparas,Y,timelags)
  omega<-e12%*%t(e12)
  S<-solve(omega)
  initialk21<-k21hat
  print(k21hat)
}
estimatedparas$k21<-k21hat


##Weird variogram alternate for 22. Time

initialmu22<-1
NumberofLags<-length(timelags)
S<-diag(NumberofLags)
for(i in 1:5){
  mu22hat<-optim(initialmu22,function(mu22){loss22timeHalf(mu22,estimatedparas,Y,timelags,S)})$par
  e22<-E22timeHalf(mu22hat,estimatedparas,Y,timelags)
  omega<-e22%*%t(e22)
  S<-solve(omega)
  initialmu22<-mu22hat
  print(mu22hat)
}
estimatedparas$mu22<-mu22hat


#And space

initiallambda22<-1
NumberofLags<-length(timelags)
S<-diag(NumberofLags)
for(i in 1:5){
  lambda22hat<-optim(initiallambda22,function(lambda22){loss22spaceHalf(lambda22,estimatedparas,Y,timelags,S)})$par
  e22<-E22spaceHalf(lambda22hat,estimatedparas,Y,spacelags)
  omega<-e22%*%t(e22)
  S<-solve(omega)
  initialmu22<-lambda22hat
  print(lambda22hat)
}
estimatedparas$lambda22<-lambda22hat

#For k22

initialk22<-1
S<-diag(NumberofLags+1)
for(i in 1:5){
  k22hat<-optim(initialk22,function(k22){lossCov22time(k22,estimatedparas,Y,timelags,S)})$par
  e22<-CovE22Time(k22hat,estimatedparas,Y,timelags)
  omega<-e22%*%t(e22)
  S<-solve(omega)
  initialk22<-k22hat
  print(k22hat)
}
estimatedparas$k22<-k22hat


estimatedparas<-aslist(c(1,1,1,1,1,1,1,1,1))
Y<-readRDS("Output Fields/newparasfine3.rds")
c<-Y[[5]]
deltat<-Y[[6]]
NumberofLags<-15
timelags<-((1:NumberofLags))*2*deltat
spacelags<-((1:NumberofLags))*2*deltat*c
