gammaijTime<-function(i,j,Tlag,Y){
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-colMeans(0.5*(yi[(Tlag/(2*deltat)+1):Nt,1:Nx]-yi[1:(Nt-Tlag/(2*deltat)),1:Nx])*(yj[(Tlag/(2*deltat)+1):Nt,1:Nx]-yj[1:(Nt-Tlag/(2*deltat)),1:Nx]))
  return(gam)
}
gammaijSpace<-function(i,j,Llag,Y){
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-rowMeans(0.5*(yi[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yi[1:Nt,1:(Nx-Llag/(2*deltat*c))])*(yj[1:Nt,(1+Llag/(2*deltat*c)):Nx]-yj[1:Nt,1:(Nx-Llag/(2*deltat*c))]))
  return(gam)
}

theoreticalgamma11Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  mu11<-parameters$mu11
  lambda11<-parameters$lambda11
  return(c*k11^2/(2*mu11*(mu11+c*lambda11))*(1-exp(-mu11*Tlag)))
}
theoreticalgamma12Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return((2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21)))*(1-0.5*(exp(-mu11*Tlag)+exp(-mu21*Tlag))))
}
theoreticalgamma22Time<-function(Tlag,parameters,c){
  k21<-parameters$k21
  k22<-parameters$k22
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return(c*k21^2/(2*mu21*(mu21+c*lambda21))*(1-exp(-mu21*Tlag))+(c*k22^2/(2*mu22*(mu22+c*lambda22)))*(1-exp(-mu22*Tlag)))
}
theoreticalgamma11Space<-function(Llag,parameters,c){
  k11<-parameters$k11
  mu11<-parameters$mu11
  lambda11<-parameters$lambda11
  return(c*k11^2/(2*mu11*(mu11+c*lambda11))*(1-exp(-(mu11/c+lambda11)*Llag)*(1+c*lambda11/mu11*(1-exp(-mu11/c*Llag)))))
}
theoreticalgamma12Space<-function(Llag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return(2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))
         -c*k11*k21*exp(-(mu11+mu21)*(Llag/c))/(mu11+mu21)*(
           (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*lambda11+c*lambda21)+
             exp(-lambda11*Llag)*(exp((mu11+mu21+c*lambda11-c*lambda21)*Llag/(2*c))-1)/(mu11+mu21+c*lambda11-c*lambda21)+
             exp(-lambda21*Llag)*(exp((mu11+mu21-c*lambda11+c*lambda21)*Llag/(2*c))-1)/(mu11+mu21-c*lambda11+c*lambda21)
         )
  )
}
theoreticalgamma22Space<-function(Llag,parameters,c){
  k21<-parameters$k21
  k22<-parameters$k22
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return(c*k21^2/(2*mu21*(mu21+c*lambda21))*(1-exp(-(mu21/c+lambda21)*Llag)*(1+c*lambda21/mu21*(1-exp(-mu21/c*Llag))))+
           c*k22^2/(2*mu22*(mu22+c*lambda22))*(1-exp(-(mu22/c+lambda22)*Llag)*(1+c*lambda22/mu22*(1-exp(-mu22/c*Llag)))))
}

theoreticalvariograms<-function(variogram,parameters,lags,c){
  NumberofLags<-length(lags)
  output<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    output[i]<-variogram(lags[i],parameters,c)
  }
  return(output)
}

datavariogram<-function(timelags,spacelags,Y){
  if(length(timelags)!=length(spacelags)){
    print("Need same number of space and time lags")
  }
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  NumberofLags<-length(timelags)
  g11hatt<-matrix(0,NumberofLags,Nx)
  g12hatt<-matrix(0,NumberofLags,Nx)
  g22hatt<-matrix(0,NumberofLags,Nx)
  g11hats<-matrix(0,NumberofLags,Nt)
  g12hats<-matrix(0,NumberofLags,Nt)
  g22hats<-matrix(0,NumberofLags,Nt)
  for(i in 1:NumberofLags){
    g11hatt[i,]<-gammaijTime(1,1,timelags[i],Y)
    g12hatt[i,]<-gammaijTime(1,2,timelags[i],Y)
    g22hatt[i,]<-gammaijTime(2,2,timelags[i],Y)
    g11hats[i,]<-gammaijSpace(1,1,spacelags[i],Y)
    g12hats[i,]<-gammaijSpace(1,2,spacelags[i],Y)
    g22hats[i,]<-gammaijSpace(2,2,spacelags[i],Y)
  }
  return(list(g11hatt=g11hatt,g12hatt=g12hatt,g22hatt=g22hatt,g11hats=g11hats,g12hats=g12hats,g22hats=g22hats))
}
aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}

E<-function(theta,gammahat,variogram,c,lags){
  (gammahat-theoreticalvariograms(variogram,aslist(theta),lags,c))/colMeans(gammahat)
}
onelossfunction<-function(E,S){
  t(rowSums(E))%*%S%*%rowSums(E)
}
getnewS<-function(E){
  return(solve(E%*%t(E)))
}
totallossfunction<-function(theta,datavars,Ss,c,timelags,spacelags){
  S11<-Ss$S11
  S12<-Ss$S12
  S22<-Ss$S22
  S11space<-Ss$S11space
  S12space<-Ss$S12space
  S22space<-Ss$S22space
  lossvector<-c(onelossfunction(E(theta,datavars$g11hatt,theoreticalgamma11Time,c,timelags),S11),
                onelossfunction(E(theta,datavars$g12hatt,theoreticalgamma12Time,c,timelags),S12),
                onelossfunction(E(theta,datavars$g22hatt,theoreticalgamma22Time,c,timelags),S22),
                onelossfunction(E(theta,datavars$g11hats,theoreticalgamma11Space,c,spacelags),S11space),
                onelossfunction(E(theta,datavars$g12hats,theoreticalgamma12Space,c,spacelags),S12space),
                onelossfunction(E(theta,datavars$g22hats,theoreticalgamma22Space,c,spacelags),S22space))
  return(t(lossvector)%*%(lossvector))
}
normed<-function(S){S/norm(S)}
getnewSs<-function(theta,datavars,c,timelags,spacelags){
  S11<-normed(getnewS(E(theta,datavars$g11hatt,theoreticalgamma11Time,c,timelags)))
  S12<-normed(getnewS(E(theta,datavars$g12hatt,theoreticalgamma12Time,c,timelags)))
  S22<-normed(getnewS(E(theta,datavars$g22hatt,theoreticalgamma22Time,c,timelags)))
  S11space<-normed(getnewS(E(theta,datavars$g11hats,theoreticalgamma11Space,c,spacelags)))
  S12space<-normed(getnewS(E(theta,datavars$g12hats,theoreticalgamma12Space,c,spacelags)))
  S22space<-normed(getnewS(E(theta,datavars$g22hats,theoreticalgamma22Space,c,spacelags)))
  return(list(S11=S11,S12=S12,S22=S22,S11space=S11space,S12space=S12space,S22space=S22space))
}
error<-function(theta,Y){sum(abs(theta-as.numeric(Y[[7]])))}

paraestimate<-function(Y,Tlag,Llag){
  y1<-Y[[1]]
  y2<-Y[[2]]
  Nx<-dim(y1)[2]
  Nt<-dim(y1)[1]
  Lx<-Y[[3]]
  Lt<-Y[[4]]
  c<-Y[[5]]
  #deltat<-Lt/(2*(Nt-1))
  #deltax<-Lx/(2*(Nx-1))
  deltat=Y[[6]]
  deltax=c*deltat
  if(c*deltat!=deltax){
    print("check dimesions")#Check must have c*deltat=deltax
  }
  if(missing(Tlag)){  #must be multiuple of 2*deltat
    Tlag=2*deltat
  }
  if(missing(Llag)){#must be multiuple of 2*deltax
    Llag=2*deltax
  }
  Cijhat<-function(i,j,Tlag,Llag){#i,j=1,2
    yi<-Y[[i]]
    yj<-Y[[j]]
    tpos<-(Tlag*(Nt-1))/Lt+1
    xpos<-(Llag*(Nx-1))/Lx+1
    return(cov(as.vector(yi[1:(Nt+1-tpos),1:(Nx+1-xpos)]),as.vector(yj[(tpos):Nt,(xpos):Nx])))
  }
  
  #mu11,mu21 estimate
  mu11hat<-1/Tlag*log(Cijhat(1,1,0,0)/Cijhat(1,1,Tlag,0))
  mu21hat<-1/Tlag*log(Cijhat(1,2,0,0)/Cijhat(1,2,Tlag,0))
  
  #lambda11 estimate
  alpha<-mu11hat*Llag/c
  gamma<-alpha*exp(alpha)/(exp(alpha)-1)
  lambda11hat<-((-lambertWm1(-Cijhat(1,1,0,Llag)/Cijhat(1,1,0,0)*exp(alpha)*gamma*exp(-gamma)))-gamma)/Llag
  
  #Estimation of k11^2
  k11hat2<-2*mu11hat*(mu11hat+c*lambda11hat)*Cijhat(1,1,0,0)/c
  
  #Estimation of lambda21
  torootl21<-function(l21){
    1/2*(Cijhat(1,2,0,0))*(mu11hat+mu21hat+c*(lambda11hat+l21))*exp(-(mu11hat+mu21hat)*Llag/c)*(
      (exp(-lambda11hat*Llag)+exp(-l21*Llag))/(mu11hat+mu21hat+c(lambda11hat+l21))
      +exp(-lambda11hat*Llag)*((exp((mu11hat+mu21hat+c(lambda11hat-l21))*Llag/(2*c))-1)/((mu11hat+mu21hat+c(lambda11hat-l21))))
      +exp(-l21*Llag)*((exp((mu11hat+mu21hat+c(-lambda11hat+l21))*Llag/(2*c))-1)/((mu11hat+mu21hat+c(-lambda11hat+l21))))
    )-Cijhat(1,2,0,Llag)
  }
  if(torootl21(0)>0&&torootl21(10)<0){
    lambda21hat<-uniroot(torootl21,c(0,10))$root #Not sure how to choose the bounds
  }
  else lambda21hat<-0
  
  #Estimation of k11k12
  k11k21hat<-Cijhat(1,2,0,0)*(mu11hat+mu21hat)*(mu11hat+mu21hat+c*(lambda11hat+lambda21hat))/(2*c)
  
  #Estimation of mu22
  k21hat2<-k11k21hat^2/k11hat2
  mu22hat<- -1/Tlag*log((2*mu21hat*(mu21hat+c*lambda21hat)*Cijhat(2,2,Tlag,0)-c*k21hat2*exp(-mu21hat*Tlag))/(2*mu21hat*(mu21hat+c*lambda21hat)*Cijhat(2,2,0,0)-c*k21hat2))
  
  #Estimation of lambda22
  torootl22<-function(l22){
    (Cijhat(2,2,0,0)-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*exp(-2*mu22hat*Llag/c)*(c/(Llag*mu22hat))*(exp(mu22hat*Llag/c)-1)*exp(-l22*Llag)*(l22*Llag+mu22hat*Llag/c*(exp(mu22hat*Llag/c)/(exp(mu22hat*Llag/c)-1)))-Cijhat(2,2,0,Llag)+(c*k21hat2)/(2*mu21hat^2)*exp(-(2*mu21hat/c+lambda21hat)*Llag)*(exp(mu21hat*Llag/c)-c*lambda21hat/(mu21hat+c*lambda21hat))
  }
  if(torootl22(0)>0&&torootl22(10)<0){
    lambda22hat<-uniroot(torootl22,c(0,10))$root
  }
  else lambda22hat<-0
  
  #estimaiton of k22
  k22hat2<-(Cijhat(2,2,0,0)-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*2*mu22hat*(mu22hat+c*lambda22hat)/c
  
  paraests<-list(Tlag=Tlag,Llag=Llag,k11est=(k11hat2)^0.5,k21est=k11k21hat/(k11hat2)^0.5,k22est=k22hat2^0.5,mu11est=mu11hat,mu21est=mu21hat,mu22est=mu22hat,lambda11est=lambda11hat,lambda21est=lambda21hat,lambda22est=lambda22hat)
  return(paraests)
}
getparasM5<-function(Y,K,initialestimate){
  c<-Y[[5]]
  deltat<-Y[[6]]
  timelags<-((1:K))*2*deltat
  spacelags<-((1:K))*2*deltat*c
  initialSs<-list(S11=diag(K),S12=diag(K),S22=diag(K),S11space=diag(K),S12space=diag(K),S22space=diag(K))
  datavars<-datavariogram(timelags,spacelags,Y)
  oldtheta<-optim(initialestimate,function(theta){totallossfunction(theta,datavars,initialSs,c,timelags,spacelags)},method="BFGS",control=list(maxit=10000))$par
  dif<-1
  oldtheta<-abs(oldtheta)
  while(dif>0.001){
    optimisation<-optim(oldtheta,function(theta){totallossfunction(theta,datavars,getnewSs(oldtheta,datavars,c,timelags,spacelags),c,timelags,spacelags)},method="BFGS")
    newtheta<-optimisation$par
    dif<-max(abs(oldtheta-newtheta))
    oldtheta<-abs(newtheta)
  }
  return(newtheta)
}



direct<-matrix(0,100,9)
lsmethod<-matrix(0,100,9)
lsmethodInit<-matrix(0,100,9)
for(i in 1:100){
  Y<-readRDS(paste("Output Fields/diffparas05depth3v",i,".rds",sep=""))
  direct[i,]<-as.numeric(paraestimate(Y)[3:11])
  lsmethod[i,]<-getparasM5(Y,10,rep(1,9))
  lsmethodInit[i,]<-getparasM5(Y,10,direct[i,])
}

violinplot(direct,Y[[7]])
violinplot(lsmethod,Y[[7]])
violinplot(lsmethodInit,Y[[7]])

differencematrix<-lsmethodInit-lsmethod

initialestimate<-paraestimate(Y)[3:11]
Y<-readRDS(paste("Output Fields/diffparas05depth3v",49,".rds",sep=""))
Y<-readRDS(paste("Output Fields/diffparas05depth3v",46,".rds",sep=""))
Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",10,".rds",sep=""))
