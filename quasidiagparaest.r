Y<-readRDS(paste("Output Fields/diffparas05depth3v",1,".rds",sep=""))
Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",1,".rds",sep=""))

K<-10
c<-Y[[5]]
deltat<-Y[[6]]
timelags<-((1:K))*2*deltat
spacelags<-((1:K))*2*deltat*c

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
  gammahat-theoreticalvariograms(variogram,aslist(theta),lags,c)
}
onelossfunction<-function(E,S){
  t(rowSums(E))%*%S%*%rowSums(E)
}
getnewS<-function(E){
  return(solve(E%*%t(E)))
}
totallossfunction<-function(theta,datavars,Ss){
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
getnewSs<-function(theta,datavars){
  S11<-getnewS(E(theta,datavars$g11hatt,theoreticalgamma11Time,c,timelags))
  S12<-getnewS(E(theta,datavars$g12hatt,theoreticalgamma12Time,c,timelags))
  S22<-getnewS(E(theta,datavars$g22hatt,theoreticalgamma22Time,c,timelags))
  S11space<-getnewS(E(theta,datavars$g11hats,theoreticalgamma11Space,c,spacelags))
  S12space<-getnewS(E(theta,datavars$g12hats,theoreticalgamma12Space,c,spacelags))
  S22space<-getnewS(E(theta,datavars$g22hats,theoreticalgamma22Space,c,spacelags))
  return(list(S11=S11,S12=S12,S22=S22,S11space=S11space,S12space=S12space,S22space=S22space))
}
initialSs<-list(S11=diag(K),S12=diag(K),S22=diag(K),S11space=diag(K),S12space=diag(K),S22space=diag(K))
datavars<-datavariogram(timelags,spacelags,Y)

theta1<-optim(rep(0.4,9),function(theta){totallossfunction(theta,datavars,initialSs)})$par
theta2<-optim(rep(0.4,9),function(theta){totallossfunction(theta,datavars,getnewSs(theta1,datavars))})$par
theta3<-optim(rep(0.4,9),function(theta){totallossfunction(theta,datavars,getnewSs(theta2,datavars))})$par
theta4<-optim(rep(0.4,9),function(theta){totallossfunction(theta,datavars,getnewSs(theta3,datavars))})$par



levelplot(getnewSs(theta1,datavars)$S22space)

optim(c(1,1,1,1,1,1,1,1,1),function(theta){onelossfunction(E(theta,datavariogram(timelags,spacelags,Y)$g11hatt,theoreticalgamma11Time,c,timelags),S)})

theoreticalvariograms(theoreticalgamma11Time,Y[[7]],timelags,Y[[5]])
datavariogram(timelags,spacelags,Y)$g11hatt-theoreticalvariograms(theoreticalgamma11Time,Y[[7]],timelags,Y[[5]])
E11<-E(c(1,1,1,1,1,1,1,1,1),datavariogram(timelags,spacelags,Y)$g11hatt,theoreticalgamma11Time,c,timelags)
rowSums(E11)

E(aslist(c(1,1,1,1,1,1,1,1,1)),datavariogram(timelags,spacelags,Y)$g11hatt,theoreticalgamma11Time,c,timelags)




theoreticalvariogramsTime<-function(parameters,timelags,c){
  NumberofLags<-length(timelags)
  g11t<-rep(0,NumberofLags)
  g12t<-rep(0,NumberofLags)
  g22t<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    g11t[i]<-theoreticalgamma11Time(timelags[i],parameters,c)
    g12t[i]<-theoreticalgamma12Time(timelags[i],parameters,c)
    g22t[i]<-theoreticalgamma22Time(timelags[i],parameters,c)
  }
  return(list(g11t=g11t,g12t=g12t,g22t=g22t))
}
theoreticalvariogramsSpace<-function(parameters,spacelags,c){
  NumberofLags<-length(spacelags)
  g11s<-rep(0,NumberofLags)
  g12s<-rep(0,NumberofLags)
  g22s<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    g11s[i]<-theoreticalgamma11Space(spacelags[i],parameters,c)
    g12s[i]<-theoreticalgamma12Space(spacelags[i],parameters,c)
    g22s[i]<-theoreticalgamma22Space(spacelags[i],parameters,c)
  }
  return(c(g11s,g12s,g22s))
}
Etime<-function(theta){datavarstime-theoreticalvariogramsTime(aslist(theta),timelags,1)}
Espace<-function(theta){datavarsspace-theoreticalvariogramsTime(aslist(theta),spacelags,1)}
WeightmatTime<-function(theta){
  Nx=dim(Etime(theta))[2]
  solve(1/Nx*Etime(theta)%*%t(Etime(theta)))}
WeightmatSpace<-function(theta){
  Nt=dim(Etime(theta))[2]
  solve(1/Nt*Espace(theta)%*%t(Espace(theta)))}
optimise<-function(initialvalue,datavars,timelags,spacelags,WT,WS,c){
  datavarstime<-datavars$timedata
  datavarsspace<-datavars$spacedata
  datatimerow<-rowMeans(datavarstime)
  dataspacerow<-rowMeans(datavarsspace)
  Etime<-function(theta,c){datavarstime-theoreticalvariogramsTime(aslist(theta),timelags,c)}
  Espace<-function(theta,c){datavarsspace-theoreticalvariogramsTime(aslist(theta),spacelags,c)}
  WeightmatTime<-function(theta,c){
    Nx=dim(Etime(theta,c))[2]
    solve(1/Nx*Etime(theta,c)%*%t(Etime(theta,c)))}
  WeightmatSpace<-function(theta,c){
    Nt=dim(Etime(theta,c))[2]
    solve(1/Nt*Espace(theta,c)%*%t(Espace(theta,c)))}
  if(missing(WT)){WT<-WeightmatTime(initialvalue,c)}
  if(missing(WS)){WS<-WeightmatSpace(initialvalue,c)}
  as.numeric(optim(initialvalue,function(vecparameters) t(dataspacerow-theoreticalvariogramsSpace(aslist(vecparameters),spacelags,c))%*%WS%*%(dataspacerow-theoreticalvariogramsSpace(aslist(vecparameters),spacelags,c))+t(datatimerow-theoreticalvariogramsTime(aslist(vecparameters),timelags,c))%*%WT%*%(datatimerow-theoreticalvariogramsTime(aslist(vecparameters),timelags,c)))$par)
}#optimisation function

twostepparas<-function(Y,NumberofLags,tol){
  c<-Y[[5]]
  deltat<-Y[[6]]
  timelags<-((1:NumberofLags))*2*deltat
  spacelags<-((1:NumberofLags))*2*deltat*c
  datavars<-datavariogram(timelags,spacelags,Y)
  
  initialvalue<-as.numeric(paraestimate(Y,c=c))[3:11]
  theta0<-optimise(initialvalue,datavars,timelags,spacelags,WT=ones(NumberofLags*3),WS=ones(NumberofLags*3),c=c)
  parameterEstimateList<-rbind(aslist(initialvalue),aslist(theta0))
  diff<-sum(abs(initialvalue-theta0))
  theta<-theta0
  while(diff>tol){
    thetaNew<-optimise(theta,datavars,timelags,spacelags,c=c)
    diff<-sum(abs(thetaNew-theta))
    theta<-thetaNew
    parameterEstimateList<-rbind(parameterEstimateList,aslist(thetaNew))
  }
  return(parameterEstimateList)
}