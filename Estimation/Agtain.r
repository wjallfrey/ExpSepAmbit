#Unnormalised cross-variogram

Y<-readRDS("Output Fields/newparasfine3.rds")


###JUST FOR TIME
gammaijTime<-function(i,j,Tlag,Y){
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  gam<-colMeans(0.5*(yi[(Tlag/(2*deltat)+1):Nt,1:Nx]-yi[1:(Nt-Tlag/(2*deltat)),1:Nx])*(yj[(Tlag/(2*deltat)+1):Nt,1:Nx]-yj[1:(Nt-Tlag/(2*deltat)),1:Nx]))
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

deltat<-Y[[6]]
NumberofLags<-20
timelags<-((1:NumberofLags))*2*deltat

theoreticalvariograms<-function(parameters,timelags,c){
  NumberofLags<-length(timelags)
  g11t<-rep(0,NumberofLags)
  g12t<-rep(0,NumberofLags)
  g22t<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    g11t[i]<-theoreticalgamma11Time(timelags[i],parameters,c)
    g12t[i]<-theoreticalgamma12Time(timelags[i],parameters,c)
    g22t[i]<-theoreticalgamma22Time(timelags[i],parameters,c)
  }
  return(c(g11t,g12t,g22t))
}

datavariogram<-function(timelags,Y){
  NumberofLags<-length(timelags)
  Nx<-dim(Y[[1]])[2]
  g11hatt<-matrix(0,NumberofLags,Nx)
  g12hatt<-matrix(0,NumberofLags,Nx)
  g22hatt<-matrix(0,NumberofLags,Nx)
  for(i in 1:NumberofLags){
    g11hatt[i,]<-gammaijTime(1,1,timelags[i],Y)
    g12hatt[i,]<-gammaijTime(1,2,timelags[i],Y)
    g22hatt[i,]<-gammaijTime(2,2,timelags[i],Y)
  }
  return(rbind(g11hatt,g12hatt,g22hatt))
}

#compare rowMeans(datavar..) to theoreticalvariograms
dim(datavariogram(timelags,Y))

aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}

datas<-rowMeans(datavariogram(timelags,Y))

thetanought<-optim(initialvalue,function(vecparameters) sum(datas-theoreticalvariograms(aslist(vecparameters),timelags,1))^2)$par

E<-datavariogram(timelags,Y)-theoreticalvariograms(aslist(thetanought),timelags,1)
Sigma<-E%*%t(E)
Weightmat<-solve(Sigma)

thetaone<-optim(initialvalue,function(vecparameters) t(datas-theoreticalvariograms(aslist(vecparameters),timelags,1))%*%Weightmat%*%((datas-theoreticalvariograms(aslist(vecparameters),timelags,1))))$par

sum((thetanought-as.numeric(Y[[7]]))^2)
sum((thetaone-as.numeric(Y[[7]]))^2)
sum((thetatwo-as.numeric(Y[[7]]))^2)

E<-datavariogram(timelags,Y)-theoreticalvariograms(aslist(thetaone),timelags,1)
Sigma<-E%*%t(E)
Weightmat<-solve(Sigma)

thetatwo<-optim(thetaone,function(vecparameters) t(datas-theoreticalvariograms(aslist(vecparameters),timelags,1))%*%Weightmat%*%((datas-theoreticalvariograms(aslist(vecparameters),timelags,1))))$par
