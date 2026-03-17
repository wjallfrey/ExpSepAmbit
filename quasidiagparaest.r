Y<-readRDS(paste("Output Fields/diffparas05depth3v",10,".rds",sep=""))
Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",1,".rds",sep=""))


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
normed<-function(S){S/norm(S,type="F")}
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

getparasM5<-function(Y,K){
  c<-Y[[5]]
  deltat<-Y[[6]]
  timelags<-((1:K))*2*deltat
  spacelags<-((1:K))*2*deltat*c
  initialSs<-list(S11=diag(K),S12=diag(K),S22=diag(K),S11space=diag(K),S12space=diag(K),S22space=diag(K))
  datavars<-datavariogram(timelags,spacelags,Y)
  oldtheta<-optim(rep(1,9),function(theta){totallossfunction(theta,datavars,initialSs,c,timelags,spacelags)},method="BFGS",control=list(maxit=10000))$par
  dif<-1
  while(dif>0.001){
    optimisation<-optim(oldtheta,function(theta){totallossfunction(theta,datavars,getnewSs(oldtheta,datavars,c,timelags,spacelags),c,timelags,spacelags)},method="BFGS")
    newtheta<-optimisation$par
    dif<-max(abs(oldtheta-newtheta))
    oldtheta<-abs(newtheta)
  }
  return(newtheta)
}
getparasM5logexp<-function(Y,K){
  c<-Y[[5]]
  deltat<-Y[[6]]
  timelags<-((1:K))*2*deltat
  spacelags<-((1:K))*2*deltat*c
  initialSs<-list(S11=diag(K),S12=diag(K),S22=diag(K),S11space=diag(K),S12space=diag(K),S22space=diag(K))
  datavars<-datavariogram(timelags,spacelags,Y)
  oldlogtheta<-optim(rep(0,9),function(logtheta){totallossfunction(exp(logtheta),datavars,initialSs,c,timelags,spacelags)},method="BFGS",control=list(maxit=10000))$par
  dif<-1
  while(dif>0.001){
    #optimisation<-optim(oldtheta,function(theta){totallossfunction(theta,datavars,getnewSs(oldtheta,datavars,c,timelags,spacelags),c,timelags,spacelags)},method="BFGS")
    optimisationexp<-optim(oldlogtheta,function(logtheta){totallossfunction(exp(logtheta),datavars,getnewSs(exp(oldlogtheta),datavars,c,timelags,spacelags),c,timelags,spacelags)},method="BFGS")
    newlogtheta<-optimisationexp$par
    dif<-max(abs(oldlogtheta-newlogtheta))
    oldlogtheta<-newlogtheta
  }
  return(exp(newlogtheta))
}

#look at autodiff function

#THE BEST ONE yayviolin
paraslist04parA3<-matrix(0,1000,9)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",i,".rds",sep=""))
  pars<-getparasM5(Y,10)
  paraslist04parA3[i,]<-pars
  if(i%%100==0){print(paste("at",i%/%10,"%"))}
}
violinplot(paraslist04parA3,Y[[7]])

#NOT YET RUN
paraslist04parA3expmethod<-matrix(0,1000,9)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",i,".rds",sep=""))
  pars<-getparasM5logexp(Y,10)
  paraslist04parA3expmethod[i,]<-pars
  if(i%%100==0){print(paste("at",i%/%10,"%",Sys.time()))}
}
violinplot(paraslist04parA3,Y[[7]])

#GOOD BUT SOME ZEROS - check code and multiple initialising points
paraslistfineA3<-matrix(0,100,9)
for(i in 1:100){
  Y<-readRDS(paste("Output Fields/fineambit3par04res01depth3v",i,".rds",sep=""))
  pars<-getparasM5logexp(Y,10)
  paraslistfineA3[i,]<-pars
  print(pars)
}
violinplot(paraslistfineA3,Y[[7]])

#WAY OFF ON A FEW OF THEM
paraslistdiffpar<-matrix(0,1000,9)
#errors<-rep(0,1000)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/diffparas05depth3v",i,".rds",sep=""))
  optimedpars<-getparasM5(Y,10)
  paraslistdiffpar[i,]<-optimedpars
  #errors[i]<-error(optimedpars,Y)
  if(i%%100==0){print(paste("at",i%/%10,"%"))}
}
violinplot(paraslistdiffpar[1:5,],Y[[7]])

#WILL THIS BE BETTER - ALSO NOT YET RUN
paraslistdiffparexpmethod<-matrix(0,1000,9)
#errors<-rep(0,1000)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/diffparas05depth3v",i,".rds",sep=""))
  optimedpars<-getparasM5logexp(Y,10)
  paraslistdiffparexpmethod[i,]<-optimedpars
  #errors[i]<-error(optimedpars,Y)
  if(i%%100==0){print(paste("at",i%/%10,"%"))}
}
violinplot(paraslistdiffparexpmethod,Y[[7]])




