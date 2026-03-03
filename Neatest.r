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

theoreticalgamma11Time<-function(Tlag,parameters,c){C11Time(0,parameters,c)-C11Time(Tlag,parameters,c)}
theoreticalgamma12Time<-function(Tlag,parameters,c){C12Time(0,parameters,c)-0.5*(C12Time(Tlag,parameters,c)+C21Time(Tlag,parameters,c))}
theoreticalgamma22Time<-function(Tlag,parameters,c){C22Time(0,parameters,c)-C22Time(Tlag,parameters,c)}
theoreticalgamma11Space<-function(Llag,parameters,c){C11Space(0,parameters,c)-C11Space(Llag,trueparas,c)}
theoreticalgamma12Space<-function(Llag,parameters,c){C12Space(0,parameters,c)-C12Space(Llag,trueparas,c)}
theoreticalgamma22Space<-function(Llag,parameters,c){C22Space(0,parameters,c)-C22Space(Llag,trueparas,c)}

aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}

ETime<-function(i,j,timelags,parameters,Y,c){
  
}

#No weights
E11Time<-function(parameters,Y,timelags){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  #deltat<-Y[[6]]
  #timelags<-((1:NumberofLags))*2*deltat
  Nx<-dim(Y[[1]])[2]
  E11timemat<-matrix(0,NumberofLags,Nx)
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E11timemat[i,]<-gammaijTimeHAT(1,1,Tlag,Y)-theoreticalgamma11Time(Tlag,parameters,c)
  }
  return(E11timemat)
}
E12Time<-function(parameters,Y,timelags){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  E12timemat<-matrix(0,NumberofLags,Nx)
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E12timemat[i,]<-gammaijTimeHAT(1,2,Tlag,Y)-theoreticalgamma12Time(Tlag,parameters,c)
  }
  return(E12timemat)
}
E22Time<-function(parameters,Y,timelags){
  NumberofLags<-length(timelags)
  c<-Y[[5]]
  Nx<-dim(Y[[1]])[2]
  E22timemat<-matrix(0,NumberofLags,Nx)
  for(i in 1:NumberofLags){
    Tlag<-timelags[i]
    E22timemat[i,]<-gammaijTimeHAT(2,2,Tlag,Y)-theoreticalgamma22Time(Tlag,parameters,c)
  }
  return(E22timemat)
}

optim(c(1,1,1),function(a){sum((E11Time(aslist(c(a[1],0,0,a[2],0,0,a[3],0,0)),Y,timelags))^2)})


#Testing workspace
c<-Y[[5]]
deltat<-Y[[6]]
NumberofLags<-15
timelags<-((1:NumberofLags))*2*deltat
spacelags<-((1:NumberofLags))*2*deltat*c
View(E11Time(Y[[7]],Y,timelags))

levelplot(solve(((E11Time(Y[[7]],Y,timelags))%*%t(E11Time(Y[[7]],Y,timelags)))))
levelplot(solve(((E12Time(Y[[7]],Y,timelags))%*%t(E12Time(Y[[7]],Y,timelags)))))

rowSums(E11Time(Y[[7]],Y,timelags)^2)
