library(ggplot2)
library(pracma)
library(lamW)
library(ggpubr)
library(reshape2)

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
  return(c(g11t,g12t,g22t))
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
  return(list(timedata=rbind(g11hatt,g12hatt,g22hatt),spacedata=rbind(g11hats,g12hats,g22hats)))
}
aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}
paraestimate<-function(Y,Tlag,Llag,c){
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
}

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

ploterrorconvergence<-function(Y,NumberofLags,tol){
  paraslist<-twostepparas(Y,NumberofLags,tol)
  twostepdata<-as.data.frame(matrix(as.numeric(paraslist),nrow(paraslist),9))
  colnames(twostepdata)<-colnames(paraslist)
  errorsdataframe<-twostepdata-t(matrix(as.numeric(Y[[7]]),9,nrow(paraslist)))
  errorsdataframe$iter<-1:nrow(errorsdataframe)
  meltederrors<-melt(errorsdataframe,id="iter",variable="parameter")
  return(ggplot(meltederrors)+geom_line(aes(x=iter,y=value,colour = parameter)))
}


Y<-readRDS("Output Fields/highrespointone.rds1.rds")
Y<-readRDS("Output Fields/newparasfine3.rds")
Y<-readRDS("Output Fields/highrespointone.rds2.rds")
Y<-readRDS("Output Fields/highrespointone.rds3.rds")
Y<-readRDS("Output Fields/highrespointone.rds4.rds")
Y<-readRDS("Output Fields/highrespointone.rds5.rds")
NumberofLags<-10

paraslist<-twostepparas(Y,10,0.001)
as.numeric(Y[[7]])

ploterrorconvergence(Y,10,0.001)


#PLOT
twostepdata<-as.data.frame(matrix(as.numeric(paraslist),nrow(paraslist),9))
colnames(twostepdata)<-colnames(paraslist)
errorsdataframe<-twostepdata-t(matrix(as.numeric(Y[[7]]),9,15))

twostepdata$iter<-1:nrow(twostepdata)
meltedtwostepdata<-melt(twostepdata,id="iter")
ggplot(meltedtwostepdata)+geom_line(aes(x=iter,y=value,colour = variable))

errorsdataframe$iter<-1:nrow(errorsdataframe)
meltederrors<-melt(errorsdataframe,id="iter")
ggplot(meltederrors)+geom_line(aes(x=iter,y=value,colour = variable))







