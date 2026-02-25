library(ggplot2)
library(pracma)
library(lamW)
library(ggpubr)
setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")

oonefine1<-readRDS("Output Fields/highrespointone.rds1.rds")
Y<-oonefine1
paras<-Y[[7]]

gammaijTime<-function(i,j,Tlag,Y){
 yi<-Y[[i]]
 yj<-Y[[j]]
 deltat<-Y[[6]]
 Nt<-dim(yi)[1]
 Nx<-dim(yi)[2]
 1/2*mean((yi[(Tlag/(2*deltat)+1):Nt,1:Nx] - yj[1:(Nt-Tlag/(2*deltat)),1:Nx])^2)
}
gammaijSpace<-function(i,j,Llag,Y){
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  1/2*mean((yi[1:Nt,(1+Llag/(2*deltat*c)):Nx] - yj[1:Nt,1:(Nx-Llag/(2*deltat*c))])^2)
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
  return((c*k11^2/(2*mu11*(mu11+c*lambda11))+c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22)))*0.5-(2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21)))*exp(-mu11*Tlag))
}
theoreticalgamma21Time<-function(Tlag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return((c*k11^2/(2*mu11*(mu11+c*lambda11))+c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22)))*0.5-(2*c*k11*k21)/((mu11+mu21)*(mu11+mu21+c*(lambda11+lambda21)))*exp(-mu21*Tlag))
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
  return((c*k11^2/(2*mu11*(mu11+c*lambda11))+c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22)))*0.5
         -c*k11*k21*exp(-(mu11+mu21)*(Llag/c))/(mu11+mu21)*(
           (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*lambda11+c*lambda21)+
             exp(-lambda11*Llag)*(exp((mu11+mu21+c*lambda11-c*lambda21)*Llag/(2*c))-1)/(mu11+mu21+c*lambda11-c*lambda21)+
             exp(-lambda21*Llag)*(exp((mu11+mu21-c*lambda11+c*lambda21)*Llag/(2*c))-1)/(mu11+mu21-c*lambda11+c*lambda21)
         )
         )
}
theoreticalgamma21Space<-function(Llag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  return((c*k11^2/(2*mu11*(mu11+c*lambda11))+c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22)))*0.5
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

theoreticalvariograms<-function(parameters,timelags,spacelags,c){
  if(length(timelags)!=length(spacelags)){
    print("Need same number of space and time lags")
  }
  NumberofLags<-length(timelags)
  g11t<-rep(0,NumberofLags)
  g11s<-rep(0,NumberofLags)
  g12t<-rep(0,NumberofLags)
  g12s<-rep(0,NumberofLags)
  g21t<-rep(0,NumberofLags)
  g21s<-rep(0,NumberofLags)
  g22t<-rep(0,NumberofLags)
  g22s<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    g11t[i]<-theoreticalgamma11Time(timelags[i],parameters,c)
    g12t[i]<-theoreticalgamma12Time(timelags[i],parameters,c)
    g21t[i]<-theoreticalgamma21Time(timelags[i],parameters,c)
    g22t[i]<-theoreticalgamma22Time(timelags[i],parameters,c)
    g11s[i]<-theoreticalgamma11Space(spacelags[i],parameters,c)
    g12s[i]<-theoreticalgamma12Space(spacelags[i],parameters,c)
    g21s[i]<-theoreticalgamma21Space(spacelags[i],parameters,c)
    g22s[i]<-theoreticalgamma22Space(spacelags[i],parameters,c)
  }
  return(rbind(g11t,g12t,g21t,g22t,g11s,g12s,g21s,g22s))
}

datavariogram<-function(timelags,spacelags,Y){
  if(length(timelags)!=length(spacelags)){
    print("Need same number of space and time lags")
  }
  NumberofLags<-length(timelags)
  g11hatt<-rep(0,NumberofLags)
  g11hats<-rep(0,NumberofLags)
  g12hatt<-rep(0,NumberofLags)
  g12hats<-rep(0,NumberofLags)
  g21hatt<-rep(0,NumberofLags)
  g21hats<-rep(0,NumberofLags)
  g22hatt<-rep(0,NumberofLags)
  g22hats<-rep(0,NumberofLags)
  for(i in 1:NumberofLags){
    g11hatt[i]<-gammaijTime(1,1,timelags[i],Y)
    g12hatt[i]<-gammaijTime(1,2,timelags[i],Y)
    g21hatt[i]<-gammaijTime(2,1,timelags[i],Y)
    g22hatt[i]<-gammaijTime(2,2,timelags[i],Y)
    g11hats[i]<-gammaijSpace(1,1,spacelags[i],Y)
    g12hats[i]<-gammaijSpace(1,2,spacelags[i],Y)
    g21hats[i]<-gammaijSpace(2,1,spacelags[i],Y)
    g22hats[i]<-gammaijSpace(2,2,spacelags[i],Y)
  }
  return(rbind(g11hatt,g12hatt,g21hatt,g22hatt,g11hats,g12hats,g21hats,g22hats))
}

aslist<-function(vecparameters){return(list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9]))}
errorfunction<-function(timelags,spacelags,datavars,vecparameters){
  #parameters<-list(k11=vecparameters[1],k21=vecparameters[2],k22=vecparameters[3],mu11=vecparameters[4],mu21=vecparameters[5],mu22=vecparameters[6],lambda11=vecparameters[7],lambda21=vecparameters[8],lambda22=vecparameters[9])
  parameters<-aslist(vecparameters)
  c<-Y[[5]]
  return(datavars-as.vector(theoreticalvariograms(parameters,timelags,spacelags,c)))
}

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

getestimateGN<-function(Y,NumberofLags){
  c<-Y[[5]]
  deltat<-Y[[6]]
  timelags<-((1:NumberofLags)-1)*2*deltat
  spacelags<-((1:NumberofLags)-1)*2*deltat*c
  initialvalue<-as.numeric(paraestimate(Y,c=c))[3:11]
  datavars<-as.vector(datavariogram(timelags,spacelags,Y))
  return(gaussNewton(initialvalue,function(vecparameters) errorfunction(timelags,spacelags,datavars,vecparameters))$xs)
}

getestimateGN(Y,30)



#Plotting to be redone
Y<-readRDS("Output Fields/newparasfine3.rds")

NumberofLags<-10
c<-Y[[5]]
deltat<-Y[[6]]
timelags<-((1:NumberofLags)-1)*2*deltat
spacelags<-((1:NumberofLags)-1)*2*deltat*c

trueparas<-Y[[7]]
estimatedparas<-aslist(getestimateGN(Y,NumberofLags))


theoreticalvariograms(estimatedparas,timelags,spacelags,c)
datavariogram(timelags,spacelags,Y)
variogramestdata<-data.frame(timelags,spacelags,t(datavariogram(timelags,spacelags,Y)),t(theoreticalvariograms(estimatedparas,timelags,spacelags,c)))
truevars<-data.frame(t(theoreticalvariograms(trueparas,timelags,spacelags,c)))
p11t<-ggplot(variogramestdata,aes(timelags,g11t))+geom_line()+geom_point(aes(timelags,g11hatt))+geom_line(aes(timelags,truevars$g11t),col="blue")+xlab("d_t")+ylab("gamma11T(d_t)")
p12t<-ggplot(variogramestdata,aes(timelags,g12t))+geom_line()+geom_point(aes(timelags,g12hatt))+geom_line(aes(timelags,truevars$g12t),col="blue")+xlab("d_t")+ylab("gamma12T(d_t)")
p21t<-ggplot(variogramestdata,aes(timelags,g21t))+geom_line()+geom_point(aes(timelags,g21hatt))+geom_line(aes(timelags,truevars$g21t),col="blue")+xlab("d_t")+ylab("gamma21T(d_t)")
p22t<-ggplot(variogramestdata,aes(timelags,g22t))+geom_line()+geom_point(aes(timelags,g22hatt))+geom_line(aes(timelags,truevars$g22t),col="blue")+xlab("d_t")+ylab("gamma22T(d_t)")
p11s<-ggplot(variogramestdata,aes(timelags,g11s))+geom_line()+geom_point(aes(timelags,g11hats))+geom_line(aes(timelags,truevars$g11s),col="blue")+xlab("d_x")+ylab("gamma11S(d_x)")
p12s<-ggplot(variogramestdata,aes(timelags,g12s))+geom_line()+geom_point(aes(timelags,g12hats))+geom_line(aes(timelags,truevars$g12s),col="blue")+xlab("d_x")+ylab("gamma12S(d_x)")
p21s<-ggplot(variogramestdata,aes(timelags,g21s))+geom_line()+geom_point(aes(timelags,g21hats))+geom_line(aes(timelags,truevars$g21s),col="blue")+xlab("d_x")+ylab("gamma21S(d_x)")
p22s<-ggplot(variogramestdata,aes(timelags,g22s))+geom_line()+geom_point(aes(timelags,g22hats))+geom_line(aes(timelags,truevars$g22s),col="blue")+xlab("d_x")+ylab("gamma22S(d_x)")

ggarrange(p11t,p12t,p21t,p22t)
ggarrange(p11s,p12s,p21s,p22s)



