library(lamW)

Y<-readRDS("Output Fields/multitest1.rds")

trueparas<-Y[[7]]

actualmoments<-function(Tlag,Llag,Y){
  Lx<-Y[[3]]
  Lt<-Y[[4]]
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  Cijhat<-function(i,j,Tlag,Llag){#i,j=1,2
    yi<-Y[[i]]
    yj<-Y[[j]]
    tpos<-(Tlag*(Nt-1))/Lt+1
    xpos<-(Llag*(Nx-1))/Lx+1
    return(cov(as.vector(yi[1:(Nt+1-tpos),1:(Nx+1-xpos)]),as.vector(yj[(tpos):Nt,(xpos):Nx])))
  }
  return(list(c1100=Cijhat(1,1,0,0),c11T0=Cijhat(1,1,Tlag,0),c110L=Cijhat(1,1,0,Llag),c1200=Cijhat(1,2,0,0),c12T0=Cijhat(1,2,Tlag,0),c120L=Cijhat(1,2,0,Llag),c21T0=Cijhat(2,1,Tlag,0),c210L=Cijhat(2,1,0,Llag),c2200=Cijhat(2,2,0,0),c22T0=Cijhat(2,2,Tlag,0),c220L=Cijhat(2,2,0,Llag)))
}

theoreticalmoments<-function(Tlag,Llag,parameters,c){
  k11<-parameters$k11
  k21<-parameters$k21
  k22<-parameters$k22
  mu11<-parameters$mu11
  mu21<-parameters$mu21
  mu22<-parameters$mu22
  lambda11<-parameters$lambda11
  lambda21<-parameters$lambda21
  lambda22<-parameters$lambda22
  c1100<-c*k11^2/(2*mu11*(mu11+c*lambda11))
  c11T0<-c*k11^2/(2*mu11*(mu11+c*lambda11))*exp(-Tlag*mu11)
  c110L<-c*k11^2/(2*mu11^2)*exp(-2*mu11*Llag/c-lambda11*Llag)*(exp(mu11*Llag/c)-c*lambda11/(mu11+c*lambda11))
  c1200<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))
  c12T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu21*Tlag)
  c120L<-c*k11*k21*exp(-(mu11+mu21)*Llag/c)/(mu11+mu21)*(
    (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*(lambda11+lambda21))+
      exp(-lambda11*Llag)/(mu11+mu21+c*(lambda11-lambda21))*(exp((mu11+mu21+c*(lambda11-lambda21))*Llag/(2*c))-1)+
      exp(-lambda21*Llag)/(mu11+mu21+c*(-lambda11+lambda21))*(exp((mu11+mu21+c*(-lambda11+lambda21))*Llag/(2*c))-1)
  )
  c21T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu11*Tlag)
  c210L<-c*k11*k21*exp(-(mu11+mu21)*Llag/c)/(mu11+mu21)*(
    (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*(lambda11+lambda21))+
      exp(-lambda11*Llag)/(mu11+mu21+c*(lambda11-lambda21))*(exp((mu11+mu21+c*(lambda11-lambda21))*Llag/(2*c))-1)+
      exp(-lambda21*Llag)/(mu11+mu21+c*(-lambda11+lambda21))*(exp((mu11+mu21+c*(-lambda11+lambda21))*Llag/(2*c))-1)
  )
  c2200<-c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22))
  c22T0<-c*k21^2/(2*mu21*(mu21+c*lambda21))*exp(-mu21*Tlag)+c*k22^2/(2*mu22*(mu22+c*lambda22))*exp(-mu22*Tlag)
  c220L<-c*k21^2/(2*mu21^2)*exp(-2*mu21*Llag/c-lambda21*Llag)*(exp(mu21*Llag/c)-c*lambda21/(mu21+c*lambda21))+
    c*k22^2/(2*mu22^2)*exp(-2*mu22*Llag/c-lambda22*Llag)*(exp(mu22*Llag/c)-c*lambda22/(mu22+c*lambda22))
  return(list(c1100=c1100,c11T0=c11T0,c110L=c110L,c1200=c1200,c12T0=c12T0,c120L=c120L,c21T0=c21T0,c210L=c210L,c2200=c2200,c22T0=c22T0,c220L=c220L))
}
paraestfrommoments<-function(Tlag,Llag,moments,c){
  C1100=moments$c1100
  C11T0=moments$c11T0
  C110L=moments$c110L
  C1200=moments$c1200
  C12T0=moments$c12T0
  C120L=moments$c120L
  C21T0=moments$c21T0
  C210L=moments$c210L
  C2200=moments$c2200
  C22T0=moments$c22T0
  C220L=moments$c220L
  #mu11,mu21 estimate
  mu11hat<-1/Tlag*log(C1100/C11T0)
  mu21hat<-1/Tlag*log(C1200/C12T0)
  
  #lambda11 estimate
  alpha<-mu11hat*Llag/c
  gamma<-alpha*exp(alpha)/(exp(alpha)-1)
  lambda11hat<-((-lambertWm1(-C110L/C1100*exp(alpha)*gamma*exp(-gamma)))-gamma)/Llag
  
  #Estimation of k11^2
  k11hat2<-2*mu11hat*(mu11hat+c*lambda11hat)*C1100/c
  
  #Estimation of lambda21
  torootl21<-function(l21){
    1/2*(C1200)*(mu11hat+mu21hat+c*(lambda11hat+l21))*exp(-(mu11hat+mu21hat)*Llag/c)*(
      (exp(-lambda11hat*Llag)+exp(-l21*Llag))/(mu11hat+mu21hat+c(lambda11hat+l21))
      +exp(-lambda11hat*Llag)*((exp((mu11hat+mu21hat+c(lambda11hat-l21))*Llag/(2*c))-1)/((mu11hat+mu21hat+c(lambda11hat-l21))))
      +exp(-l21*Llag)*((exp((mu11hat+mu21hat+c(-lambda11hat+l21))*Llag/(2*c))-1)/((mu11hat+mu21hat+c(-lambda11hat+l21))))
    )-C120L
  }
  lambda21hat<-uniroot(torootl21,c(0,10))$root #Not sure how to choose the bounds
  
  #Estimation of k11k12
  k11k21hat<-C1200*(mu11hat+mu21hat)*(mu11hat+mu21hat+c*(lambda11hat+lambda21hat))/(2*c)
  
  #Estimation of mu22
  k21hat2<-k11k21hat^2/k11hat2
  mu22hat<- -1/Tlag*log((2*mu21hat*(mu21hat+c*lambda21hat)*C22T0-c*k21hat2*exp(-mu21hat*Tlag))/(2*mu21hat*(mu21hat+c*lambda21hat)*C2200-c*k21hat2))
  
  #Estimation of lambda22
  torootl22<-function(l22){
    (C2200-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*exp(-2*mu22hat*Llag/c)*(c/(Llag*mu22hat))*(exp(mu22hat*Llag/c)-1)*exp(-l22*Llag)*(l22*Llag+mu22hat*Llag/c*(exp(mu22hat*Llag/c)/(exp(mu22hat*Llag/c)-1)))-C220L+(c*k21hat2)/(2*mu21hat^2)*exp(-(2*mu21hat/c+lambda21hat)*Llag)*(exp(mu21hat*Llag/c)-c*lambda21hat/(mu21hat+c*lambda21hat))
  }
  lambda22hat<-uniroot(torootl22,c(0,10))$root
  
  #estimaiton of k22
  k22hat2<-(C2200-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*2*mu22hat*(mu22hat+c*lambda22hat)/c
  
  paraests<-list(Tlag=Tlag,Llag=Llag,k11est=(k11hat2)^0.5,k21est=k11k21hat/(k11hat2)^0.5,k22est=k22hat2^0.5,mu11est=mu11hat,mu21est=mu21hat,mu22est=mu22hat,lambda11est=lambda11hat,lambda21est=lambda21hat,lambda22est=lambda22hat)
  return(paraests)
}

errorMoment<-function(dataMoments,theoryMoments,absolute,max){
  if(missing(max)){
    max=F
  }
  if(missing(absolute)){
    absolute=T
  }
  errors<-(as.numeric(dataMoments)-as.numeric(theoryMoments))/as.numeric(dataMoments)
  if(max){return(max(abs(errors)))}
  else if(absolute){return(abs(errors))}
}

actualmoments(1,1,Y)
paraestfrommoments(1,1,actualmoments(1,1,Y),1)

theoreticalmoments(1,1,paraestfrommoments(1,1,actualmoments(1,1,Y),1),1)


rbind(actualmoments(1,1,Y),theoreticalmoments(1,1,paraestfrommoments(1,1,actualmoments(1,1,Y),1),1))


library(pracma)
gaussNewton(c(0,0),function(k) moments(0.1,2,k[1],k[2]))

c<-1
deltat<-0.1
deltax<-c*deltat

Tlags<-1:10*(2*deltat)
Llags<-1:10*(2*deltax)

moments<-function(Tlags,Llags,parameters){
  output<-rep(0,0)
  k11<-parameters[1]
  k21<-parameters[2]
  k22<-parameters[3]
  mu11<-parameters[4]
  mu21<-parameters[5]
  mu22<-parameters[6]
  lambda11<-parameters[7]
  lambda21<-parameters[8]
  lambda22<-parameters[9]
  c1100<-c*k11^2/(2*mu11*(mu11+c*lambda11))
  c1200<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))
  c2200<-c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22))
  output<-c(c1100,c1200,c2200)
  for(Tlag in Tlags){
    c11T0<-c*k11^2/(2*mu11*(mu11+c*lambda11))*exp(-Tlag*mu11)
    c12T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu21*Tlag)
    c21T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu11*Tlag)
    c22T0<-c*k21^2/(2*mu21*(mu21+c*lambda21))*exp(-mu21*Tlag)+c*k22^2/(2*mu22*(mu22+c*lambda22))*exp(-mu22*Tlag)
    output<-c(output,c11T0,c12T0,c21T0,c22T0)
  }
  for(Llag in Llags){
    c110L<-c*k11^2/(2*mu11^2)*exp(-2*mu11*Llag/c-lambda11*Llag)*(exp(mu11*Llag/c)-c*lambda11/(mu11+c*lambda11))
    c120L<-c*k11*k21*exp(-(mu11+mu21)*Llag/c)/(mu11+mu21)*(
      (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*(lambda11+lambda21))+
        exp(-lambda11*Llag)/(mu11+mu21+c*(lambda11-lambda21))*(exp((mu11+mu21+c*(lambda11-lambda21))*Llag/(2*c))-1)+
        exp(-lambda21*Llag)/(mu11+mu21+c*(-lambda11+lambda21))*(exp((mu11+mu21+c*(-lambda11+lambda21))*Llag/(2*c))-1)
    )
    c210L<-c*k11*k21*exp(-(mu11+mu21)*Llag/c)/(mu11+mu21)*(
      (exp(-lambda11*Llag)+exp(-lambda21*Llag))/(mu11+mu21+c*(lambda11+lambda21))+
        exp(-lambda11*Llag)/(mu11+mu21+c*(lambda11-lambda21))*(exp((mu11+mu21+c*(lambda11-lambda21))*Llag/(2*c))-1)+
        exp(-lambda21*Llag)/(mu11+mu21+c*(-lambda11+lambda21))*(exp((mu11+mu21+c*(-lambda11+lambda21))*Llag/(2*c))-1)
    )
    c220L<-c*k21^2/(2*mu21^2)*exp(-2*mu21*Llag/c-lambda21*Llag)*(exp(mu21*Llag/c)-c*lambda21/(mu21+c*lambda21))+
      c*k22^2/(2*mu22^2)*exp(-2*mu22*Llag/c-lambda22*Llag)*(exp(mu22*Llag/c)-c*lambda22/(mu22+c*lambda22))
    output<-c(output,c110L,c120L,c210L,c220L)
  }
  return(output)
}
datamoments<-function(Tlags,Llags,Y){
  Lx<-Y[[3]]
  Lt<-Y[[4]]
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  Cijhat<-function(i,j,Tlag,Llag){#i,j=1,2
    yi<-Y[[i]]
    yj<-Y[[j]]
    tpos<-(Tlag*(Nt-1))/Lt+1
    xpos<-(Llag*(Nx-1))/Lx+1
    return(cov(as.vector(yi[1:(Nt+1-tpos),1:(Nx+1-xpos)]),as.vector(yj[(tpos):Nt,(xpos):Nx])))
  }
  
  output<-c(Cijhat(1,1,0,0),Cijhat(1,2,0,0),Cijhat(2,2,0,0))
  for(Tlag in Tlags){
    output<-c(output,Cijhat(1,1,Tlag,0),Cijhat(1,2,Tlag,0),Cijhat(2,1,Tlag,0),Cijhat(2,2,Tlag,0))
  }
  for(Llag in Llags){
    output<-c(output,Cijhat(1,1,0,Llag),Cijhat(1,2,0,Llag),Cijhat(2,1,0,Llag),Cijhat(2,2,0,Llag))
  }
  return(output)
}

datamoments(Tlags,Llags,Y)
moments(Tlags,Llags,as.numeric(trueparas))
errors<-function(Tlags,Llags,Y,parameters){return((datamoments(Tlags,Llags,Y)-moments(Tlags,Llags,parameters)))}

initialvalue<-
gaussNewton(rep(0,9),function(parameters) errors(Tlags,Llags,Y,parameters))



