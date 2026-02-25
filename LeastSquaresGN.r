library(pracma)
oonefine1<-readRDS("Output Fields/highrespointone.rds1.rds")
oonefine2<-readRDS("Output Fields/highrespointone.rds2.rds")
oonefine3<-readRDS("Output Fields/highrespointone.rds3.rds")
oonefine4<-readRDS("Output Fields/highrespointone.rds4.rds")
oonefine5<-readRDS("Output Fields/highrespointone.rds5.rds")



moments<-function(Tlags,Llags,parameters,c){
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

errors<-function(Tlags,Llags,Y,parameters,datamoms){return(datamoms-moments(Tlags,Llags,parameters,Y[[5]]))}

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

getparameterestimate<-function(Y){
  c<-Y[[5]]
  deltat<-Y[[6]]
  deltax<-c*deltat
  Tlags<-1:30*(2*deltat)
  Llags<-1:30*(2*deltax)
  initialvalue<-as.numeric(paraestimate(Y,Tlags[1],Llags[1],c))[3:11]
  datamoms<-datamoments(Tlags,Llags,Y)
  lsparas<-gaussNewton(initialvalue,function(parameters) errors(Tlags,Llags,Y,parameters,datamoms))$xs
  return(lsparas)
}

paras<-getparameterestimate(Y)

paras1<-getparameterestimate(oonefine1)
paras2<-getparameterestimate(oonefine2)
paras3<-getparameterestimate(oonefine3)
paras4<-getparameterestimate(oonefine4)
paras5<-getparameterestimate(oonefine5)

rbind(as.numeric(oonefine1[[7]])-paras1,as.numeric(oonefine2[[7]])-paras2,as.numeric(oonefine3[[7]])-paras3,as.numeric(oonefine4[[7]])-paras4,as.numeric(oonefine5[[7]])-paras5)
sum(abs(rbind(as.numeric(oonefine1[[7]])-paras1,as.numeric(oonefine2[[7]])-paras2,as.numeric(oonefine3[[7]])-paras3,as.numeric(oonefine4[[7]])-paras4,as.numeric(oonefine5[[7]])-paras5)))
sum(abs(rbind(as.numeric(oonefine1[[7]])-as.numeric(paraestimate(oonefine1)[3:11]),as.numeric(oonefine2[[7]])-as.numeric(paraestimate(oonefine2)[3:11]),as.numeric(oonefine3[[7]])-as.numeric(paraestimate(oonefine3)[3:11]),as.numeric(oonefine4[[7]])-as.numeric(paraestimate(oonefine4)[3:11]),as.numeric(oonefine5[[7]])-as.numeric(paraestimate(oonefine5)[3:11]))))

