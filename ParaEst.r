setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")
library(lattice)
library(lamW)
library(gridExtra)

#Plotting
plotfun<-function(Y){
  y1<-Y[[1]]
  y2<-Y[[2]]
  Nx<-dim(y1)[2]
  Nt<-dim(y1)[1]
  Lx<-Y[[3]]
  Lt<-Y[[4]]
  Xdomain<-(((-(Nx-1)/2)):((Nx-1)/2))/((Nx-1)/Lx)
  Tdomain<-0:(Nt-1)*(Lt/(Nt-1))
  grid<-expand.grid(t = Tdomain,x = Xdomain)
  grid$y1vec<-as.vector(y1)
  grid$y2vec<-as.vector(y2)
  cdomscale1<-ceiling(max(max(y1),-min(y1))*10)/10
  cdomscale2<-ceiling(max(max(y2),-min(y2))*10)/10
  plty1<-levelplot(y1vec~x*t,grid,main="Y1",at=c((-35:35)/35*cdomscale1))
  plty2<-levelplot(y2vec~x*t,grid,main="Y2",at=c((-35:35)/35*cdomscale2))
  grid.arrange(plty1,plty2)
}

#Parameter estimation
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
  lambda21hat<-uniroot(torootl21,c(0,10))$root #Not sure how to choose the bounds
  
  #Estimation of k11k12
  k11k21hat<-Cijhat(1,2,0,0)*(mu11hat+mu21hat)*(mu11hat+mu21hat+c*(lambda11hat+lambda21hat))/(2*c)
  
  #Estimation of mu22
  k21hat2<-k11k21hat^2/k11hat2
  mu22hat<- -1/Tlag*log((2*mu21hat*(mu21hat+c*lambda21hat)*Cijhat(2,2,Tlag,0)-c*k21hat2*exp(-mu21hat*Tlag))/(2*mu21hat*(mu21hat+c*lambda21hat)*Cijhat(2,2,0,0)-c*k21hat2))
  
  #Estimation of lambda22
  torootl22<-function(l22){
    (Cijhat(2,2,0,0)-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*exp(-2*mu22hat*Llag/c)*(c/(Llag*mu22hat))*(exp(mu22hat*Llag/c)-1)*exp(-l22*Llag)*(l22*Llag+mu22hat*Llag/c*(exp(mu22hat*Llag/c)/(exp(mu22hat*Llag/c)-1)))-Cijhat(2,2,0,Llag)+(c*k21hat2)/(2*mu21hat^2)*exp(-(2*mu21hat/c+lambda21hat)*Llag)*(exp(mu21hat*Llag/c)-c*lambda21hat/(mu21hat+c*lambda21hat))
  }
  lambda22hat<-uniroot(torootl22,c(0,10))$root
  
  #estimaiton of k22
  k22hat2<-(Cijhat(2,2,0,0)-(c*k21hat2)/(2*mu21hat*(mu21hat+c*lambda21hat)))*2*mu22hat*(mu22hat+c*lambda22hat)/c
  
  paraests<-list(Tlag=Tlag,Llag=Llag,k11est=(k11hat2)^0.5,k21est=k11k21hat/(k11hat2)^0.5,k22est=k22hat2^0.5,mu11est=mu11hat,mu21est=mu21hat,mu22est=mu22hat,lambda11est=lambda11hat,lambda21est=lambda21hat,lambda22est=lambda22hat)
  return(paraests)
}


Y<-readRDS("Yfield20.rds")
plotfun(Y)
parametertable<-rbind(NULL,list(Tlag=NULL,Llag=NULL,k11est=1,k21est=0.5,k22est=2,mu11est=0.1,mu21est=0.4,mu22est=0.2,lambda11est=0.3,lambda21est=0.6,lambda22est=0.5))
deltat<-Y[[6]]
c<-Y[[5]]
deltax<-c*deltat
for(i in 1:5){
  for(j in 1:5){
    parametertable<-rbind(parametertable,paraestimate(Y,Tlag=2*deltat*i,Llag=2*deltax*j))
  }
}
paraframe<-data.frame(parametertable)

par(mfrow=c(3,1))
plot(as.numeric(paraframe$k11est),type="l",ylab = "k estimates",ylim=c(0,3))
lines(as.numeric(paraframe$k21est))
lines(as.numeric(paraframe$k22est))
plot(as.numeric(paraframe$mu11est),type="l",ylab = "mu estimates",ylim=c(0,0.5))
lines(as.numeric(paraframe$mu21est))
lines(as.numeric(paraframe$mu22est))
plot(as.numeric(paraframe$lambda11est),type="l",ylab = "lambda estimates",ylim=c(0,5))
lines(as.numeric(paraframe$lambda21est))
lines(as.numeric(paraframe$lambda22est))

y1<-Y[[1]]
y2<-Y[[2]]


cov(as.vector(y2[1:(Nt-1),1:100]),as.vector(y1[2:Nt,1:100]))
sum<-0
for(i in 1:100){
  sum<-sum+cov(as.vector(y2[1:(Nt-1),i]),as.vector(y1[2:Nt,i]))
}
print(sum/100)

theoreticalmoments<-function(Tlag,Llag,k11,k21,k22,mu11,mu21,mu22,lambda11,lambda21,lambda22){
  c1100<-c*k11^2/(2*mu11*(mu11+c*lambda11))
  c11T0<-c*k11^2/(2*mu11*(mu11+c*lambda11))*exp(-Tlag*mu11)
  c110L<-c*k11^2/(2*mu11^2)*exp(-2*mu11*Llag/c-lambda11*Llag)*(exp(mu11*Llag/c)-c*lambda11/(mu11+c*lambda11))
  c1200<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))
  c12T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu21*Tlag)
  c120L<-0
  c21T0<-2*c*k11*k21/((mu11+mu21)*(mu11+mu21+c*lambda11+c*lambda21))*exp(-mu11*Tlag)
  c210L<-0
  c2200<-c*k21^2/(2*mu21*(mu21+c*lambda21))+c*k22^2/(2*mu22*(mu22+c*lambda22))
  c22T0<-c*k21^2/(2*mu21*(mu21+c*lambda21))*exp(-mu21*Tlag)+c*k22^2/(2*mu22*(mu22+c*lambda22))*exp(-mu22*Tlag)
  c220L<-0
  return(list(c1100=c1100,c11T0=c11T0,c110L=c110L,c1200=c1200,c12T0=c12T0,c120L=c120L,c21T0=c21T0,c210L=c210L,c2200=c2200,c22T0=c22T0,c220L=c220L))
}
actualmoments<-function(Tlag,Llag,Y){
  Cijhat<-function(i,j,Tlag,Llag){#i,j=1,2
    yi<-Y[[i]]
    yj<-Y[[j]]
    tpos<-(Tlag*(Nt-1))/Lt+1
    xpos<-(Llag*(Nx-1))/Lx+1
    return(cov(as.vector(yi[1:(Nt+1-tpos),1:(Nx+1-xpos)]),as.vector(yj[(tpos):Nt,(xpos):Nx])))
  }
  return(list(C1100=Cijhat(1,1,0,0),C11T0=Cijhat(1,1,Tlag,0),C110L=Cijhat(1,1,0,Llag),C1200=Cijhat(1,2,0,0),C12T0=Cijhat(1,2,Tlag,0),C120L=Cijhat(1,2,0,Llag),C21T0=Cijhat(2,1,Tlag,0),C210L=Cijhat(2,1,0,Llag),C2200=Cijhat(2,2,0,0),C22T0=Cijhat(2,2,Tlag,0),C220L=Cijhat(2,2,0,Llag)))
}

actualmoments2<-function(Tlag,Llag,Y){
  Cijhat<-function(i,j,Tlag,Llag){#i,j=1,2
    yi<-Y[[i]]
    yj<-Y[[j]]
    if(Llag==0){
      tpos<-(Tlag*(Nt-1))/Lt+1
      covs<-rep(0,Nx)
      for(j in 1:(Nx)){
        covs[j]<-cov(as.vector(yi[1:(Nt+1-tpos),j]),as.vector(yj[(tpos):Nt,j]))
      }
      return(sum(covs)/Nx)
    }
    else if(Tlag==0){
      xpos<-(Llag*(Nx-1))/Lx+1
      covs<-rep(0,Nt)
      for(i in 1:(Nt)){
        covs[i]<-cov(as.vector(yi[i,1:(Nx+1-xpos)]),as.vector(yj[i,(xpos):Nx]))
      }
      return(sum(covs)/Nt)
    }
  }
  return(list(C1100=Cijhat(1,1,0,0),C11T0=Cijhat(1,1,Tlag,0),C110L=Cijhat(1,1,0,Llag),C1200=Cijhat(1,2,0,0),C12T0=Cijhat(1,2,Tlag,0),C120L=Cijhat(1,2,0,Llag),C21T0=Cijhat(2,1,Tlag,0),C210L=Cijhat(2,1,0,Llag),C2200=Cijhat(2,2,0,0),C22T0=Cijhat(2,2,Tlag,0),C220L=Cijhat(2,2,0,Llag)))
}



