library(ggplot2)
library(reshape2)
library(lattice)
library(gridExtra)
library(gganimate)

Cijhat0<-function(i,j,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  return(mean(yi*yj))
}
CijhatTime<-function(i,j,Tlag,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  return(colMeans(yi[1:(Nt-Tlag/(2*deltat)),1:(Nx)]*yj[(Tlag/(2*deltat)+1):Nt,(1):Nx]))
}
CijhatSpace<-function(i,j,Llag,Y){#i,j=1,2
  yi<-Y[[i]]
  yj<-Y[[j]]
  deltat<-Y[[6]]
  c<-Y[[5]]
  Nt<-dim(yi)[1]
  Nx<-dim(yi)[2]
  return(rowMeans(yi[1:(Nt),1:(Nx-Llag/(2*deltat*c))]*yj[(1):Nt,(1+Llag/(2*deltat*c)):Nx]))
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

C11Time(0,Y[[7]],Y[[5]])==0.25
C22Time(0,Y[[7]],Y[[5]])==0.5
C12Time(0,Y[[7]],Y[[5]])==0.25

MCdata1<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
MCdata2<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
MCdata3<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())

for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/parapointfour05v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata1[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
  Y<-readRDS(paste("Output Fields/parapointfour05depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata2[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
  Y<-readRDS(paste("Output Fields/diffparas05depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata3[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
}
meltedMCdata1<-melt(MCdata1,id="iter")
meltedMCdata2<-melt(MCdata2,id="iter")
meltedMCdata3<-melt(MCdata3,id="iter")

ggplot(meltedMCdata1,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()
ggplot(meltedMCdata2,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()
ggplot(meltedMCdata3,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)#+geom_violin()

paraslist<-matrix(0,1000,9)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/parapointfour05v",i,".rds",sep=""))
  paraslist[i,]<-method4(Y,10)
}
trueparas<-Y[[7]]
allparadata<-data.frame(paraslist)
colnames(allparadata)=c("k11","k21","k22","mu11","mu21","mu22","lambda11","lambda21","lambda22")
meltedallparadata<-melt(allparadata)
kviolin<-ggplot(meltedallparadata,aes(x=variable,y=abs(value)))+geom_violin()+scale_x_discrete(limits=c("k11","k21","k22"))+#geom_boxplot(width=0.1)+
  annotate("point",x=1,y=trueparas$k11,col="blue")+
  annotate("point",x=2,y=trueparas$k21,col="blue")+
  annotate("point",x=3,y=trueparas$k22,col="blue")
muviolin<-ggplot(meltedallparadata,aes(x=variable,y=value))+geom_violin()+scale_x_discrete(limits=c("mu11","mu21","mu22"))+#geom_boxplot(width=0.1)+
  annotate("point",x=1,y=trueparas$mu11,col="blue")+
  annotate("point",x=2,y=trueparas$mu21,col="blue")+
  annotate("point",x=3,y=trueparas$mu22,col="blue")
lambdaviolin<-ggplot(meltedallparadata,aes(x=variable,y=value))+geom_violin()+scale_x_discrete(limits=c("lambda11","lambda21","lambda22"))+#geom_boxplot(width=0.1)+
  annotate("point",x=1,y=trueparas$lambda11,col="blue")+
  annotate("point",x=2,y=trueparas$lambda21,col="blue")+
  annotate("point",x=3,y=trueparas$lambda22,col="blue")
ggarrange(kviolin,muviolin,lambdaviolin,nrow=1)


K<-10
c<-Y[[5]]
deltat<-Y[[6]]
timelags<-((1:K))*2*deltat
timelags0<-c(0,timelags)
truec11time<-as.numeric(lapply(timelags0,function(tlag){C11Time(tlag,Y[[7]],Y[[5]])}))
truec12time<-as.numeric(lapply(timelags0,function(tlag){C12Time(tlag,Y[[7]],Y[[5]])}))
truec21time<-as.numeric(lapply(timelags0,function(tlag){C21Time(tlag,Y[[7]],Y[[5]])}))
truec22time<-as.numeric(lapply(timelags0,function(tlag){C22Time(tlag,Y[[7]],Y[[5]])}))
C11Timedata<-matrix(0,1000,K+1)
C12Timedata<-matrix(0,1000,K+1)
C21Timedata<-matrix(0,1000,K+1)
C22Timedata<-matrix(0,1000,K+1)
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",i,".rds",sep=""))
  for(j in 1:(K+1)){
    Tlag<-timelags0[j]
    C11Timedata[i,j]<-mean(CijhatTime(1,1,Tlag,Y))
    C12Timedata[i,j]<-mean(CijhatTime(1,2,Tlag,Y))
    C21Timedata[i,j]<-mean(CijhatTime(2,1,Tlag,Y))
    C22Timedata[i,j]<-mean(CijhatTime(2,2,Tlag,Y))
  }
}
c11timedataframe<-data.frame(C11Timedata)
c12timedataframe<-data.frame(C12Timedata)
c21timedataframe<-data.frame(C21Timedata)
c22timedataframe<-data.frame(C22Timedata)
colnames(c11timedataframe)<-timelags0
colnames(c12timedataframe)<-timelags0
colnames(c21timedataframe)<-timelags0
colnames(c22timedataframe)<-timelags0
meltedc11timedata<-melt(c11timedataframe)
meltedc12timedata<-melt(c12timedataframe)
meltedc21timedata<-melt(c21timedataframe)
meltedc22timedata<-melt(c22timedataframe)
c11tplt<-ggplot()+geom_boxplot(data=meltedc11timedata,aes(x=variable,y=value))+
  geom_point(data=NULL,aes(x=timelags0+1,y=truec11time),col="blue",shape="cross")+
  geom_line(data=NULL,aes(x=timelags0+1,y=truec11time),col="blue")+
  ggtitle("C_11^T(L) boxplots over 1000 realisations of field. Expected values in blue.")+
  xlab("L")+ylab("C_11 Time")
c12tplt<-ggplot()+geom_boxplot(data=meltedc12timedata,aes(x=variable,y=value))+
  geom_point(data=NULL,aes(x=timelags0+1,y=truec12time),col="blue",shape="cross")+
  geom_line(data=NULL,aes(x=timelags0+1,y=truec12time),col="blue")+
  ggtitle("C_12^T(L) boxplots over 1000 realisations of field. Expected values in blue.")+
  xlab("L")+ylab("C_12 Time")
c21tplt<-ggplot()+geom_boxplot(data=meltedc21timedata,aes(x=variable,y=value))+
  geom_point(data=NULL,aes(x=timelags0+1,y=truec21time),col="blue",shape="cross")+
  geom_line(data=NULL,aes(x=timelags0+1,y=truec21time),col="blue")+
  ggtitle("C_21^T(L) boxplots over 1000 realisations of field. Expected values in blue.")+
  xlab("L")+ylab("C_21 Time")
c22tplt<-ggplot()+geom_boxplot(data=meltedc22timedata,aes(x=variable,y=value))+
  geom_point(data=NULL,aes(x=timelags0+1,y=truec22time),col="blue",shape="cross")+
  geom_line(data=NULL,aes(x=timelags0+1,y=truec22time),col="blue")+
  ggtitle("C_22^T(L) boxplots over 1000 realisations of field. Expected values in blue.")+
  xlab("L")+ylab("C_22 Time")
ggarrange(c11tplt,c12tplt,c21tplt,c22tplt)

ggplot()+geom_boxplot(data=meltedc22timedata,aes(x=variable,y=value))+scale_x_discrete(limits=c("10"))+annotate("point",x=1,y=truec22time[11],col="blue",shape="cross")


MCdata4<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/ambit2!par04res05depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata4[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
}
meltedMCdata4<-melt(MCdata4,id="iter")
ggplot(meltedMCdata4,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()


MCdata5<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/ambit3par04res05depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata5[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
}
meltedMCdata5<-melt(MCdata5,id="iter")
ggplot(meltedMCdata5,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()

MCdata6<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/courseambit3par04res1depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata6[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
}
meltedMCdata6<-melt(MCdata6,id="iter")
ggplot(meltedMCdata6,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()

MCdata7<-data.frame(iter=character(),meanY1=character(),meanY2=character(),varY1=character(),varY2=character(),c12=character())
for(i in 1:1000){
  Y<-readRDS(paste("Output Fields/courseparas04res1depth3v",i,".rds",sep=""))
  m1<-mean(Y[[1]])
  m2<-mean(Y[[2]])
  var1centred<-(var(as.vector(Y[[1]]))-C11Time(0,Y[[7]],Y[[5]]))/C11Time(0,Y[[7]],Y[[5]])
  var2centred<-(var(as.vector(Y[[2]]))-C22Time(0,Y[[7]],Y[[5]]))/C22Time(0,Y[[7]],Y[[5]])
  C12centred<-(cov(as.vector(Y[[1]]),as.vector(Y[[2]]))-C12Time(0,Y[[7]],Y[[5]]))/C12Time(0,Y[[7]],Y[[5]])
  MCdata7[i,]<-c(i,m1,m2,var1centred,var2centred,C12centred)
}
meltedMCdata7<-melt(MCdata7,id="iter")
ggplot(meltedMCdata7,aes(x=variable,y=as.numeric(value)))+geom_boxplot(width=0.5)+ylim(-0.2,0.2)#+geom_violin()
