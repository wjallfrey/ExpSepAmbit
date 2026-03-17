#MC

for(i in 1:50){
  Y<-readRDS(paste("Output Fields/newparasfine05v",i,".rds",sep=""))
  ggsave(filename=paste("Plots/errorconvergence05v",i,".png",sep=""),ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))
  print(ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))
  }

ggsave(filename=paste("Plots/errorconvergence05v",i,".png",sep=""),ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))

errorComparison<-matrix(0,40,3)
allparaests<-matrix(0,50,9)

for(i in 1:50){
  Y<-readRDS(paste("Output Fields/newparasfine05v",i,".rds",sep=""))
  #errorComparison[i,]<-c(errors1(Y),errors2(Y,5),errors2(Y,10),errors3(Y,5),errors3(Y,10),errors4(Y,5),errors4(Y,10))
  #errorComparison[i,]<-c(errors1(Y),errors4(Y,5),errors4(Y,10))
  allparaests[i,]<-method4(Y,5)
  #print(c(mean(Y[[1]]),mean(Y[[2]])))
}

for(i in 1:50){
  Y<-readRDS(paste("Output Fields/newnewparas05v",i,".rds",sep=""))
  allparaests[i,]<-method4(Y,5)
}

allparaests<-matrix(0,50,9)
errorComparison<-matrix(0,50,9)
listofcovariance<-matrix(0,50,3)
for(i in 1:50){
  Y<-readRDS(paste("Output Fields/bulk01v",i,".rds",sep=""))
  #lis2[i]<-(mean(normalisedgammaijTimeHAT(1,1,0.2,Y)))
  #allparaests[i,]<-method4(Y,10)
  #allparaests[i,]<-as.numeric(getparametersOAAT(Y,10))
  listofcovariance[i,1]<-Cijhat0(1,1,Y)
  listofcovariance[i,2]<-Cijhat0(1,2,Y)
  listofcovariance[i,3]<-Cijhat0(2,2,Y)
  #errorComparison[i,]=c(errors1(Y),errors4(Y,2),errors4(Y,4),errors4(Y,6),errors4(Y,8),errors4(Y,10),errors4(Y,12),errors4(Y,14),errors4(Y,16))
}


trueparas<-Y[[7]]
allparadata<-data.frame(allparaests)
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

violinplot<-function(paralist,trueparas){
  allparadata<-data.frame(paralist)
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
}


da<-data.frame(t((t(listofcovariance)-c(C11Time(0,Y[[7]],Y[[5]]),C12Time(0,Y[[7]],Y[[5]]),C22Time(0,Y[[7]],Y[[5]])))/c(C11Time(0,Y[[7]],Y[[5]]),C12Time(0,Y[[7]],Y[[5]]),C22Time(0,Y[[7]],Y[[5]]))))
colnames(da)<-c("c11","c12","c22")
meltedcovdata<-melt(da)
ggplot(meltedcovdata,aes(x=variable,y=value))+geom_boxplot()+ggtitle("Normalised and centred estimation of the (cross-)covariances")


ggplot(data.frame(lis2))+geom_boxplot(aes(x=lis2))+annotate("point",x=C12Time(0,Y[[7]],Y[[5]]),y=0,col="blue")



colMeans(allparaests)


errors<-t(allparaests)-as.numeric(Y[[7]])

i<-10
Y<-readRDS(paste("Output Fields/bulk01v",i,".rds",sep=""))
twostepparas(Y,5,0.001)
plotvariograms(twostepparas(Y,20,0.001,T))
plotfield(Y)

as.numeric(Y[[7]])
errorsData<-data.frame(errorComparison)
colnames(errorsData)<-c("m1","m4(k=5)","m4(k=10)")



method4(Y,10)
ploterrorconvergence(twostepparas(Y,4,0.001),Y[[7]])


mean(pmin(errorComparison[,2],errorComparison[,3]))
mean(errorComparison[,1])
