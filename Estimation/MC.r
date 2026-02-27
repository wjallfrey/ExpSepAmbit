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
for(i in 41:50){
  Y<-readRDS(paste("Output Fields/bulk01v",i,".rds",sep=""))
  allparaests[i,]<-method4(Y,15)
  errorComparison[i,]=c(errors1(Y),errors4(Y,2),errors4(Y,4),errors4(Y,6),errors4(Y,8),errors4(Y,10),errors4(Y,12),errors4(Y,14),errors4(Y,16))
}


trueparas<-Y[[7]]
allparadata<-data.frame(allparaests)
colnames(allparadata)=c("k11","k21","k22","mu11","mu21","mu22","lambda11","lambda21","lambda22")
meltedallparadata<-melt(allparadata)
kviolin<-ggplot(meltedallparadata,aes(x=variable,y=value))+geom_violin()+scale_x_discrete(limits=c("k11","k21","k22"))+#geom_boxplot(width=0.1)+
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


colMeans(allparaests)




i<-17
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
