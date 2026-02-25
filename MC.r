#MC

for(i in 1:40){
  Y<-readRDS(paste("Output Fields/newparasfine05v",i,".rds",sep=""))
  ggsave(filename=paste("Plots/errorconvergence05v",i,".png",sep=""),ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))
  print(ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))
  }

ggsave(filename=paste("Plots/errorconvergence05v",i,".png",sep=""),ploterrorconvergence(twostepparas(Y,10,0.001),as.numeric(Y[[7]])))

errorComparison<-matrix(0,40,3)
allparaests<-matrix(0,40,9)
for(i in 1:40){
  Y<-readRDS(paste("Output Fields/newparasfine05v",i,".rds",sep=""))
  #errorComparison[i,]<-c(errors1(Y),errors2(Y,5),errors2(Y,10),errors3(Y,5),errors3(Y,10),errors4(Y,5),errors4(Y,10))
  errorComparison[i,]<-c(errors1(Y),errors4(Y,5),errors4(Y,10))
  allparaests[i,]<-method4(Y,5)
}
colMeans(allparaests)
median(allparaests[,8])



i<-15
errors4(Y,10)

as.numeric(Y[[7]])
errorsData<-data.frame(errorComparison)
colnames(errorsData)<-c("m1","m4(k=5)","m4(k=10)")



method4(Y,10)
ploterrorconvergence(twostepparas(Y,4,0.001),Y[[7]])


mean(pmin(errorComparison[,2],errorComparison[,3]))
mean(errorComparison[,1])
