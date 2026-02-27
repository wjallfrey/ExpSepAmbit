setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")
library(lattice)
library(gridExtra)
library(ggplot2)
library(gganimate)
library(reshape2)

Y<-readRDS("Output Fields/bulk01v17.rds")

plotfield<-function(Y){
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

trajectoryOfY1<-function(Y){
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  plot(Y[[1]][,1],type='l',col="black",main="Y1 at x=-Lx/2,x=0,x=Lx/2 in black red blue")
  lines(Y[[1]][,(Nx-1)/2],col="red")
  lines(Y[[1]][,Nx],col="blue")
}
trajectoryOfY1Y2<-function(Y){
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  plot(Y[[1]][,(Nx-1)/2],type='l',col="black",main="Y1 in black and Y2 in blue over time",ylim=c(min(Y[[1]][,(Nx-1)/2],Y[[2]][,(Nx-1)/2]),max(Y[[1]][,(Nx-1)/2],Y[[2]][,(Nx-1)/2])))
  lines(Y[[2]][,(Nx-1)/2],col="blue")
}
gifintime<-function(Y,nameforfile){
  if(!missing(nameforfile)){if(typeof(nameforfile)!="character"){
    return("Need string for file name")
  }}
  Nx<-dim(Y[[1]])[2]
  Nt<-dim(Y[[1]])[1]
  Lx<-Y[[3]]
  Lt<-Y[[4]]
  Y1data<-data.frame(t=0:(Nt-1)*Lt/(Nt-1),Y[[1]])
  Y2data<-data.frame(t=0:(Nt-1)*Lt/(Nt-1),Y[[2]])
  dfnew<-data.frame(melt(Y1data,id="t"),Y2=melt(Y2data,id="t")[,3])
  Ydata<-dfnew
  Ydata[,2]<-(as.numeric(dfnew[,2])-1)/(Nx-1)*Lx-Lx/2
  colnames(Ydata)<-c("t","x","Y1","Y2")
  animation<-ggplot(Ydata)+geom_line(aes(y=Y1,x=x,group=t))+geom_line(aes(y=Y2,x=x,group=t),col="blue")+labs(title="Y at time t = {round(frame_time)}")+xlab("x")+ylab("Y(t,x)")+transition_time(t)
  animate(animation,nframes=501,width=960,height=480)
  if(!missing(nameforfile)){
    anim_save(paste(nameforfile,".gif",sep=""))
  }
}

plotvariograms<-function(fulldata){
  paralist<-fulldata$paralist
  trueparas<-fulldata$trueparas
  datavars<-fulldata$datavars
  timelags<-fulldata$timelags
  spacelags<-fulldata$spacelags
  c<-fulldata$c
  NumberofLags<-length(timelags)
  estimatedparas<-aslist(as.numeric(tail(paralist,1)))
  fullvariogramdata<-data.frame(timelags,spacelags,matrix(c(rowMeans(datavars$timedata),rowMeans(datavars$spacedata),theoreticalvariogramsTime(estimatedparas,timelags,c),theoreticalvariogramsSpace(estimatedparas,spacelags,c),theoreticalvariogramsTime(trueparas,timelags,c),theoreticalvariogramsSpace(trueparas,spacelags,c)),NumberofLags,18))
  colnames(fullvariogramdata)=c("timelags","spacelags","g11hatt","g12hatt","g22hatt","g11hats","g12hats","g22hats","g11t","g12t","g22t","g11s","g12s","g22s","trueg11t","trueg12t","trueg22t","trueg11s","trueg12s","trueg22s")
  p11t<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g11t))+geom_point(aes(timelags,g11hatt))+geom_line(aes(timelags,trueg11t),col="blue")+xlab("d_t")+ylab("gamma11T(d_t)")
  p12t<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g12t))+geom_point(aes(timelags,g12hatt))+geom_line(aes(timelags,trueg12t),col="blue")+xlab("d_t")+ylab("gamma12T(d_t)")
  p22t<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g22t))+geom_point(aes(timelags,g22hatt))+geom_line(aes(timelags,trueg22t),col="blue")+xlab("d_t")+ylab("gamma22T(d_t)")
  p11s<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g11s))+geom_point(aes(timelags,g11hats))+geom_line(aes(timelags,trueg11s),col="blue")+xlab("d_t")+ylab("gamma11S(d_t)")
  p12s<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g12s))+geom_point(aes(timelags,g12hats))+geom_line(aes(timelags,trueg12s),col="blue")+xlab("d_t")+ylab("gamma12S(d_t)")
  p22s<-ggplot(fullvariogramdata)+geom_line(aes(timelags,g22s))+geom_point(aes(timelags,g22hats))+geom_line(aes(timelags,trueg22s),col="blue")+xlab("d_t")+ylab("gamma22S(d_t)")
  ggarrange(p11t,p12t,p12t,p22t,p11s,p12s,p12s,p22s,nrow=4,ncol=2)
}

fulldata<-twostepparas(Y,10,0.001,T)
plotvariograms(fulldata)

plot(rowMeans((datavars$timedata-rowMeans(datavars$timedata))^2))
plot(rowMeans((datavars$timedata-rowMeans(datavars$timedata))^2)/rowMeans(datavars$timedata))
plot(rowMeans((datavars$timedata-theoreticalvariogramsTime(aslist(as.numeric(tail(paralist,1))),timelags,c))^2)/theoreticalvariogramsTime(aslist(as.numeric(tail(paralist,1))),timelags,c)^2)
plot(rowMeans((datavars$spacedata-theoreticalvariogramsSpace(aslist(as.numeric(tail(paralist,1))),timelags,c))^2)/theoreticalvariogramsSpace(aslist(as.numeric(tail(paralist,1))),spacelags,c)^2)
plot(rowMeans((datavars$spacedata-rowMeans(datavars$spacedata))^2)/rowMeans(datavars$spacedata)^2)


