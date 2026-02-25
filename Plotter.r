setwd("/Users/willallfrey/Documents/R/ExpSepAmbit")
library(lattice)
library(gridExtra)

oonefine1<-readRDS("Output Fields/highrespointone.rds1.rds")
oonefine2<-readRDS("Output Fields/highrespointone.rds2.rds")
oonefine3<-readRDS("Output Fields/highrespointone.rds3.rds")
oonefine4<-readRDS("Output Fields/highrespointone.rds4.rds")
oonefine5<-readRDS("Output Fields/highrespointone.rds5.rds")

Y<-readRDS("Output Fields/highrespointone.rds2.rds")

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

plotfield(Y)

trajectoryOfY1(oonefine1)
trajectoryOfY1Y2(oonefine1)

library(ggplot2)
library(gganimate)
library(reshape2)

Y<-oonefine1
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
anim_save("timeplot.gif")

