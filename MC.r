#MC

for(i in 1:10){
  Y<-readRDS(paste("Output Fields/newparasfine05v",i,".rds",sep=""))
  png(ploterrorconvergence(Y,10,0.001))
}

ggsave(filename="Plots/errorconvergence.png",ploterrorconvergence(Y,10,0.001))

