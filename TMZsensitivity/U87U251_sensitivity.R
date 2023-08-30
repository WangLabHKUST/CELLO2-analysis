#change working directory
library(rstudioapi)
present_folder = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(present_folder)
print(getwd())

library(drc)
inputfile<-"U87U251.drc.dat"
data0<-read.table(inputfile,sep="\t",header=T,quote=NULL, stringsAsFactors = F)

data = data0[startsWith(data0$Cellline,"U251"),]
modelC3 <- drm(Viability~Conc,Cellline, fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")), data=data)
pdf('Fig5D_sensitivity.pdf', width = 4, height=4)
plot(modelC3, ylim = c(0,1.),type = 'bars', col = c("#fdbf6f","#ff7f00"), lty=c(3,1), lwd = 2, legend = F, pch = c(20,17),cex=0.75,
     ylab="Viability", xlab="Temozolomide (µM)", broken = T, add= F)

legend('bottomleft',legend = c('U251','U251TR'), pch=c(17,20), col = c("#fdbf6f","#ff7f00"),lty=c(1,3),bty='n', cex = 0.65)
dev.off()

dt= data0[startsWith(data0$Cellline,"U87"),]#
dt = dt[dt$Conc<=1000,]
modelC3 <- drm(Viability~Conc,Cellline, fct=LL.4(names=c("Slope", "Lower", "Upper", "ED50")), data=dt)

plot(modelC3, ylim = c(0,1.),type = 'bars',
     col = c("#a6cee3","#1f78b4"), lty=c(3,1,1,3,1), lwd = 2, legend = F, pch = c(20,17,18,19,16),cex=0.75,
     ylab="Viability", xlab="Temozolomide (µM)", broken = T,legendPos = c(100,0.3))
legend('bottomleft',legend = c('U87TR','U87'), pch=c(17,20), col = c("#1f78b4","#a6cee3"),lty=c(1,3),bty='n', cex = 0.65)

