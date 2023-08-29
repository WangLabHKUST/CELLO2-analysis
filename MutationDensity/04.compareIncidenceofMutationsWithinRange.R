setwd("/Users/kevinmu/Documents/pangliomaevolution/MYC/HMnearMYC")
options(stringsAsFactors = F, scipen = 12)


#
prefix = "U87.H3K4me3.SRX129079.20"
#prefix = "U87.MYC.SRX129069.05"
prefix = "U87.H3K9me3.SRX1989885.05"
prefix="U87.RPII.SRX100504.20"
prefix="U87.DNase.SRX132055.05"
prefix="U87.MAX.SRX129085.20"
NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df)){#
  
  di = df$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df$ratioNHM = df$NHM/totNHM
df$ratioMYCHM = df$HMMYC/totHMMYC
df$rationonMYC = df$HMnoMYC/totHMnoMYC

plot(log10(df$distance), log10(df$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df$distance), log10(df$rationonMYC), type = 'l', col = 'orange')
points(log10(df$distance), log10(df$ratioMYCHM), type = 'l', col = 'red')



####
#compare  U87.H3K4me3MYC and U87.H3K4me3noMYC
####
#      ##
#    #   #   #
#   #     ###  #
# #             #
##               #
#
##################
prefix="U87.H3K4me3MYC"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df)){#
  
  di = df$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df$ratioNHM = df$NHM/totNHM
df$ratioMYCHM = df$HMMYC/totHMMYC
df$rationonMYC = df$HMnoMYC/totHMnoMYC

plot(log10(df$distance), log10(df$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df$distance), log10(df$rationonMYC), type = 'l', col = 'orange')
points(log10(df$distance), log10(df$ratioMYCHM), type = 'l', col = 'red')
df$MYCHMoverNHM = df$ratioMYCHM/df$ratioNHM
df$noMYCHMoverNHM = df$rationonMYC/df$ratioNHM

plot(log10(df$distance), log2(df$MYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Log2(HM/NHM)")
points(log10(df$distance), log2(df$noMYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'orange')


prefix="U87.H3k4me3noMYC"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df2 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df2)){#
  
  di = df2$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df2$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df2$ratioNHM = df2$NHM/totNHM
df2$ratioMYCHM = df2$HMMYC/totHMMYC
df2$rationonMYC = df2$HMnoMYC/totHMnoMYC

# plot(log10(df2$distance), log10(df2$ratioNHM), type = 'l', main = prefix,col = 'black',
#     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
# points(log10(df2$distance), log10(df2$rationonMYC), type = 'l', col = 'orange')
# points(log10(df2$distance), log10(df2$ratioMYCHM), type = 'l', col = 'red')
df2$MYCHMoverNHM = df2$ratioMYCHM/df2$ratioNHM
df2$noMYCHMoverNHM = df2$rationonMYC/df2$ratioNHM

points(log10(df2$distance), log2(df2$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red')
points(log10(df2$distance), log2(df2$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')


####
prefix="U87.MYCnoH3K4me3"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df3 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df3)){#
  
  di = df3$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df3$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df3$ratioNHM = df3$NHM/totNHM
df3$ratioMYCHM = df3$HMMYC/totHMMYC
df3$rationonMYC = df3$HMnoMYC/totHMnoMYC

plot(log10(df3$distance), log10(df3$rationonMYC), type = 'l', main = prefix,col = 'orange',
    xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log10(df3$ratioNHM), type = 'l', col = 'black')
points(log10(df3$distance), log10(df3$ratioMYCHM), type = 'l', col = 'red')
df3$MYCHMoverNHM = df3$ratioMYCHM/df3$ratioNHM
df3$noMYCHMoverNHM = df3$rationonMYC/df3$ratioNHM

plot(log10(df3$distance), log2(df3$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log2(df3$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')



###
#other TF
#MED1: U87.MED1noH3K4me3.bed; U87.H3k4me3noMED1.bed
#
prefix="U87.H3K4me3MED1"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df)){#
  
  di = df$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df$ratioNHM = df$NHM/totNHM
df$ratioMYCHM = df$HMMYC/totHMMYC
df$rationonMYC = df$HMnoMYC/totHMnoMYC

plot(log10(df$distance), log10(df$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df$distance), log10(df$rationonMYC), type = 'l', col = 'orange')
points(log10(df$distance), log10(df$ratioMYCHM), type = 'l', col = 'red')
df$MYCHMoverNHM = df$ratioMYCHM/df$ratioNHM
df$noMYCHMoverNHM = df$rationonMYC/df$ratioNHM

plot(log10(df$distance), log2(df$MYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Log2(HM/NHM)", ylim = c(0,1))
points(log10(df$distance), log2(df$noMYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'orange')


prefix="U87.H3k4me3noMED1"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df2 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df2)){#
  
  di = df2$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df2$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df2$ratioNHM = df2$NHM/totNHM
df2$ratioMYCHM = df2$HMMYC/totHMMYC
df2$rationonMYC = df2$HMnoMYC/totHMnoMYC

# plot(log10(df2$distance), log10(df2$ratioNHM), type = 'l', main = prefix,col = 'black',
#     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
# points(log10(df2$distance), log10(df2$rationonMYC), type = 'l', col = 'orange')
# points(log10(df2$distance), log10(df2$ratioMYCHM), type = 'l', col = 'red')
df2$MYCHMoverNHM = df2$ratioMYCHM/df2$ratioNHM
df2$noMYCHMoverNHM = df2$rationonMYC/df2$ratioNHM

points(log10(df2$distance), log2(df2$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red')
points(log10(df2$distance), log2(df2$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')


####
prefix="U87.MED1noH3K4me3"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df3 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df3)){#
  
  di = df3$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df3$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df3$ratioNHM = df3$NHM/totNHM
df3$ratioMYCHM = df3$HMMYC/totHMMYC
df3$rationonMYC = df3$HMnoMYC/totHMnoMYC

plot(log10(df3$distance), log10(df3$rationonMYC), type = 'l', main = prefix,col = 'orange',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log10(df3$ratioNHM), type = 'l', col = 'black')
points(log10(df3$distance), log10(df3$ratioMYCHM), type = 'l', col = 'red')
df3$MYCHMoverNHM = df3$ratioMYCHM/df3$ratioNHM
df3$noMYCHMoverNHM = df3$rationonMYC/df3$ratioNHM

plot(log10(df3$distance), log2(df3$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log2(df3$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')

###
#other TF
#BRD4: U87.H3K4me3BRD4.bed   U87.H3k4me3noBRD4.bed
#
prefix="U87.H3K4me3BRD4"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df)){#
  
  di = df$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df$ratioNHM = df$NHM/totNHM
df$ratioMYCHM = df$HMMYC/totHMMYC
df$rationonMYC = df$HMnoMYC/totHMnoMYC

plot(log10(df$distance), log10(df$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df$distance), log10(df$rationonMYC), type = 'l', col = 'orange')
points(log10(df$distance), log10(df$ratioMYCHM), type = 'l', col = 'red')
df$MYCHMoverNHM = df$ratioMYCHM/df$ratioNHM
df$noMYCHMoverNHM = df$rationonMYC/df$ratioNHM

plot(log10(df$distance), log2(df$MYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Log2(HM/NHM)")
points(log10(df$distance), log2(df$noMYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'orange')


prefix="U87.H3k4me3noBRD4"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df2 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df2)){#
  
  di = df2$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df2$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df2$ratioNHM = df2$NHM/totNHM
df2$ratioMYCHM = df2$HMMYC/totHMMYC
df2$rationonMYC = df2$HMnoMYC/totHMnoMYC

# plot(log10(df2$distance), log10(df2$ratioNHM), type = 'l', main = prefix,col = 'black',
#     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
# points(log10(df2$distance), log10(df2$rationonMYC), type = 'l', col = 'orange')
# points(log10(df2$distance), log10(df2$ratioMYCHM), type = 'l', col = 'red')
df2$MYCHMoverNHM = df2$ratioMYCHM/df2$ratioNHM
df2$noMYCHMoverNHM = df2$rationonMYC/df2$ratioNHM

points(log10(df2$distance), log2(df2$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red')
points(log10(df2$distance), log2(df2$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')


####
prefix="U87.BRD4noH3K4me3"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df3 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df3)){#
  
  di = df3$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df3$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df3$ratioNHM = df3$NHM/totNHM
df3$ratioMYCHM = df3$HMMYC/totHMMYC
df3$rationonMYC = df3$HMnoMYC/totHMnoMYC

plot(log10(df3$distance), log10(df3$rationonMYC), type = 'l', main = prefix,col = 'orange',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log10(df3$ratioNHM), type = 'l', col = 'black')
points(log10(df3$distance), log10(df3$ratioMYCHM), type = 'l', col = 'red')
df3$MYCHMoverNHM = df3$ratioMYCHM/df3$ratioNHM
df3$noMYCHMoverNHM = df3$rationonMYC/df3$ratioNHM

plot(log10(df3$distance), log2(df3$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log2(df3$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')


####
prefix="U87.H3K4me3BRD4noMYC"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df3 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df3)){#
  
  di = df3$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df3$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df3$ratioNHM = df3$NHM/totNHM
df3$ratioMYCHM = df3$HMMYC/totHMMYC
df3$rationonMYC = df3$HMnoMYC/totHMnoMYC

plot(log10(df3$distance), log10(df3$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log10(df3$rationonMYC), type = 'l', col = 'orange')
points(log10(df3$distance), log10(df3$ratioMYCHM), type = 'l', col = 'red')
df3$MYCHMoverNHM = df3$ratioMYCHM/df3$ratioNHM
df3$noMYCHMoverNHM = df3$rationonMYC/df3$ratioNHM

plot(log10(df3$distance), log2(df3$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)", main = prefix)
points(log10(df3$distance), log2(df3$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red')


###
#other TF
#REST: U87.H3K4me3REST.bed   U87.H3k4me3noREST.bed
#
prefix="U87.H3K4me3REST"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df)){#
  
  di = df$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df$ratioNHM = df$NHM/totNHM
df$ratioMYCHM = df$HMMYC/totHMMYC
df$rationonMYC = df$HMnoMYC/totHMnoMYC

plot(log10(df$distance), log10(df$ratioNHM), type = 'l', main = prefix,col = 'black',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df$distance), log10(df$rationonMYC), type = 'l', col = 'orange')
points(log10(df$distance), log10(df$ratioMYCHM), type = 'l', col = 'red')
df$MYCHMoverNHM = df$ratioMYCHM/df$ratioNHM
df$noMYCHMoverNHM = df$rationonMYC/df$ratioNHM

plot(log10(df$distance), log2(df$MYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Log2(HM/NHM)")
points(log10(df$distance), log2(df$noMYCHMoverNHM), type = 'l',lwd=2, main = prefix,col = 'orange')


prefix="U87.H3k4me3noREST"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df2 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df2)){#
  
  di = df2$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df2$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df2$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df2$ratioNHM = df2$NHM/totNHM
df2$ratioMYCHM = df2$HMMYC/totHMMYC
df2$rationonMYC = df2$HMnoMYC/totHMnoMYC

# plot(log10(df2$distance), log10(df2$ratioNHM), type = 'l', main = prefix,col = 'black',
#     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
# points(log10(df2$distance), log10(df2$rationonMYC), type = 'l', col = 'orange')
# points(log10(df2$distance), log10(df2$ratioMYCHM), type = 'l', col = 'red')
df2$MYCHMoverNHM = df2$ratioMYCHM/df2$ratioNHM
df2$noMYCHMoverNHM = df2$rationonMYC/df2$ratioNHM

points(log10(df2$distance), log2(df2$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red')
points(log10(df2$distance), log2(df2$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')


####
prefix="U87.RESTnoH3K4me3"

NPEAKS=nrow(read.delim(paste0(prefix,".bed"), header = F))
dts = c(1:99,seq(from=1,to=100, by=1)*100, seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000)
df3 = data.frame(distance =dts, NHM=0, HMMYC=0, HMnoMYC =0 )
for (i in 1:nrow(df3)){#
  
  di = df3$distance[i]
  print(di)
  fn = paste0("nonHM.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  
  df3$NHM[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn = paste0("HM.MYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
  fn =  paste0("HM.noMYCgain.mutations.chr.",prefix,".dist",di,".txt")
  f = read.delim(fn, header = F, stringsAsFactors = F)
  
  df3$HMnoMYC[i] = nrow(f)/NPEAKS#mean(table(f$V10)) #
  
}

totNHM = 14819/10000
totHMMYC = 26521/10000
totHMnoMYC = 69149/10000

df3$ratioNHM = df3$NHM/totNHM
df3$ratioMYCHM = df3$HMMYC/totHMMYC
df3$rationonMYC = df3$HMnoMYC/totHMnoMYC

plot(log10(df3$distance), log10(df3$rationonMYC), type = 'l', main = prefix,col = 'orange',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log10(df3$ratioNHM), type = 'l', col = 'black')
points(log10(df3$distance), log10(df3$ratioMYCHM), type = 'l', col = 'red')
df3$MYCHMoverNHM = df3$ratioMYCHM/df3$ratioNHM
df3$noMYCHMoverNHM = df3$rationonMYC/df3$ratioNHM

plot(log10(df3$distance), log2(df3$MYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'red',
     xlab = "Distance from peak (log10)", ylab = "Number of mutations (log10)")
points(log10(df3$distance), log2(df3$noMYCHMoverNHM), type = 'l', lty=2,lwd=2,col = 'orange')






