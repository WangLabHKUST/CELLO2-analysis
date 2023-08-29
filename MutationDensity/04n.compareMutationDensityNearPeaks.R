setwd("~/Dropbox/communter/HMWGS/")
options(stringsAsFactors = F, scipen = 12)
seeds = c(123,42,5,666,7890,1,2,22,3,333)
library(data.table)
library(reshape2)
library(ggplot2)
library(scales)
#N = 150000
#sample 150000 mutations from each group
#repeat 5 times
#ignore if already done
notdone = FALSE
if (notdone){
  s0 = read.delim("PairedGlioma.WGS.mutations.HM.ix.bed", header = F)
  s1 = read.delim("PairedGlioma.WGS.mutations.HMMYC.ix.bed", header = F)
  s2 = read.delim("PairedGlioma.WGS.mutations.HMnoMYC.ix.bed", header = F)
  s3 = read.delim("PairedGlioma.WGS.mutations.NHM.ix.bed", header = F)
  
  
  for (sdi in seeds){
    set.seed(sdi)
    
    N0 = 150000
    s0i = s0[sample(1:nrow(s0),N0),]
    f0i = paste0("PairedGlioma.WGS.mutations.HM.ix.sample150000.seed",sdi,".bed")
    write.table(s0i, file = f0i,row.names = F, col.names = F, quote = F, sep = "\t" )
    
    N1 = 150000#nrow(s1)/6 #6 is the number of HM with MYC gain
    s1i = s1[sample(1:nrow(s1),N1),]
    f1i = paste0("PairedGlioma.WGS.mutations.HMMYC.ix.sample150000.seed",sdi,".bed")
    write.table(s1i, file = f1i,row.names = F, col.names = F, quote = F, sep = "\t" )
    
    N2 = 150000#nrow(s2)/8 #8 is the number of HM with no MYC gain
    s2i = s2[sample(1:nrow(s2),N2),]
    f2i = paste0("PairedGlioma.WGS.mutations.HMnoMYC.ix.sample150000.seed",sdi,".bed")
    write.table(s2i, file = f2i,row.names = F, col.names = F, quote = F, sep = "\t" )
    
    N3 = 150000#nrow(s3)/64 #64 is the number of NHM patients
    s3i = s3[sample(1:nrow(s3),N3),]
    f3i = paste0("PairedGlioma.WGS.mutations.NHM.ix.sample150000.seed",sdi,".bed")
    write.table(s3i, file = f3i,row.names = F, col.names = F, quote = F, sep = "\t" )
  }
}
#\\sample 150000 mutations from each group
#\\repeat 5 times
#\\ignore if already done
#
#prefix="U87.MYC.SRX129069.05"
#prefix="U87.MAX.SRX129085.20"
prefix="U87.H3K9me3.SRX1989885.05"
prefix = "U87.DNase.SRX132056.05"

prefix="U87.H3K4me3.SRX129079.20"
prefix = "U87.H3K27ac.SRX129073.05"
prefix = "U87.RPII.SRX100504.20"
prefix="ENCFF870JRW" #h3k36me3

#prefix = "ENCFF870JRW.activemediangt1"
#prefix = "ENCFF870JRW.activemedianle1"

#prefix='H3K4me3MYCDNase'
# prefix='H3K4me3MAXDNase'
# prefix='RPII_and_MYC5k_DNase'
# prefix="U87.H3K27ac.SRX129073.05"
# 
# prefix="ENCFF870JRW" #glial cell line, H3K36me3, https://www.encodeproject.org/experiments/ENCSR136NUH/
# prefix = "ENCFF687HQC"#glial cell line, H3K4me3, https://www.encodeproject.org/experiments/ENCSR433PUR/
# prefix="ENCFF059VOH" #glial cell line, H3K9me3, https://www.encodeproject.org/experiments/ENCSR104BSN/
# prefix = "ENCFF856JTO" #glial cell line H3K4me1, 1
# prefix="test.H3K36me3MYC5k"
# prefix="test.H3K36me35kMYC"
prefix="U87.DNaseH3K27ac.withMYCwa"
prefix="U87.DNaseH3K27ac.withoutMYCA"

prefix="U87.DNaseH3K4me3.withMYCwa"
prefix="U87.DNaseH3K4me3.withoutMYCA"

prefix = "U87.RPII_and_MYC1k"
prefix = "U87.RPII_no_MYC1k"

prefix="ENCFF870JRW5k.MYCbind"
prefix= "ENCFF870JRW5k.nonMYCbind"

prefix = "U87.DNase.MYCbind"
prefix="U87.DNase.nonMYCbind"
#--
df = read.delim(paste0("MutationDensity.per150000.",prefix,".txt"))
#df = df[df$bin>=100,]
df$Density = df$Muts*1000000/df$Width
dts = c(1:9,seq(from=1,to=9, by=1)*10,seq(from=1,to=9, by=1)*100, seq(from=1,to=10, by=1)*1000,seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000,5000000,10000000)
df2 = data.frame(distance =dts, 
                 HM_mn=0,HM_lb = 0, HM_ub = 0,
                 NHM_mn=0,NHM_lb = 0, NHM_ub = 0,
                 HMMYC_mn =0, HMMYC_lb =0,HMMYC_ub =0,
                 HMnoMYC_mn =0, HMnoMYC_lb =0, HMnoMYC_ub =0)
for (i in 1:nrow(df2)){
  tmpi = df[df$bin==df2$distance[i],]
  
  d0 = tmpi$Density[tmpi$Type=="HM"]
  df2$HM_mn[i] = mean(d0)
  df2$HM_lb[i] = mean(d0)-sd(d0)
  df2$HM_ub[i] = mean(d0)+sd(d0)
  
  d1 = tmpi$Density[tmpi$Type=="NHM"]
  df2$NHM_mn[i] = mean(d1)
  df2$NHM_lb[i] = mean(d1)-sd(d1)
  df2$NHM_ub[i] = mean(d1)+sd(d1)
  
  d2 = tmpi$Density[tmpi$Type=="HMMYC"]
  df2$HMMYC_mn[i] = mean(d2)
  df2$HMMYC_lb[i] = mean(d2) - sd(d2)
  df2$HMMYC_ub[i] = mean(d2) + sd(d2)
  
  d3 = tmpi$Density[tmpi$Type=="HMnoMYC"]
  df2$HMnoMYC_mn[i] = mean(d3)
  df2$HMnoMYC_lb[i] = mean(d3) - sd(d3)
  df2$HMnoMYC_ub[i] = mean(d3) + sd(d3)
}

##----compare HM versus NHM

df2_mean = reshape2::melt(df2[,c('distance','NHM_mn','HM_mn')],id.vars = 'distance')
df2_ci = data.frame(distance = rep(df2$distance,2), mutations = rep(c('NHM','HM'), each = nrow(df2)),
                    lo = c(df2$NHM_lb, df2$HM_lb), hi = c(df2$NHM_ub, df2$HM_ub))

prefixn = prefix #"Glia_H3K4me1"
scaleFactor = 150000/(2855879637/1000000)
p <- ggplot()+
  geom_line(aes(x =distance, y = value/scaleFactor, color = variable), data = df2_mean,lwd=1)+
  scale_color_manual(values = c("#1f78b4","#e31a1c"))+
  geom_ribbon(aes(x = distance,ymin = lo/scaleFactor, ymax= hi/scaleFactor, group = mutations, fill = mutations), data = df2_ci, alpha=0.7)+
  scale_fill_manual(values = c("#e31a1c","#1f78b4"))+
  #geom_errorbar(aes(x = distance,ymin=lo, ymax=hi,group = mutations, fill = mutations),data = df2_ci, width=.1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),#limits=c(10^2,2*10^6),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #             labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 1,lty=2)+
  #lims(y=c(0.5,2))+
  #expand_limits(x = c(1,NA))+
  theme_classic()+#annotation_logticks()+
  theme(axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  labs(x = "Distance from the ChIP peaks (bp)", color = "",y = "Normalized mutation density", title = prefixn)

p
formatC(wilcox.test(value ~ variable, data = df2_mean)$p.value, format = "e", digits = 1)
#ggsave(filename = paste0("MutationDensity.HMvsNHMper150000.",prefixn,".pdf"),plot = p,width = 3.05,height = 1.6,units = 'in',dpi = 300)
##---- separate MYCgainHM and nonMYCgainHM
df2_mean = reshape2::melt(df2[,c('distance','NHM_mn','HMMYC_mn','HMnoMYC_mn')],id.vars = 'distance')
df2_ci = data.frame(distance = rep(df2$distance,3), mutations = rep(c('NHM','HMMYC','HMnoMYC'), each = nrow(df2)),
                    lo = c(df2$NHM_lb, df2$HMMYC_lb, df2$HMnoMYC_lb), hi = c(df2$NHM_ub, df2$HMMYC_ub, df2$HMnoMYC_ub))

prefixn = prefix #"Glia_H3K4me1"
p <- ggplot()+
  geom_line(aes(x =distance, y = value/scaleFactor, color = variable), data = df2_mean,lwd=1)+
  scale_color_manual(values = c("#1f78b4","red","orange"))+
  geom_ribbon(aes(x = distance,ymin = lo/scaleFactor, ymax= hi/scaleFactor, group = mutations, fill = mutations), data = df2_ci, alpha=0.7)+
  scale_fill_manual(values = c("red","orange","#1f78b4"))+
  #geom_errorbar(aes(x = distance,ymin=lo, ymax=hi,group = mutations, fill = mutations),data = df2_ci, width=.1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),#limits=c(10^2,2*10^6),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #             labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 1,lty=2)+
  #lims(y=c(.5,2.5))+
  expand_limits(x = 1)+
  theme_classic()+#annotation_logticks()+
  theme(axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  labs(x = "Distance from the ChIP peaks (bp)", color = "",y = "Normalized mutation density", title = prefixn)

p

#ggsave(filename = paste0("MutationDensity.per150000.",prefixn,".pdf"),plot = p,width = 3.5,height = 2.5,units = 'in',dpi = 300)

#----
#---
# SF = 150000/(2855879637/1000000)
# q <- ggplot()+
#   geom_line(aes(x =distance, y = value/SF, color = variable), data = df2_mean,lwd=1)+
#   scale_color_manual(values = c("#b3cde3","red","orange"))+
#   geom_ribbon(aes(x = distance,ymin = lo/SF, ymax= hi/SF, group = mutations, fill = mutations), data = df2_ci, alpha=0.7)+
#   scale_fill_manual(values = c("red","orange","#b3cde3"))+
#   #geom_errorbar(aes(x = distance,ymin=lo, ymax=hi,group = mutations, fill = mutations),data = df2_ci, width=.1)+
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   scale_y_log10(breaks = trans_breaks("log2", function(x) 2^x),
#                labels = trans_format("log2", math_format(2^.x))) +
#   geom_hline(yintercept = 1,lty=2)+
#   #lims(y=c(20,100))+
#   theme_classic()+#annotation_logticks()+
#   labs(x = "Distance from the CHIP peak (bp)", color = "",y = "Mutation density (/Mb)", title = prefix)
# q
prefixes = c("U87.DNaseH3K27ac.withMYCwa","U87.DNaseH3K27ac.withoutMYCA")

prefixes = c("U87.DNaseH3K4me3.withMYCwa","U87.DNaseH3K4me3.withoutMYCA")

prefixes = c("U87.RPII_and_MYC1k","U87.RPII_no_MYC1k")

prefixes = c("ENCFF870JRW5k.MYCbind","ENCFF870JRW5k.nonMYCbind")

prefixes =c("U87.DNase.MYCbind", "U87.DNase.nonMYCbind")
#--
df = read.delim(paste0("MutationDensity.per150000.",prefixes[1],".txt"))
#df = df[df$bin>=100,]
df$Density = df$Muts*1000000/df$Width
dts = c(1:9,seq(from=1,to=9, by=1)*10,seq(from=1,to=9, by=1)*100, seq(from=1,to=10, by=1)*1000,seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000,5000000,10000000)
df2 = data.frame(distance =dts, 
                 HM_mn=0,HM_lb = 0, HM_ub = 0,
                 NHM_mn=0,NHM_lb = 0, NHM_ub = 0,
                 HMMYC_mn =0, HMMYC_lb =0,HMMYC_ub =0,
                 HMnoMYC_mn =0, HMnoMYC_lb =0, HMnoMYC_ub =0)
for (i in 1:nrow(df2)){
  tmpi = df[df$bin==df2$distance[i],]
  
  d0 = tmpi$Density[tmpi$Type=="HM"]
  df2$HM_mn[i] = mean(d0)
  df2$HM_lb[i] = mean(d0)-sd(d0)
  df2$HM_ub[i] = mean(d0)+sd(d0)
  
  d1 = tmpi$Density[tmpi$Type=="NHM"]
  df2$NHM_mn[i] = mean(d1)
  df2$NHM_lb[i] = mean(d1)-sd(d1)
  df2$NHM_ub[i] = mean(d1)+sd(d1)
  
  d2 = tmpi$Density[tmpi$Type=="HMMYC"]
  df2$HMMYC_mn[i] = mean(d2)
  df2$HMMYC_lb[i] = mean(d2) - sd(d2)
  df2$HMMYC_ub[i] = mean(d2) + sd(d2)
  
  d3 = tmpi$Density[tmpi$Type=="HMnoMYC"]
  df2$HMnoMYC_mn[i] = mean(d3)
  df2$HMnoMYC_lb[i] = mean(d3) - sd(d3)
  df2$HMnoMYC_ub[i] = mean(d3) + sd(d3)
}

scaleFactor = 150000/(2855879637/1000000)


df2_mean_cmycbind = reshape2::melt(df2[,c('distance','HMMYC_mn')],id.vars = 'distance')
df2_mean_cmycbind$bind="cMYc-bind"
df2_ci_cmycbind = data.frame(distance = df2$distance, mutations = 'HMMYC',
                             lo = c(df2$HMMYC_lb), hi = df2$HMMYC_ub)
df2_ci_cmycbind$bind="cMYc-bind"

#
df = read.delim(paste0("MutationDensity.per150000.",prefixes[2],".txt"))
#df = df[df$bin>=100,]
df$Density = df$Muts*1000000/df$Width
dts = c(1:9,seq(from=1,to=9, by=1)*10,seq(from=1,to=9, by=1)*100, seq(from=1,to=10, by=1)*1000,seq(from=20,to=100,by=20)*1000,seq(from=50,to=200,by=50)*10000,5000000,10000000)
df2 = data.frame(distance =dts, 
                 HM_mn=0,HM_lb = 0, HM_ub = 0,
                 NHM_mn=0,NHM_lb = 0, NHM_ub = 0,
                 HMMYC_mn =0, HMMYC_lb =0,HMMYC_ub =0,
                 HMnoMYC_mn =0, HMnoMYC_lb =0, HMnoMYC_ub =0)
for (i in 1:nrow(df2)){
  tmpi = df[df$bin==df2$distance[i],]
  
  d0 = tmpi$Density[tmpi$Type=="HM"]
  df2$HM_mn[i] = mean(d0)
  df2$HM_lb[i] = mean(d0)-sd(d0)
  df2$HM_ub[i] = mean(d0)+sd(d0)
  
  d1 = tmpi$Density[tmpi$Type=="NHM"]
  df2$NHM_mn[i] = mean(d1)
  df2$NHM_lb[i] = mean(d1)-sd(d1)
  df2$NHM_ub[i] = mean(d1)+sd(d1)
  
  d2 = tmpi$Density[tmpi$Type=="HMMYC"]
  df2$HMMYC_mn[i] = mean(d2)
  df2$HMMYC_lb[i] = mean(d2) - sd(d2)
  df2$HMMYC_ub[i] = mean(d2) + sd(d2)
  
  d3 = tmpi$Density[tmpi$Type=="HMnoMYC"]
  df2$HMnoMYC_mn[i] = mean(d3)
  df2$HMnoMYC_lb[i] = mean(d3) - sd(d3)
  df2$HMnoMYC_ub[i] = mean(d3) + sd(d3)
}



df2_mean_nocmycbind = reshape2::melt(df2[,c('distance','HMMYC_mn')],id.vars = 'distance')
df2_mean_nocmycbind$bind="non-cMYc-bind"
df2_ci_nocmycbind = data.frame(distance = df2$distance, mutations = 'HMMYC',
                               lo = c(df2$HMMYC_lb), hi = df2$HMMYC_ub)
df2_ci_nocmycbind$bind="non-cMYc-bind"

#
df2_mean= rbind(df2_mean_cmycbind,df2_mean_nocmycbind)
df2_ci = rbind(df2_ci_cmycbind,df2_ci_nocmycbind)

prefixn = strsplit(prefixes[1],split = "\\.")[[1]][2]
p <- ggplot()+
  geom_line(aes(x =distance, y = value/scaleFactor, color = bind), data = df2_mean,lwd=1)+
  scale_color_manual(values = c("red","orange"))+
  geom_ribbon(aes(x = distance,ymin = lo/scaleFactor, ymax= hi/scaleFactor, group = bind, fill = bind), data = df2_ci, alpha=0.5)+
  scale_fill_manual(values = c("red","orange"))+
  #geom_errorbar(aes(x = distance,ymin=lo, ymax=hi,group = mutations, fill = mutations),data = df2_ci, width=.1)+
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),#limits=c(10^0,2*10^6),
                labels = trans_format("log10", math_format(10^.x))) +
  #scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
  #             labels = trans_format("log10", math_format(10^.x))) +
  geom_hline(yintercept = 1,lty=2)+
  lims(y=c(.5,2.5))+ #
  #expand_limits(x = 1)+
  theme_classic()+#annotation_logticks()+
  theme(axis.text.x = element_text(colour = "black"),axis.text.y = element_text(colour = "black"))+
  labs(x = "Distance from the ChIP peaks (bp)", color = "",y = "Normalized mutation density", title = prefixn)

p
formatC(wilcox.test(value ~ bind, data = df2_mean, paired = T )$p.value, format = "e", digits = 2)
#ggsave(filename = paste0("MutationDensity.per150000.",prefixn,"cMycbindornot.pdf"),plot = p,width = 3.5,height = 2.5,units = 'in',dpi = 300)
