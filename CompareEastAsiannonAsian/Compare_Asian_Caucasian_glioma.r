#change working directory
library(rstudioapi)
present_folder = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(present_folder)
print(getwd())

options(stringsAsFactors = F)
library(ggplot2)
library(cowplot)
ini = read.delim('Pairedgliomas.genomics_initial.txt', na.strings = c('NA','#N/A'))
rec = read.delim('Pairedgliomas.genomics_recurrence.txt', na.strings = c('NA','#N/A'))
c710 = read.delim('panglioma.chr710.tsv')
ini$chr7gain10loss = c710$new_chr710_ini[match(ini$Patient_ID, c710$Patient_ID)]
rec$chr7gain10loss = c710$new_chr710_rec[match(rec$Patient_ID, c710$Patient_ID)]

ini = ini[order(ini$Patient_ID),]; rec = rec[order(rec$Patient_ID),]
ini$TP53del[ini$Cohorts=='Yale'] = ini$MDM2amp[ini$Cohorts=='Yale'] = ini$MDM4amp[ini$Cohorts=='Yale'] =ini$CDKN2Adel[ini$Cohorts=='Yale'] =ini$CDK4amp[ini$Cohorts=='Yale'] =ini$RB1del[ini$Cohorts=='Yale'] =ini$EGFRamp[ini$Cohorts=='Yale'] =ini$PDGFRAamp[ini$Cohorts=='Yale'] =ini$METamp[ini$Cohorts=='Yale'] =ini$NF1del[ini$Cohorts=='Yale'] =ini$PTENdel[ini$Cohorts=='Yale'] =NA
rec$TP53del[rec$Cohorts=='Yale'] = rec$MDM2amp[rec$Cohorts=='Yale'] = rec$MDM4amp[rec$Cohorts=='Yale'] =rec$CDKN2Adel[rec$Cohorts=='Yale'] =rec$CDK4amp[rec$Cohorts=='Yale'] =rec$RB1del[rec$Cohorts=='Yale'] =rec$EGFRamp[rec$Cohorts=='Yale'] =rec$PDGFRAamp[rec$Cohorts=='Yale'] =rec$METamp[rec$Cohorts=='Yale'] =rec$NF1del[rec$Cohorts=='Yale'] =rec$PTENdel[rec$Cohorts=='Yale'] =NA
stopifnot(identical(ini$Patient_ID, rec$Patient_ID))
names(ini)[names(ini)=="chr17p_NLOH"] = "chr17p_CNLOH"
names(rec)[names(rec)=="chr17p_NLOH"] = "chr17p_CNLOH"
gns = c("IDH1.2","TP53","chr17p_CNLOH","ATRX","codel","CIC","FUBP1", 
        "MDM2amp","MDM4amp","CDKN2Adel","CDK4amp","RB1","RB1del",
        "EGFR","EGFRamp","PDGFRA","PDGFRAamp","METamp","F3T3","ZM","METex14","EGFRvIII","NF1","NF1del",
        "PTEN","PIK3CA","PIK3CG","PIK3R1",
        "chr7gain10loss", "MYC_gain")#
#br=read.delim('~/Documents/pangliomaevolution/Rev1/Novelty/TreatmentEffect/glioma.branch.evolution.txt',sep=' ')
inim = ini[,gns]
recm = rec[,gns]
p1m = inim
for (i in 1:nrow(p1m)){
  for (j in 1:ncol(p1m)){
    a = inim[i,j]; b = recm[i,j]
    if (is.na(a) | is.na(b)){
      p1m[i,j]=NA
    }else if (a==1 & b ==1){
      p1m[i,j] = 3
    } else if (a==1 & b ==0){
      p1m[i,j] = 1
    }else if (a==0 & b ==1){
      p1m[i,j] = 2
    }else if (a==0 & b ==0) {
      p1m[i,j] = 0
    }else{
      print(paste('WARNING: weired value detected at coordinate',i,j))
    }
  }
}
p1a = data.frame(Subtype = paste0(ini$IDH,ini$X1p19qcodel),
                 Grade_I = ini$WHO2016_grade_initial,
                 Grade_R = ini$WHO2016_grade_recurrence,
                 TMZ_R = ifelse(is.na(rec$TMZ_R),NA,ifelse(rec$TMZ_R==1,4,0)),
                 MGMTfusion = ifelse(is.na(rec$MGMTfusion),NA,ifelse(rec$MGMTfusion==1,4,0)),
                 Hypermutation = ifelse(is.na(rec$HM_R),NA,ifelse(rec$HM_R=="YES"|rec$HM_R=="Yes",4,0)),
                 MMR = ifelse(is.na(rec$MMR),NA,ifelse(rec$MMR==1,4,0)),
                 #MYC_gain = ifelse(is.na(ini$MYC_gain),NA,ifelse(ini$MYC_gain==1,4,0)),
                 stringsAsFactors = F)

rownames(p1a) = ini$Patient_ID
p1a$Subtype[p1a$Subtype=='IDHwtnoncodel']='IDHwt'
p1a$Subtype[p1a$Subtype=='IDHmutnoncodel']='IDHnon'
p1a$Subtype[p1a$Subtype=='IDHmutcodel']='IDHcodel'


clin = ini
clin$Race[is.na(clin$Race)] = 'Others'
clin$Race[clin$Race=="non-Asian"] = "Others"
clin$Race[which(clin$Race=='Asian')] = 'EastAsian'

names(clin)[which(names(clin)=='X')] = 'Subtype'
p1a$Cohort = clin$Cohorts[match(rownames(p1a), clin$Patient_ID)]
p1a$New =ifelse(p1a$Cohort %in% c('CGGA','Tiantan','CUHK','SMC'),"New",'Published')
p1a$Race = clin$Race[match(rownames(p1a),clin$Patient_ID)]

#        compare driver mutations in Asian and non-Asian by subtype         #
library(ggplot2)
library(ggpubr)
library(ggbeeswarm)


lowPurity = c("PS003","PS006","PS011","PS017","PS019","PS028","PS065","PS069","PS095","PS097","PS103",
              "PS104","PS119","PS123","PS269","PS274","PS278","PS284","PS287")
stopifnot(identical(ini$Patient_ID, rownames(p1a)))
idx = which((!(is.na(ini$IDH1.2))) & (!(is.na(p1a$Race))));
p1a = p1a[idx,]; p1m = p1m[idx,] #all Asian

#idx = which(ini$Cohorts=='GLSS'); p1a = p1a[idx,]; p1m = p1m[idx,] #GLASS
#idx = which(startsWith(rownames(p1a),"PSX")); p1a = p1a[idx,]; p1m = p1m[idx,] #highqual203



p1m = as.data.frame(t(p1m)) 
idx = order(p1a$Subtype,p1a$Grade_R, p1a$Grade_I, decreasing = T)
p1a = p1a[idx,]
p1m = p1m[,idx]
#p1mmyc = as.data.frame(t(p1m));# p1m = p1m[rownames(p1m)!="MYC_gain",] #MYC_gain is moved to annotation

#
# scdf = read.delim('Asian_nonAsian_score_20220324.txt')
# scdf = scdf[order(scdf$Asian_Cauca_Dif_Founding.Events),]
# scdf = scdf[!scdf$gene %in% c("IDH1.2","codel"),]

scdf=data.frame(gene = c("chr7gain10loss","MYC_gain","PI3K_pathway","PIK3R1","PIK3CA","PIK3CG","PTEN","RTK_pathway","NF1del","NF1","EGFRvIII","METex14","ZM","F3T3","METamp","PDGFRAamp","PDGFRA","EGFRamp","EGFR","cellcycle_pathway","RB1del","RB1","CDKN2Adel","CDK4amp","MDM2amp","MDM4amp","Hypermutation","MMR","FUBP1","CIC","p53_pathway","chr17p_CNLOH","ATRX","TP53"))
##
# IDHwt
##
table( p1a$Subtype,p1a$Race)
p1mt0 = as.data.frame(t(p1m))
gns = c("CIC","FUBP1","ATRX","TP53","chr17p_CNLOH","MDM2amp","MDM4amp",#"IDH1.2","codel",
        "CDKN2Adel","CDK4amp","RB1del","RB1",
        "EGFRamp","EGFR","EGFRvIII","PDGFRAamp","PDGFRA","METamp","ZM","METex14","F3T3",
        "NF1","NF1del","PTEN","PIK3CA","PIK3CG","PIK3R1","chr7gain10loss",
        "MYC_gain","Hypermutation","MMR")
gns2 = c("CIC","FUBP1","ATRX","TP53","chr17p_CNLOH","MDM2amp","MDM4amp","p53_pathway", #"IDH1.2","codel",
         "CDKN2Adel","CDK4amp","RB1del","RB1","cellcycle_pathway",
         "EGFRamp","EGFR","EGFRvIII","PDGFRAamp","PDGFRA","METamp","ZM","METex14","F3T3","NF1","NF1del","RTK_pathway",
         "PTEN","PIK3CA","PIK3CG","PIK3R1", "PI3K_pathway",
         "chr7gain10loss",
         "MYC_gain","Hypermutation","MMR")
rc = "EastAsian";sbt = "IDHwt"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "EastAsian"
#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "Others";sbt = "IDHwt"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}


p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "Others"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)

p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st_idhwt = p1st
#
library(reshape2)
tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='Others']=-tmp$value[tmp$race=='Others']
tmp$race = factor(tmp$race, levels = c("Others",'EastAsian'))

tmp_idhwt = tmp
gnorder = scdf$gene#p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
gnorder0 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}

df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhwt = df2
df2 = reshape2::melt(df2, id.vars = c('gene'))
library(ggplot2)
p1<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene)& tmp$status!='wtp',])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c("#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p1_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.05), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p1,p1_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
#write.table(p1st_idhwt, file = 'p1st_idhwt.txt',quote = F,sep = "\t")
#
# now IDHmut-noncodel
#
rc = "EastAsian";sbt = "IDHnon"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "EastAsian"
#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "Others";sbt = "IDHnon"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}

p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "Others"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)


p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st_idhnon = p1st
#

tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='Others']=-tmp$value[tmp$race=='Others']
tmp$race = factor(tmp$race, levels = c("Others",'EastAsian'))
tmp_idhnon = tmp
gnorder = scdf$gene #p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
#gnorder1 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}
df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhnon = df2
df2 = reshape2::melt(df2, id.vars = c('gene'))

p2<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene)& tmp$status!='wtp',])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c("#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0, 0, 0), "cm"))

p2_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.05), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p2,p2_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
#write.table(p1st_idhnon, file = 'p1st_idhnon.txt',quote = F, sep = "\t")
#
# now IDHmut-codel
#
rc = "EastAsian";sbt = "IDHcodel"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] 
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_Asian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_Asian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_Asian) = names(p1mt)
for (i in 1:nrow(p1st_Asian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_Asian[i,1] = sum(x==0)
  p1st_Asian[i,2] = sum(x==1)
  p1st_Asian[i,3] = sum(x==2)
  p1st_Asian[i,4] = sum(x==3)
}

p1st_Asian$gene = rownames(p1st_Asian)
p1st_Asian$race = "EastAsian"

#hm: only consider TMZ treated cases
p1st_Asian$wt[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_Asian$wt[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_Asian$rec[p1st_Asian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#

p1st_Asian$alt_ini = (p1st_Asian$ini+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)
p1st_Asian$alt_rec = (p1st_Asian$rec+p1st_Asian$shr)/(p1st_Asian$ini+p1st_Asian$shr+p1st_Asian$wt+p1st_Asian$rec)

#
rc = "Others";sbt = "IDHcodel"
p1mt = p1mt0[p1a$Race==rc & p1a$Subtype==sbt,] #ignore MYC gain
p1asub = p1a[p1a$Race==rc & p1a$Subtype==sbt,]
p1mt = cbind(p1mt,p1asub[,c("Hypermutation","MMR")])
p1mt$MMR[which(p1mt$MMR==4)]=2; p1mt$Hypermutation[which(p1mt$Hypermutation==4)]=2
p53ini = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('ATRX','TP53','MDM2amp','MDM4amp')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$p53_pathway = p53ini+2*p53rec


p53ini = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('CDKN2Adel','CDK4amp','RB1del','RB1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$cellcycle_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('EGFRamp','EGFR','EGFRvIII','PDGFRAamp','PDGFRA','METamp','ZM','METex14','F3T3','NF1','NF1del')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$RTK_pathway = p53ini+2*p53rec

p53ini = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] %%2,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53rec = apply(p1mt[,c('PTEN','PIK3CA','PIK3CG','PIK3R1')] >1,1,FUN = function(x) ifelse(sum(is.na(x))>0 & sum(x,na.rm = T)==0, NA, sum(x,na.rm = T)))
p53ini[p53ini>0]=1;  p53rec[p53rec>0]=1
p1mt$PI3K_pathway = p53ini+2*p53rec
p1st_nonAsian = data.frame(matrix(nrow = ncol(p1mt), ncol = 4))
names(p1st_nonAsian) = c("wt", "ini", "rec", 'shr')
rownames(p1st_nonAsian) = names(p1mt)
for (i in 1:nrow(p1st_nonAsian)){
  x = p1mt[,i]
  x = x[!is.na(x)]
  p1st_nonAsian[i,1] = sum(x==0)
  p1st_nonAsian[i,2] = sum(x==1)
  p1st_nonAsian[i,3] = sum(x==2)
  p1st_nonAsian[i,4] = sum(x==3)
}

p1st_nonAsian$gene = rownames(p1st_nonAsian)
p1st_nonAsian$race = "Others"
#hm: only consider TMZ treated cases
p1st_nonAsian$wt[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='Hypermutation'] = table(p1asub$Hypermutation[p1asub$TMZ_R==4])[2]
p1st_nonAsian$wt[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[1]
p1st_nonAsian$rec[p1st_nonAsian$gene=='MMR'] = table(p1asub$MMR[p1asub$TMZ_R==4])[2]
#
p1st_nonAsian$alt_ini = (p1st_nonAsian$ini+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)
p1st_nonAsian$alt_rec = (p1st_nonAsian$rec+p1st_nonAsian$shr)/(p1st_nonAsian$ini+p1st_nonAsian$shr+p1st_nonAsian$wt+p1st_nonAsian$rec)

p1st = rbind(p1st_Asian, p1st_nonAsian)
p1st$inip = p1st$ini/(p1st$wt+p1st$shr+p1st$ini)
p1st$recp = p1st$rec/(p1st$wt+p1st$shr+p1st$rec)
p1st$shrp = p1st$shr/(p1st$wt+p1st$shr+p1st$ini)
p1st$wtp = p1st$wt/(p1st$wt+p1st$shr+p1st$ini)
p1st[is.na(p1st)]=0
p1st_idhcod = p1st
#

tmp = melt(p1st[,c('wtp','inip','recp','shrp','gene','race')], id.vars = c('gene','race'),variable.name = 'status')
tmp$xlab = paste(tmp$gene,tmp$race, sep = '.') #lost during evolution

tmp$value[tmp$race=='Others']=-tmp$value[tmp$race=='Others']
tmp$race = factor(tmp$race, levels = c("Others",'EastAsian'))
tmp_idhcod = tmp
gnorder = scdf$gene #p1st$gene[order(p1st$alt_rec[p1st$race=='non-Asian'] - p1st$alt_rec[p1st$race=='Asian'],decreasing = T)]
tmp$gene = factor(tmp$gene, levels = gnorder)
#gnorder2 = gnorder
#tmp = tmp[!tmp$gene %in% c('IDH1.2','codel','CIC','FUBP1','METex14','ZM','F3T3','TERTp'),]
#tmp$gene = factor(tmp$gene)

#
df2 = data.frame(gene = rev(levels(tmp$gene)), p_alt = 1, p_evol = 1,p_overall = 1, stringsAsFactors = F)
for (ix in 1:nrow(df2)){ #comapre Asian and non-Asian
  tb = p1st[p1st$gene==df2$gene[ix],1:4]
  tb$alt = tb[,2]+tb[,3]+tb[,4]
  p = fisher.test(tb[,c(1,5)])$p.value 
  df2$p_alt[ix] = signif(p,2)
  tb$statuschange = tb[,2]+tb[,3]
  tb$nochange = tb[,1]+tb[,4]
  q = fisher.test(tb[,c(7,6)])$p.value
  df2$p_evol[ix] = signif(q,2)
  r = fisher.test(tb[,c(1,4,6)])$p.value
  r = fisher.test(tb[,1:4])$p.value
  df2$p_overall[ix] = signif(r,2)
}
df2$gene = factor(df2$gene, levels = levels(tmp$gene))
df2_idhcod=df2
df2 = reshape2::melt(df2, id.vars = c('gene'))

p3<-ggplot(mapping = aes(x = gene, y = 100*value),data = tmp[!is.na(tmp$gene)& tmp$status!='wtp',])+
  geom_bar(mapping = aes(fill = status),stat = 'identity', show.legend = F)+
  geom_hline(yintercept = 0, size= 0.25)+
  scale_fill_manual(values = c("#d02a7c","#010078","#e2b449"))+
  scale_y_continuous(breaks=c(-100,-50,0,50,100),labels = c(100,50,0,50,100),limits = c(-100,100))+
  coord_flip()+#facet_grid( ~ race,scales = 'free_x')+
  theme_bw()+labs(x = '', y = '',title = paste0(sbt,'') )+
  theme( axis.text = element_text(colour = 'black',size=6),axis.ticks.y = element_blank(),#panel.grid  = element_blank(),
         plot.title = element_text(hjust = 0.5),plot.margin = unit(c(0, 0.2, 0, 0), "cm"))

p3_2<- ggplot(mapping = aes(x = gene, y = 90,label = value),data = df2[df2$variable=='p_overall',])+
  geom_text(size=2,position = 'identity',aes(color = value<0.05), show.legend = F)+
  scale_color_manual(values = c('black','red'))+
  coord_flip()+#facet_grid( ~ variable,scales = 'free_x')+
  theme_bw()+labs(y = '', x = '' )+
  theme(axis.ticks = element_blank(),axis.text = element_blank(),
        panel.grid.minor.x  = element_blank(),panel.grid.major.x  = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "cm"))

cowplot::plot_grid(p3,p3_2,align = 'h',rel_widths = c(0.8,0.2),nrow = 1)
#write.table(p1st_idhcod, file = 'p1st_idhcod.txt',quote =F, sep = "\t")
#

pdf('Fig1c_EastAsian_others_comparison.pdf', width = 12, height = 5)
cowplot::plot_grid(p1, p2+theme(axis.text.y = element_blank()),
                   p3+theme(axis.text.y = element_blank()),
                   align = 'h',nrow = 1, rel_widths = c(0.35,0.25,0.25))

dev.off()
