#change working directory
library(rstudioapi)
present_folder = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(present_folder)
print(getwd())


mycsnp = read.delim('rs55705857_frequency.tsv')
mycsnp =mycsnp[  mycsnp$Population!='Global'& mycsnp$Population!='Other'& mycsnp$Population!='Total',]
mycsnp = mycsnp[grepl('East Asia', mycsnp$Population) | grepl('Europe', mycsnp$Population),]
library(ggplot2)


ggplot(mycsnp, aes(x = paste(Population,Study), y = 100*AF2))+
  geom_bar(stat = 'identity', aes(fill = Population !='East Asian'))+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_brewer(palette = 7, type = 'qual')+
  theme_classic()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  labs(y= 'Risk allele frequency (%)', x = '', fill =  'Population')
#ggsave(file = 'mycsnpAF.pdf', width = 3, height = 3.9)
