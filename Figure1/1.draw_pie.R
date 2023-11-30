rm(list=ls())
setwd("D:/FusionGene2/1.statistics/8.pie_fusion")
library(tidyverse)
library(plyr)
library(ggsci)
library(ggplot2)
library(gridExtra)
library(patchwork)

Fusion<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957_with_pheno.txt")


fusionType<-Fusion %>% dplyr::select(fusionName,Fusion_Type) %>% distinct()
dt<-as.data.frame(table(fusionType$Fusion_Type))
colnames(dt)<-c("Type","Freq")
dt$ratio<-dt$Freq/sum(dt$Freq)


fusionNumberPersample<-Fusion %>% dplyr::select(fusionName,sample,Fusion_Type) %>% distinct()
fusionNumberPersample$id<-1
fusionNumberPersample2<-fusionNumberPersample %>% dplyr::select(sample,Fusion_Type) %>% 
  group_by(Fusion_Type,sample)  %>% dplyr::summarise(dplyr::n())

cscRNA<-fusionNumberPersample %>% filter(Fusion_Type=="cscRNA")
nrow(cscRNA)/length(unique(cscRNA$sample))

DNA<-fusionNumberPersample %>% filter(Fusion_Type=="DNA level fusion")
nrow(DNA)/length(unique(DNA$sample))


RNA<-fusionNumberPersample %>% filter(Fusion_Type=="RNA level fusion")
nrow(RNA)/length(unique(RNA$sample))

dt$perSampleNumber<-c(nrow(cscRNA)/length(unique(cscRNA$sample)),nrow(DNA)/length(unique(DNA$sample)),nrow(RNA)/length(unique(RNA$sample)))


dt$label<-paste0(paste0(dt$Type,"\n",
                        round(dt$ratio,3)*100,"%\n",
                        dt$Freq,"events","\n",
                        "(",round(dt$perSampleNumber,1),"/Sample",")"))


library(ggplot2)

#画圆环图
library(ggsci)#+scale_color_npg()+scale_fill_npg()

mycolor1<-pal_npg("nrc", alpha = 0.5)(8)
mycolor2<-pal_aaas(alpha = 0.4)(8)

cols<-list("DNA level fusion"="#bc9584",
           "RNA level fusion"="#FA7F6F",
           "cscRNA"="#82B0D2")
library(ggpubr)
p<-ggpie(dt, "ratio",
         label = "label",                               
         fill = "Type",                           
         color = "white",                              
         palette = cols)
p
pdf("pie_fusion_new.pdf",width=5,height = 5)
print(p)
dev.off()


