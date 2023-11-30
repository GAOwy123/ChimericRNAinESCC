rm(list=ls())
##画ESCC与TCGA,CCLE交集的融合基因图，首先找到recurrent最高的基因
library(tidyverse)
library(dplyr)
setwd("D:/FusionGene2/1.statistics/3.publicFusion")

Fusion<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957.txt")#ESCC_All_fusion_12133.txt
Fusion<-Fusion %>% filter(Fusion_Type !="cscRNA") %>% dplyr::select(fusionName,index) %>% unique()

dim(Fusion)#23797     3
Fusion_escc_recurrent<-Fusion %>% group_by(fusionName) %>% dplyr::summarise(Recurrent = n()) %>% dplyr::arrange(desc(Recurrent))

Fusion_escc_recurrent<-Fusion_escc_recurrent %>% filter(Recurrent>3)


#load TCGA fusion
TCGA<-read_tsv("D:/FusionGene/16.publicFusions/FusionGDB/TCGA_fusion_information_on_hg19.txt",col_names = F)
TCGA2<-TCGA %>% dplyr::select(X5,X9,X3,X4) %>% unique()
colnames(TCGA2)<-c("G5","G3","Cancer Type","Sample")
TCGA2$fusionName<-paste0(TCGA2$G5,"--",TCGA2$G3)
TCGA<-TCGA %>% dplyr::select(X5,X9,X4) %>% unique()
TCGA_recurrent<-TCGA %>% group_by(X5,X9) %>% dplyr::summarise(Recurrent = n()) %>% dplyr::arrange(desc(Recurrent))

# Fusion_escc_recurrent$FusionGene<-paste0(Fusion_escc_recurrent$G5,"--",Fusion_escc_recurrent$G3)
TCGA_recurrent$fusionName<-paste0(TCGA_recurrent$X5,"--",TCGA_recurrent$X9)
escc.phaseII_TCGA_jiaoji<-left_join(Fusion_escc_recurrent,TCGA_recurrent,by=c("fusionName"="fusionName"))

escc.phaseII_TCGA_jiaoji<-escc.phaseII_TCGA_jiaoji %>% dplyr::select(fusionName,Recurrent.x,Recurrent.y) 
colnames(escc.phaseII_TCGA_jiaoji)<-c("fusionName","ESCC Recurrent","TCGA Recurrent")
write_tsv(escc.phaseII_TCGA_jiaoji,"output/escc.phaseII_recurrent_larger_than_3_TCGA.txt") #109

escc.phaseII_TCGA_jiaoji<-escc.phaseII_TCGA_jiaoji %>% dplyr::arrange(desc(`ESCC Recurrent`),desc(`TCGA Recurrent`))

######
escc.phaseII_TCGA_jiaoji$`TCGA Recurrent`[is.na(escc.phaseII_TCGA_jiaoji$`TCGA Recurrent`)]<-0
escc.phaseII_TCGA_jiaoji$`ESCC Recurrent`<-(-escc.phaseII_TCGA_jiaoji$`ESCC Recurrent`)



##====filter jiaoji fusiongene
# escc<-Fusion %>% filter(fusionName %in% escc.phaseII_TCGA_jiaoji$fusionName) %>% dplyr::select(fusionName,index) %>% unique()
tcga<-TCGA2 %>% filter(fusionName %in% escc.phaseII_TCGA_jiaoji$fusionName) %>% dplyr::select(fusionName,`Cancer Type`,Sample) %>% unique()
tcga<-tcga %>% group_by(fusionName,`Cancer Type`) %>% summarise(recurrent=dplyr::n())


jiaoji<-left_join(escc.phaseII_TCGA_jiaoji,tcga)
jiaoji<-jiaoji[,c("fusionName","ESCC Recurrent","Cancer Type","recurrent")]
colnames(jiaoji)<-c("fusionName","ESCC Samples","Cancer Type","TCGA Samples")

# jiaoji2<-left_join(jiaoji,escc.phaseII_TCGA_jiaoji,by=c("fusionName"="fusionName"))
jiaoji$fusionName<-factor(jiaoji$fusionName,levels=rev(escc.phaseII_TCGA_jiaoji$fusionName))


tmp1<-jiaoji %>% dplyr::select(fusionName,`ESCC Samples`) %>% unique()
tmp1$`Cancer Type`<-"ESCC"
tmp1$Recurrent<-tmp1$`ESCC Samples`
tmp1<-tmp1[,c("fusionName","Cancer Type","Recurrent")]
tmp2<-jiaoji2 %>% dplyr::select(fusionName,`Cancer Type`,`TCGA Samples`) %>% unique() %>% group_by(fusionName,`Cancer Type`) %>% dplyr::summarise(Recurrent=n())
data2<-rbind(tmp1,tmp2)
write_tsv(data2,"output/escc.phaseII_recurrent_larger_than_3_TCGA2.txt") 

data2<-data2[!is.na(data2$`Cancer Type`),]

pdf("output/escc.phaseII_recurrent_larger_than_3_TCGA.pdf",width=8,height=10)
p<-ggplot(data2)+
  geom_bar(aes(x=fusionName,y=Recurrent,fill=`Cancer Type`),
           stat = 'identity',position="stack",
           width = rep(0.8,nrow(data2)))+  #,position = position_dodge(width =0)
  # scale_alpha_discrete(range = c(0.5))+#设置前一个透明度为1，后一个为0.5
  labs(x='Fusion Gene', y = 'Recurrent') +
  coord_flip()+#翻转坐标轴
  # scale_fill_manual(values = c("#CC0000","#006633"))+
  #修改坐标轴刻度和标签
  scale_y_continuous(breaks = seq(-34, 6, 2), 
                     labels = as.character(abs(seq(-34, 6, 2))), 
                     limits = c(-34, 6)) +
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "top")
# panel.background = element_blank())

#现在的版本需要这样设置移除图例啦
#  guides(fill="none",alpha="none")
print(p)
dev.off()

png("output/escc.phaseII_recurrent_larger_than_3_TCGA.png",width=480,height=1000)
print(p)
dev.off()


