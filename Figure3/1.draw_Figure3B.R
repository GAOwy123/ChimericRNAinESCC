rm(list=ls())

setwd("D:/FusionGene2/1.statistics/21.neoantigen/neoantigen/fusiongene_neoantigen/Figure3B")
library(tidyverse)
dt<-read_tsv("neoantigen.txt")
dt<-dt %>% dplyr::select(`MT Epitope Seq`,Type) %>% unique()

dt2<-dt %>% dplyr::group_by(Type) %>% summarise(count=n())
dt2$ratio<-dt2$count/sum(dt2$count)

library(ggpubr)
p<-ggpie(dt2,x="ratio",label = paste0(dt2$Type," (",round(dt2$ratio,3)*100,"%",")"),fill="Type",
      palette = c("#a7c957","#e29578","#2ec4b6","#ff9f1c","#457b9d","#f07167","#edae49"),
      lab.pos = c ( "out" )
      )+
  theme(legend.position = "bottom",legend.text = element_text(size=8))
p

pdf("Figure3B_pie.pdf",width=5,height = 5)
print(p)
dev.off()

