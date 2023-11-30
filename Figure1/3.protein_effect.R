# 展示，所有fusion events能形成蛋白值的比例有多少，
# 以及形成蛋白的fusion 产生的effect的比例是多少
rm(list=ls())
library(ggpubr)
library(tidyverse)
setwd("D:/FusionGene2/5.kinase/draw_kinase_domian")
dir.create("output/protein")

coding_protein<-read_tsv("input/fusion_protein.combine2.txt")
coding_protein<-coding_protein %>% separate(event,c("G5","B5","G3","B3"),sep="-|_",remove = F)
coding_protein<-coding_protein %>% filter(G5!=G3)

coding_protein$effect<-gsub("in-frame (with mutation)","in-frame",coding_protein$effect,fixed = T)
effect_percentage<-coding_protein %>% dplyr::select(event,effect) %>% group_by(effect) %>% summarise(Number=n())
effect_percentage$Percentage<-round(effect_percentage$Number/sum(effect_percentage$Number)*100,digits = 2)
labs<-paste0(effect_percentage$effect," (",effect_percentage$Percentage,"%)")




pdf("D:/FusionGene2/1.statistics/20.FusionProtein/output/protein_effect.pdf",width=4,height = 4)
pie(effect_percentage$Number, border="white", col=c("#9fcc58","#da87b6"), label=labs,
    main="")
dev.off()

write_tsv(effect_percentage,"D:/FusionGene2/1.statistics/20.FusionProtein/output/protein_effect.txt")
