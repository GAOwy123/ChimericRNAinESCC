rm(list=ls())
library(tidyverse)
##==============combine fusiongene neoantigen
# dir.create("D:/FusionGene2/6.peptides/neoantigen")
# setwd("D:/FusionGene2/6.peptides/neoantigen")
# setwd("/Users/gaowy/Desktop/work/neoantigen/fusiongene_neoantigen")
setwd("D:/FusionGene2/1.statistics/21.neoantigen/neoantigen/fusiongene_neoantigen")
I<-read_tsv("output/fusiongene_MHCI_neoantigen.txt")
table(I$Tools)
# cscMap FusionCatcher         JAFFA   STAR-Fusion 
# 2997          2208          2146           284
I<- I %>% filter(Tools !="cscMap")
I$group<-"Chimeric RNA"
I<-I %>% dplyr::select(group,`MHCflurry Score`,`MHCnuggetsI Score`,`MHCnuggetsI Score`,
                       `NetMHC Score`,`NetMHCcons Score`,
                       `SMM Score`,`SMMPMBEC Score`)


mutation.neoantigen<-read_tsv("D:/FusionGene2/6.peptides/mutation.neoantigen/input/pvacseq_combined_MHC_Class_I_result.txt")
mutation.neoantigen$group<-"Mutation"

mutation_I<-mutation.neoantigen %>% dplyr::select(group,`MHCflurry MT Score`,`MHCnuggetsI MT Score`,`MHCnuggetsI MT Score`,
                                                  `NetMHC MT Score`,`NetMHCcons MT Score`,
                                                  `SMM MT Score`,`SMMPMBEC MT Score`)
colnames(I)<-colnames(mutation_I)
MHCI<-rbind(I,mutation_I)

MHCI$group<-factor(MHCI$group,levels = c("Chimeric RNA","Mutation"))
library(ggpubr)
dir.create("output/affinity_in_chimericRNA_mutation_removeCscRNA")
MHCI$`SMM MT Score`<-as.numeric(MHCI$`SMM MT Score`)
MHCI$`SMMPMBEC MT Score`<-as.numeric(MHCI$`SMMPMBEC MT Score`)


for(i in colnames(MHCI)[2:ncol(MHCI)]){
  # i<-colnames(MHCI)[2]
  dt<-MHCI[,c("group",i)]
  dt[,i]<-log2(dt[,i])
  p<-ggboxplot(dt,x="group",y=i,color="black",fill="group",
               title = "MHC class I",
               palette = c("#fdd957","#efbb85"))+
    stat_compare_means(label.x = 1.3)+ylim(0,15)+xlab("")+
    theme(legend.position = "none",
          plot.title = element_text(color="black", size=14, hjust = 0.5))
  
  # p
  pdf(paste0("output/affinity_in_chimericRNA_mutation_removeCscRNA/MHCI_",i,".pdf"),width = 4,height = 5)
  print(p)
  dev.off()
}
##+=============================compare neoantigen immunogenicity in chimeric derived and mutation derived=========

II<-read_tsv("output/fusiongene_MHCII_neoantigen.txt")
II<- II %>% filter(Tools !="cscMap")

II$group<-"Chimeric RNA"
II<-II %>% dplyr::select(group,`MHCnuggetsII Score`,`NetMHCIIpan Score`)

mutation.neoantigen<-read_tsv("D:/FusionGene2/6.peptides/mutation.neoantigen/input/pvacseq_combined_MHC_Class_II_result.txt")
mutation.neoantigen$group<-"Mutation"

mutation_II<-mutation.neoantigen %>% dplyr::select(group,`MHCnuggetsII MT Score`,`NetMHCIIpan MT Score`) #,`NNalign MT Score`, `SMMalign MT Score`
colnames(II)<-colnames(mutation_II)
MHCII<-rbind(II,mutation_II)

MHCII$group<-factor(MHCII$group,levels = c("Chimeric RNA","Mutation"))


library(ggpubr)
for(i in colnames(MHCII)[2:ncol(MHCII)]){
  # i<-colnames(MHCI)[2]
  dt<-MHCII[,c("group",i)]
  dt[,i]<-log2(dt[,i])
  p<-ggboxplot(dt,x="group",y=i,color="black",fill="group",
               title = "MHC class II",
               palette = c("#fdd957","#efbb85"))+
    stat_compare_means(label.x = 1.3)+ylim(0,15)+xlab("")+
    theme(legend.position = "none",
          plot.title = element_text(color="black", size=14, hjust = 0.5))
  
  # p
  pdf(paste0("output/affinity_in_chimericRNA_mutation_removeCscRNA/MHCII_",i,".pdf"),width = 4,height = 5)
  print(p)
  dev.off()
}
