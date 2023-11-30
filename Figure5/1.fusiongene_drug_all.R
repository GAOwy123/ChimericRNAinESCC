rm(list=ls())
setwd("D:/FusionGene2/11.drug")
library(tidyverse)


## load drug
gene_drug<-read_tsv("input/oncoKB/oncokb_biomarker_drug_associations.tsv")

fusiongene_drug<-gene_drug  %>% filter(Level %in% c(1,2,3,4))
gene<-unique(fusiongene_drug$Gene)

## load fusion gene
fusionData<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957.changed_name.txt")
fusionData<-fusionData %>% filter(Fusion_Type %in% c("DNA level fusion","RNA level fusion")) 

fusionData_1<-fusionData %>% filter(gene_5 %in% gene)  %>% dplyr::select(gene_5,sample,Fusion_Type)
fusionData_2<-fusionData %>% filter(gene_3 %in% gene)  %>% dplyr::select(gene_3,sample,Fusion_Type)

colnames(fusionData_1)<-c("gene","sample","Fusion_Type")
colnames(fusionData_2)<-c("gene","sample","Fusion_Type")

fusionData<-rbind(fusionData_1,fusionData_2)
fusionData<-unique(fusionData)

#####====================================================================================
##     统计可以治疗的fusion 样本数量     ================================================
#####====================================================================================
sam_count = fusionData %>% dplyr::select(gene,sample) %>% group_by(gene) %>% dplyr::summarise(n=n()) %>% dplyr::arrange(n)
sam_count$gene = factor(sam_count$gene,levels = sam_count$gene)
order_gene<-sam_count$gene
p_sample = ggplot(sam_count, aes(x=gene, y=n, fill="#deedcd"))+
  geom_bar(stat="identity")+
  ylim(0,10)+
  coord_flip()+
  xlab("")+ylab("")+ 
  theme(legend.position = "none",plot.margin = unit(c(0,0,0,0),'lines'), #c(4,1,1,0)
                           axis.line.y = element_blank(),
                           # axis.line.x = element_blank(),
                           # axis.text.x =element_blank(),
                           axis.text.y =element_blank(),
                           axis.ticks.y = element_blank(),
                           axis.title.y=element_blank(),
                           axis.line.x = element_line(),
                           # axis.ticks.x = element_blank(),
                           panel.background = element_blank())


#====================================================================================
fusiongene_drug<-fusiongene_drug %>% filter(Gene %in% fusionData$gene)
fusiongene_drug<-fusiongene_drug %>% dplyr::select(Level, Gene,   Alterations,`Drugs (for therapeutic implications only)`) %>% distinct()

colnames(fusiongene_drug)<-c("Level", "Gene",   "Alterations", "Drugs")
res<-tibble()
for(i in 1:nrow(fusiongene_drug)){
  if(grepl(", ",fusiongene_drug[i,"Drugs"])){
    drug<-unlist(strsplit(unlist(fusiongene_drug[i,"Drugs"]),", "))
    Level<-rep(unlist(fusiongene_drug[i,"Level"]),length(drug))
    Gene<-rep(unlist(fusiongene_drug[i,"Gene"]),length(drug))
    Alterations<-rep(unlist(fusiongene_drug[i,"Alterations"]),length(drug))
    Drugs<-drug
    res<-rbind(res,tibble(Level,Gene,Alterations,Drugs))
  }else{
    Level<-unlist(fusiongene_drug[i,"Level"])
    Gene<-unlist(fusiongene_drug[i,"Gene"])
    Alterations<-unlist(fusiongene_drug[i,"Alterations"])
    Drugs<-unlist(fusiongene_drug[i,"Drugs"])
    res<-rbind(res,tibble(Level,Gene,Alterations,Drugs))
    
  }

}

dir.create("output")
write_tsv(res,"output/gene_drug.tsv")

##########==================================================================================
##########===================draw the potential fusion gene drug ======================
##########==================================================================================
res$Gene<-factor(res$Gene,levels=order_gene)

pdf("output/oncoKB_fusiongene_drug1.pdf",width=7,height = 5)
onco <-ggplot(res, aes(Drugs, Gene)) + 
  geom_tile(aes(fill=Level))+xlab("")+ylab("")+
  scale_fill_manual(values=c("1" = "#f3aab4","2" = "#f9cf99","3"="#b596c5","4"="#74c9d4"),na.value = "white")

onco<-onco + 
  theme_bw()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines'),
        axis.text.y= element_text(size = 18),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size =18,vjust = 1, hjust = 0.9, angle = 30),
        panel.grid.major=element_line(colour=NA)
  )+
  theme(legend.position = "top")
onco
dev.off()

onco <-ggplot(res, aes(Drugs, Gene)) + 
  geom_tile(aes(fill=Level))+xlab("")+ylab("")+
  scale_fill_manual(values=c("1" = "#f3aab4","2" = "#f9cf99","3"="#b596c5","4"="#74c9d4"),na.value = "white")

onco<-onco + 
  theme_bw()+
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),'lines'),
        axis.text.y= element_text(size = 18),
        axis.ticks = element_blank(),
        axis.text.x = element_text(size =18,vjust = 1, hjust = 0.9, angle = 30),
        panel.grid.major=element_line(colour=NA)
  )+
  theme(legend.position = "none")




####====================================================================================
###===================draw the potential fusion gene drug ==============================
###=====================================================================================

onco = ggplot_gtable(ggplot_build(onco))
p_sample = ggplot_gtable(ggplot_build(p_sample))


library(grid)
library(gridExtra)


m<-cowplot::plot_grid(onco,p_sample,ncol = 2, align = "vh",axis = "l",rel_widths = c(2, 1)) ##垂直水平对齐坐标轴
ggsave("output/oncoKB_fusiongene_drug.pdf",m,device = "pdf",width=15,height=5,limitsize = FALSE)

