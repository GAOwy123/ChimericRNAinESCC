###===============================================================
###=====================draw immune inflitration==================
###===============================================================
rm(list=ls())
library(tidyverse)
library(ggplot2)

rm(list=ls())
setwd("D:/FusionGene2/1.statistics/21.neoantigen/neoantigen/fusiongene_neoantigen")
# setwd("/Users/gaowy/Desktop/work/neoantigen/fusiongene_neoantigen")
dir.create("immune_inflitration_removeCscRNA")


###chimeric RNA derived neoantigen number
neoantigen<-read_tsv("output/fusiongene_MHCI_MHCII_neoantigen.txt")
colnames(neoantigen)<-c("Patient",colnames(neoantigen)[2:ncol(neoantigen)])

neoantigen<-neoantigen %>% filter(Fusion_Type !="cscRNA")


neoantigen<-neoantigen %>% dplyr::select(Patient,`MT Epitope Seq`) %>% distinct()



sample_sum<-neoantigen %>% group_by(Patient) %>% summarise(sum=dplyr::n())

phenoData<-read_tsv("D:/FusionGene2/8.survival/input/ESCC_155_phenotype_20210512.txt")
phenoData<-phenoData %>% dplyr::select(sample)
group<-right_join(sample_sum,phenoData,by=c("Patient"="sample"))
# group$sum[is.na(group$sum)]<-0
group<-group[!is.na(group$sum),]
#####select cutpoint

# group$group<-ifelse(group$sum>quantile(group$sum)[4],"High",
#                     ifelse(group$sum<quantile(group$sum)[2],"Low","Medium"))
# group<-group %>% dplyr::filter(group !="Medium")

group<-group %>% dplyr::arrange(desc(sum))

group1<-head(group,20)
group1$group<-"High"
group2<-tail(group,20)
group2$group<-"Low"
group<-rbind(group1,group2)

group<-group %>% dplyr::select(Patient,group)

colnames(group)<-c("sample","group")

estimate.result<-read.table("D:/work/eRNA/0.rawData/immune_inflitration/estimate/ESCC_estimate_score.txt",row.names = 1)
estimate.result$sample<-gsub("^X","",rownames(estimate.result))

data<-inner_join(group,estimate.result)


library(ggplot2)
library(ggpubr)
data<-data[,-c(1)]

data1<-data[,c("group","StromalScore","ImmuneScore")]
colnames(data1)<-c("group","Stromal","Immune")

data1<-gather(data1,key="type",value = "Score",-group)
#data1$group<-ifelse(data1$group=="1","type 1","type 2")
data1$group<-as.factor(data1$group)
###===============================================================
####================  draw Figure3A   ============================
###===============================================================
p <- ggboxplot(data1, x = "type", y = "Score",
               color = "group", palette=c("#FFBE7A","#8ECFC9","#f08d85","#46b566","#f1a327","#e082b2"),#palette = "jco",
               # add = "jitter",
               xlab = "")
p <-p + stat_compare_means(aes(group = group))

dir.create("immune_inflitration_removeCscRNA")
pdf("immune_inflitration_removeCscRNA/stromal_immune_score.pdf",width = 4,height=4)
print(p)
dev.off()

tiff("immune_inflitration_removeCscRNA/stromal_immune_score.tiff",width = 480,height=480)
print(p)
dev.off()
###===============================================================
####================  draw Figure3B   ============================
###===============================================================
data2<-data[,c("group","TumorPurity")]
#data2$group<-ifelse(data2$group=="1","type 1","type 2")
data2$group<-as.factor(data2$group)
p2 <- ggboxplot(data2, x = "group", y = "TumorPurity",
                color = "group", palette=c("#FFBE7A","#8ECFC9","#f08d85","#46b566","#f1a327","#e082b2"),
                # add = "jitter",
                xlab = "",ylab="Tumor purity")

p2 <- p2+stat_compare_means(aes(group = group))
pdf("immune_inflitration_removeCscRNA/TumorPurity.pdf",width = 6,height=6)
print(p2)
dev.off()

tiff("immune_inflitration_removeCscRNA/TumorPurity.tiff",width = 480,height=480)
print(p2)
dev.off()

###===============================================================
####=====  draw Figure3C   ============================
###===============================================================
# group<-read_tsv("consensusclusterplus/2_hc_spearman_consensus_cluster_group.txt",col_names = T)

estimate.result<-read.table("D:/work/eRNA/0.rawdata/immune_inflitration/estimate/ESCC_estimate_score.txt",row.names = 1)
estimate.result$sample<-gsub("^X","",rownames(estimate.result))
estimate.result<-estimate.result[,c(1,2,4,5)]

cibersort.result<-read.table("D:/work/eRNA/0.rawdata/immune_inflitration/cibersort/tumor_cibersort.result.txt",row.names = 1,sep="\t")
cibersort.result<-cibersort.result[,c(1:22)]
cibersort.result$sample<-gsub("^X","",rownames(cibersort.result))


xcell.result<-read.table("D:/work/eRNA/0.rawdata/immune_inflitration/xcell/xCell_result.tsv",row.names = 1,sep="\t")
xcell.result<-t(xcell.result)
xcell.result<-as.data.frame(xcell.result)
xcell.result<-xcell.result[,1:64]
rownames(xcell.result)<-gsub(".","-",rownames(xcell.result),fixed = T)
xcell.result$sample<-gsub("^X","",rownames(xcell.result))


ssGSEA.result<-read.table("D:/work/eRNA/0.rawdata/immune_inflitration/ssGSEA/ssGSEA_result.tsv",row.names = 1,sep="\t")
ssGSEA.result<-t(ssGSEA.result)
ssGSEA.result<-as.data.frame(ssGSEA.result)
#rownames(ssGSEA.result)<-gsub(".","-",rownames(ssGSEA.result),fixed = T)
ssGSEA.result$sample<-gsub("^X","",rownames(ssGSEA.result))

###===============================================================
#####======Figure3C  draw cibersort result ===============================
###===============================================================

data<-inner_join(estimate.result,cibersort.result)
data<-inner_join(group,data)
data<-data %>% dplyr::arrange(group)

data<-as.data.frame(data)
# data[,4:ncol(data)]<-apply(data[,4:ncol(data)],2,as.numeric)
data<-data[,c(2,6:ncol(data))]
rownames(data)<-data$sample
data1<-reshape2::melt(data,id=1)
colnames(data1)<-c("group","celltype","immune_inflitration")
p <- ggboxplot(data1, x = "celltype", y = "immune_inflitration",
               color = "group", palette=c("#FFBE7A","#8ECFC9","#f08d85","#46b566","#f1a327","#e082b2"),
               # add = "jitter",
               xlab = "",ylab="Immune inflitration",legend.title="Neoantigen Number",title = "cybersort")
p <-p +stat_compare_means(aes(group = group),label="p.signif")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust=1,colour = "black"),
        panel.grid = element_blank(),
        axis.text.y =element_text(colour = "black"),
        legend.position = "top",
        plot.title = element_text(hjust=0.5)
  )


pdf("immune_inflitration_removeCscRNA/cybersort_boxplot.pdf",width = 10,height=6)
p
dev.off()

tiff("immune_inflitration_removeCscRNA/cybersort_boxplot.tiff",width = 800,height=480)
p
dev.off()



###===============================================================
#####======Figure3C  draw xcell result ===============================
###===============================================================

data<-inner_join(estimate.result,xcell.result)
data<-inner_join(group,data)
data<-data %>% dplyr::arrange(group)

data<-as.data.frame(data)
# data[,4:ncol(data)]<-apply(data[,4:ncol(data)],2,as.numeric)
data<-data[,c(2,6:ncol(data))]
rownames(data)<-data$sample
data1<-reshape2::melt(data,id=1)
colnames(data1)<-c("group","celltype","immune_inflitration")
p <- ggboxplot(data1, x = "celltype", y = "immune_inflitration",
               color = "group", palette=c("#FFBE7A","#8ECFC9","#f08d85","#46b566","#f1a327","#e082b2"),
               # add = "jitter",
               xlab = "",ylab="Immune inflitration",legend.title="Neoantigen Number",title = "xcell")
p <-p +stat_compare_means(aes(group = group),label="p.signif")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust=1,colour = "black"),
        panel.grid = element_blank(),
        axis.text.y =element_text(colour = "black"),
        legend.position = "top",
        plot.title = element_text(hjust=0.5)
  )

pdf("immune_inflitration_removeCscRNA/xcell_boxplot.pdf",width = 12,height=5)
p
dev.off()

tiff("immune_inflitration_removeCscRNA/xcell_boxplot.tiff",width = 1200,height=500)
p
dev.off()



###===============================================================
#####======Figure3C  draw ssGSEA result ===============================
###===============================================================

data<-inner_join(estimate.result,ssGSEA.result)
data<-inner_join(group,data)
data<-data %>% dplyr::arrange(group)

data<-as.data.frame(data)
# data[,4:ncol(data)]<-apply(data[,4:ncol(data)],2,as.numeric)
data<-data[,c(2,6:ncol(data))]
rownames(data)<-data$sample
data1<-reshape2::melt(data,id=1)
colnames(data1)<-c("group","celltype","immune_inflitration")
p <- ggboxplot(data1, x = "celltype", y = "immune_inflitration",
               color = "group", palette=c("#FFBE7A","#8ECFC9","#f08d85","#46b566","#f1a327","#e082b2"),
               # add = "jitter",
               xlab = "",ylab="Immune inflitration",legend.title="Neoantigen Number",title = "ssGSEA")
p <-p +stat_compare_means(aes(group = group),label="p.signif")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust=1,colour = "black"),
        panel.grid = element_blank(),
        axis.text.y =element_text(colour = "black"),
        legend.position = "top",
        plot.title = element_text(hjust=0.5)
  )
p

pdf("immune_inflitration_removeCscRNA/ssGSEA_boxplot.pdf",width = 12,height=6)
p
dev.off()

tiff("immune_inflitration_removeCscRNA/ssGSEA_boxplot.tiff",width = 1200,height=600)
p
dev.off()



data1<-data1 %>% filter(celltype %in% c("Activated dendritic cell",
                                        "Immature dendritic cell",
                                        "Plasmacytoid dendritic cell",
                                        "Effector memeory CD8 T cell",
                                        "Effector memeory CD4 T cell",
                                        "T follicular helper cell",
                                        
                                        "Gamma delta T cell",
                                        "Natural killer T cell",
                                        "CD56bright natural killer cell",
                                        "CD56dim natural killer cell",
                                        "Memory B cell",
                                        "Neutrophil"
                                        
))

p <- ggboxplot(data1, x = "celltype", y = "immune_inflitration",
               color = "group", palette=c("#FFBE7A","#8ECFC9"),
               # add = "jitter",
               xlab = "",ylab="Immune inflitration",legend.title="Neoantigen Number",title = "ssGSEA")
p <-p +stat_compare_means(aes(group = group),label="p.signif")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 60,hjust = 1,vjust=1,colour = "black"),
        panel.grid = element_blank(),
        axis.text.y =element_text(colour = "black"),
        legend.position = "top",
        plot.title = element_text(hjust=0.5)
  )
p

pdf("immune_inflitration_removeCscRNA/select_ssGSEA_boxplot.pdf",width = 8,height=7)
p
dev.off()

tiff("immune_inflitration_removeCscRNA/select_ssGSEA_boxplot.tiff",width = 1200,height=800)
p
dev.off()

