rm(list = ls())
setwd("/Users/gaowy/Desktop/work/keygene_new/")
library(tidyverse)
library(readxl)
library(hash)

phenoData<-read_tsv("input/ESCC_155_phenotype_20210512.txt")

ref_data<-read_tsv("input/Homo_sapiens.gene_info")
ref_data<-ref_data %>% dplyr::select(Symbol,Synonyms)

fusionData<-read_tsv("share_data/ESCC_All_fusion_3957.changed_name.txt")
#fusionData<-fusionData %>% filter(Fusion_Type %in% c("DNA level fusion","RNA level fusion"))


#####load 预后相关的key gene 展示  ==========================
# gene<-read_tsv("output/drivergene_survival_pvalue.txt")
# # gene<-read_tsv("key_gene/keygene_survival_pvalue.txt")
# gene<-gene %>% filter(survival_pvalue<0.05)
# gene<-gene$gene

gene<-c("MSN","CREBBP","MYH11","ZFHX3","CASP8","EIF3E","PTPRC","HSPG2","BTG1","FN1","ATR")

changeName<-function(ref_data,query_name){
  # find_gene<-intersect(query_name,ref_data$Symbol)
  tmp<-data.frame(old_gene=query_name)
  query_name<-setdiff(query_name,ref_data$Symbol)
  library(hash)
  h<-hash()
  for(i in 1:length(ref_data$Synonyms)){
    old.gene<-unlist(strsplit(ref_data$Synonyms[i],"|",fixed = T))
    for (g in old.gene){
      if(g %in% query_name){
        h[[g]]<-ref_data$Symbol[i]
      }
    }
  }
  res<-data.frame(old_gene=keys(h),new_gene=values(h))
  
  tmp<-left_join(tmp,res)
  tmp <- tmp %>% mutate(final_gene=ifelse(is.na(new_gene),old_gene,new_gene))
  return(tmp$final_gene)
}

#
#
# expData<-read_tsv("D:/work/ESCC/runtime/0.rawData/CDS-gene_155.paired.all_331_FPKM.txt")
# expData<-read_tsv("input/escc.phaseII.counts.featureCounts.TPM.geneSymbol.txt")
expData<-read_tsv("input/escc.phaseII.hisat2.genesymbol.median.log2.TPM+1.txt")
# ZNF714_Expr<-expData %>% dplyr::filter(gene=="ZNF714")
# expData<-read_tsv("input/escc.phaseII.counts.featureCounts.TPM.geneSymbol.txt")

colnames(expData)<-c("Group.1",colnames(expData)[2:ncol(expData)])
expData2<-expData[,c("Group.1",colnames(expData)[grepl("T",colnames(expData))])]
expData2<-expData2 %>% dplyr::filter(Group.1 %in% gene)
# Gene<-changeName(ref_data,expData2$Gene)
# expData2$Gene<-Gene
# 
# gene<-intersect(expData2$Gene,gene)
# expData2<-expData2 %>% dplyr::filter(Gene %in% gene)

expr<-expData2 %>% gather("sample","expr",-Group.1)


dir.create("output/keygene/")

res<-tibble()
for(i in gene){
  # i<-gene[1]
  sample_1<-fusionData %>% filter(gene_5 %in% i) %>% dplyr::select(sample)
  sample_2<-fusionData %>% filter(gene_3 %in% i) %>% dplyr::select(sample)
  sample<-rbind(sample_1,sample_2)
  fusionsample<-unique(sample$sample)
  
  sub<-expr %>% dplyr::filter(Group.1 ==i ) %>% mutate(Fusion_Status=ifelse(sample %in% fusionsample,"Chimera","No Chimera"))
  sub$Fusion_Status<-factor(sub$Fusion_Status,levels = c("No Chimera","Chimera"))
  res<-rbind(res,sub)
  
}

write_tsv(res,"output/key_genegene_expr.txt")

##===============draw the 9 gene expression in Chimera and No chimera ====================


# res$expr<-log2(res$expr+1)
gene_order<-unique(res$Group.1)
res$Gene<-factor(res$Group.1,levels = gene)
p<-ggplot(res,aes(x=Fusion_Status,y=expr,color=Gene))+
  geom_violin(draw_quantiles=0.5)+geom_jitter(size=.4)+
  facet_grid(. ~ Gene)+
  # stat_summary(fun=median, geom="line", color="black")+
  scale_color_manual(values=c("#bad5e8","#307bb5","#d7ee91","#2ea41c",
                              "#fea8a9","#d82818","#ffc875","#f17b00","#e3c2f1",
                                "#b58db6","#f6bd60"))+
  xlab("")+ylab("Gene expression TPM (log2)")+
  theme_bw()+
  theme(axis.text = element_text(color="black"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(angle = 60,hjust = 1))

p
pdf("output/keygene/keygene_expr_violin_plot2.pdf",width=10,height = 4)
print(p)
dev.off()


library(ggpubr)
pdf("output/key_gene9_keygene_expr_violin_plot_pvalue2.pdf",width=7,height = 7)

p2<-ggviolin(res,x="Fusion_Status",y="expr",facet.by = "Gene",
             color="Gene",ncol=length(gene),
             # group="Fusion_Status",
             add = "jitter",#,
             palette = c("#bad5e8","#307bb5","#d7ee91","#2ea41c",
             "#fea8a9","#d82818","#ffc875","#f17b00","#e3c2f1",
             "#b58db6","#f6bd60")
             )+
  stat_compare_means()
print(p2)
dev.off()

tmp1<-res %>% group_by(Gene,Fusion_Status) %>% summarise(Med=median(expr))
tmp2<-res %>% group_by(Gene,Fusion_Status) %>% summarise(Min=min(expr))
tmp3<-res %>% group_by(Gene,Fusion_Status) %>% summarise(Max=max(expr))

res2<-inner_join(tmp1,tmp2)
res2<-inner_join(res2,tmp3)

as.data.frame(res2)
write_tsv(res2,"output/keygene/keygene_expression_range.txt")
# Gene Fusion_Status       Med      Min       Max
# 1     MSN    No Chimera  8.244109 5.779211  9.628863
# 2     MSN       Chimera  8.579657 7.474061  9.366737
# 3  CREBBP    No Chimera  4.606851 3.223885  5.798645
# 4  CREBBP       Chimera  4.320135 3.962183  5.156303
# 5   MYH11    No Chimera  4.214694 1.541213  9.279206
# 6   MYH11       Chimera 11.455718 8.561050 12.071107
# 7   ZFHX3    No Chimera  2.933024 1.245802  4.512790
# 8   ZFHX3       Chimera  2.090416 1.390760  2.731817
# 9   CASP8    No Chimera  4.515193 3.100247  5.692142
# 10  CASP8       Chimera  3.923126 3.233518  4.612735
# 11  EIF3E    No Chimera  5.064547 4.063072  6.421065
# 12  EIF3E       Chimera  5.305540 4.795028  5.991621
# 13  PTPRC    No Chimera  5.500119 2.200652  8.491906
# 14  PTPRC       Chimera  6.153512 5.452335  6.854689
# 15  HSPG2    No Chimera  5.932496 3.833270  8.066491
# 16  HSPG2       Chimera  5.088737 3.893241  7.504799
# 17   BTG1    No Chimera  6.564663 5.062028  9.585925
# 18   BTG1       Chimera  7.141021 5.787140  9.370384
# 19    FN1    No Chimera  9.137790 5.355699 12.434420
# 20    FN1       Chimera 10.971179 9.053345 13.238352
# 21    ATR    No Chimera  4.678999 1.933833  5.971263
# 22    ATR       Chimera  4.942261 4.706234  5.178287



####calculate the significance =============
dt<-read_tsv("output/key_genegene_expr.txt")
dt$Fusion_Status<-factor(dt$Fusion_Status,levels=c("No Chimera","Chimera"))

# gene_order<-c("PTK2","CREBBP", "MAST4", "PDPK1","ZNF680", 
#               "ZNF714", "ZFHX3",    "KANSL1",
#               "MAPK10")
dt$Gene<-factor(dt$Group.1,levels = gene_order)

gene_wilcox_pvale<-tibble()
for(i in gene_order){
  # i<-"ZNF680"
  sub<-dt %>% dplyr::filter(Group.1==i)
  res<-wilcox.test(expr~Fusion_Status,sub)
  gene_wilcox_pvale<-rbind(gene_wilcox_pvale,c(i,res$p.value))
}
colnames(gene_wilcox_pvale)<-c("gene","wilcox_pvale")
gene_wilcox_pvale$pvalue<-round(as.numeric(gene_wilcox_pvale$wilcox_pvale),2)
write_tsv(gene_wilcox_pvale,"output/keygene/keygene_expr_wilcox_pvalue.txt")
