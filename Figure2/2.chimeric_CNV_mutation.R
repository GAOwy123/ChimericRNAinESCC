rm(list=ls())
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(grid)

# setwd("/Users/gaowy/Desktop/chimericRNA/FusionGene2/7.key_gene")
setwd("D:/FusionGene2/1.statistics/keygene_new")#/Users/gaowy/Desktop/work/keygene_new/

#####load 预后相关的key gene 展示  ==========================
# gene<-read_tsv("key_gene/keygene_survival_pvalue.txt")
# gene<-gene$gene

gene<-c("MSN","CREBBP","MYH11","ZFHX3","CASP8","EIF3E","PTPRC","HSPG2","BTG1","FN1","ATR")


####load mutation 
mutationData<-read_tsv("input/laml_maftools.maf")
df_mut<-mutationData %>% filter(grepl("T$",Tumor_Sample_Barcode)) %>% 
  filter(Hugo_Symbol %in% gene) %>% 
  filter(Variant_Classification %in% c( "Missense_Mutation","Nonsense_Mutation",
                                        "In_Frame_Ins","In_Frame_Del","Frame_Shift_Ins",
                                        "Frame_Shift_Del","Splice_Site","Nonstop_Mutation","Silent" )) %>% 
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification)

colnames(df_mut)<-c("gene","sample","Mutation")


####load fusion
fusionData<-read_tsv("share_data/ESCC_All_fusion_3957.changed_name.txt")

# fusionData<-fusionData %>% filter(Fusion_Type %in% c("DNA level fusion","RNA level fusion","cscRNA")) 
fusionData_1<-fusionData %>% filter(gene_5 %in% gene)  %>% dplyr::select(gene_5,sample,Fusion_Type)
fusionData_2<-fusionData %>% filter(gene_3 %in% gene)  %>% dplyr::select(gene_3,sample,Fusion_Type)
colnames(fusionData_1)<-c("gene","sample","Fusion")
colnames(fusionData_2)<-c("gene","sample","Fusion")

df_fusion<-rbind(fusionData_1,fusionData_2)
df_fusion<-unique(df_fusion)


##=======load cnv
cnv <- read_tsv("input/CDS-gene_cnv_status.xls")
cnv <- cnv[,c("Sample_id",gene)]
df_cnv<-gather(cnv,"gene","CNV",-Sample_id)
colnames(df_cnv)<-c("sample","gene","CNV")


#left_join
df_cnv_mut_rna <- df_cnv %>% full_join(df_mut,by=c("sample",'gene'))
df_cnv_mut_fusion_rna <- df_cnv_mut_rna %>% full_join(df_fusion,by=c("sample",'gene'))

####====基因，样本排序后计算比例========================================================================================================
# gene_order<-c("MAPK10","KANSL1","PDPK1", "ZFHX3", "ZNF714" , "MAST4", "ZNF680", "CREBBP" ,  "PTK2")
# gene_order<-rev(c("CREBBP","BTG1","FN1","ATR","MSN","MYH11","ZFHX3","EIF3E","CASP8","PTPRC","HSPG2"))
# gene_order<-rev(gene_order)
# df_cnv_mut_fusion_rna<-df_cnv_mut_fusion_rna %>% filter( CNV!=0 | Mutation!="NA" | Fusion !="NA") 
# 
# gene_order<-df_cnv_mut_fusion_rna %>% dplyr::select(sample,gene) %>% distinct() %>% group_by(gene) %>% 
#   dplyr::summarise(n=n()) %>% dplyr::arrange(n) %>% dplyr::select(gene) %>% unlist()


## calculate gene ratio of mutation/cnv/fusion samples
# CNV_ratio<-df_cnv_mut_fusion_rna %>% filter(CNV!=0) %>% dplyr::select(sample,gene) %>% distinct() %>% group_by(gene)  %>% dplyr::summarise(CNV_ratio=n()) %>% arrange(desc(CNV_ratio))
# Mutation_ratio<-df_cnv_mut_fusion_rna %>% filter(Mutation!="NA") %>% dplyr::select(sample,gene) %>% distinct() %>% group_by(gene)  %>% dplyr::summarise(Mutation_ratio=n()) %>% arrange(desc(Mutation_ratio))
# Fusion_ratio<-df_cnv_mut_fusion_rna %>% filter(Fusion!="NA") %>% dplyr::select(sample,gene) %>% distinct() %>% group_by(gene)  %>% dplyr::summarise(Fusion_ratio=n()) %>% arrange(desc(Fusion_ratio))


##设定不同类型的高度
df_cnv_mut_fusion<-df_cnv_mut_fusion_rna %>% dplyr::select(sample,gene,CNV,Mutation,Fusion)
df_cnv_mut_fusion<-gather(df_cnv_mut_fusion,"type","value",-c(sample,gene))
df_cnv_mut_fusion$value[is.na(df_cnv_mut_fusion$value)]<-0

# df_cnv_mut_fusion$type2<-ifelse(df_cnv_mut_fusion$value==0,"NoShow",df_cnv_mut_fusion$type)
df_cnv_mut_fusion<- df_cnv_mut_fusion %>% mutate(hgt=case_when(type=="Fusion" ~ 0.90, 
                                                               type=="CNV" ~ 0.60, 
                                                               type=="Mutation" ~ 0.30
))


# df_cnv_mut_fusion<-full_join(df_cnv_mut_fusion,CNV_ratio)
# df_cnv_mut_fusion<-full_join(df_cnv_mut_fusion,Mutation_ratio)
# df_cnv_mut_fusion<-full_join(df_cnv_mut_fusion,Fusion_ratio)
# df_cnv_mut_fusion$CNV_ratio<-df_cnv_mut_fusion$CNV_ratio/155
# df_cnv_mut_fusion$Mutation_ratio<-df_cnv_mut_fusion$Mutation_ratio/155
# df_cnv_mut_fusion$Fusion_ratio<-df_cnv_mut_fusion$Fusion_ratio/155


df_cnv_mut_fusion<-df_cnv_mut_fusion %>% filter(value != 0)
df_cnv_mut_fusion$value<-gsub("-1","Deletion",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("1","Amplification",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("Missense_Mutation","Mutation",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("Silent","Mutation",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("Nonsense_Mutation","Mutation",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("In_Frame_Del","Mutation",df_cnv_mut_fusion$value)
df_cnv_mut_fusion$value<-gsub("Frame_Shift_Ins","Mutation",df_cnv_mut_fusion$value)

df_cnv_mut_fusion$type<-factor(df_cnv_mut_fusion$type,levels=c("Fusion","CNV","Mutation"))



CREBBP_sample<-df_cnv_mut_fusion %>% filter(gene=="CREBBP") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
MSN_sample<-df_cnv_mut_fusion %>% filter(gene=="MSN") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
MYH11_sample<-df_cnv_mut_fusion %>% filter(gene=="MYH11") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
ZFHX3_sample<-df_cnv_mut_fusion %>% filter(gene=="ZFHX3") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
ATR_sample<-df_cnv_mut_fusion %>% filter(gene=="ATR") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
EIF3E_sample<-df_cnv_mut_fusion %>% filter(gene=="EIF3E") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
CASP8_sample<-df_cnv_mut_fusion %>% filter(gene=="CASP8") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
PTPRC_sample<-df_cnv_mut_fusion %>% filter(gene=="PTPRC") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
BTG1_sample<-df_cnv_mut_fusion %>% filter(gene=="BTG1") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
FN1_sample<-df_cnv_mut_fusion %>% filter(gene=="FN1") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()
HSPG2_sample<-df_cnv_mut_fusion %>% filter(gene=="HSPG2") %>% dplyr::arrange(desc(hgt)) %>% dplyr::select(sample) %>% unlist()




fusion_df<-df_cnv_mut_fusion %>% filter(type=="Fusion")
CNV_df<-df_cnv_mut_fusion %>% filter(type=="CNV")
mutation_df<-df_cnv_mut_fusion %>% filter(type=="Mutation")

dir.create("statistic_result")
write_tsv(fusion_df,"statistic_result/fusion_df.txt")
write_tsv(CNV_df,"statistic_result/CNV_df.txt")
write_tsv(mutation_df,"statistic_result/mutation_df.txt")

library(ggsci)#+scale_color_npg()+scale_fill_npg()

mycolor1<-pal_npg("nrc", alpha = 0.5)(8)
mycolor2<-pal_aaas(alpha = 0.4)(8)

# "lightblue",mycolor2[c(2,3)]
# sample_order<-c(CREBBP_sample,setdiff(gene,CREBBP_sample))
sample_order<-unique(c(CREBBP_sample,MSN_sample,MYH11_sample,
                       ATR_sample,EIF3E_sample,
                       ZFHX3_sample,FN1_sample,
                       CASP8_sample,
                       PTPRC_sample,HSPG2_sample,
                       BTG1_sample))

df_cnv_mut_fusion$sample<-factor(df_cnv_mut_fusion$sample,levels=sample_order)

write_tsv(as.data.frame(sample_order),"output/sample_order.txt")
sample_order<-read_tsv("output/sample_order.txt")
df_cnv_mut_fusion$sample<-factor(df_cnv_mut_fusion$sample,levels=sample_order$sample_order)

gene_order<-rev(c("CREBBP","MSN","MYH11","ATR","EIF3E","ZFHX3","FN1","CASP8","PTPRC","HSPG2","BTG1"))

df_cnv_mut_fusion$gene<-factor(df_cnv_mut_fusion$gene,levels=gene_order)
p = df_cnv_mut_fusion %>%
  ggplot(aes(x= sample,y=gene,fill=value,height=hgt))+
  geom_tile()+
  geom_tile(data=df_cnv_mut_fusion %>% filter(type=="Fusion"))+
  geom_tile(data=df_cnv_mut_fusion %>% filter(type=="Mutation"))+
  geom_tile(data=df_cnv_mut_fusion %>% filter(type=="CNV"))+
  scale_fill_manual(values=c("DNA level fusion"="#bc9584",#mycolor2[2]
                             "RNA level fusion"="#FA7F6F",#mycolor2[3],
                             "cscRNA"="#82B0D2",
                             "Amplification"="#FFBE7A",
                             "Deletion"="#8ECFC9",
                             "Mutation"="#BEB8DC"
                             # "Nonsense_Mutation"="#9192ab",
                             # "In_Frame_Del"="#60af46",#9E9E9E
                             # "Frame_Shift_Ins"="#efa666",
                             # # "Splice_Site"="#eddd86",
                             # "Silent"="#9987ce"
  ),na.value = "gray")+
  xlab("")+ylab("")+theme_bw()+
  theme(
    # axis.text.x=element_text(face="italic",size=12,angle=90,hjust=2,vjust=0.9),#color="black",
    axis.text.x=element_blank(),
    axis.text.y = element_text(face="italic",color="black",size=18),
    axis.ticks = element_blank(),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    panel.border = element_rect(),
    # panel.background = element_rect(fill = "gray100"), ##背景颜色#f0f8ff
    #      legend.text=element_text(size=11),legend.spacing.x = unit(0.3, 'cm'),
    legend.position = "bottom",
    legend.title=element_blank())+
  theme(plot.title = element_text(hjust = 0.5,size=20))
p

pdf("output/keygene/keygene_fusion_mutation_CNV.pdf",width=7,height=4)
print(p)
dev.off()

