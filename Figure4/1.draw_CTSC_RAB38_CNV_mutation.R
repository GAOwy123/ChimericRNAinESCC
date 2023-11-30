rm(list=ls())
library(tidyverse)
setwd("D:/FusionGene/merge.fusion/fusiongene.phenotype")

# expr<-read_tsv("D:/FusionGene/data/expression/escc.phaseII.counts.featureCounts.TPM.geneSymbol.txt")
# expr<-expr[,c(1,grep("*T$",colnames(expr)))]

dir.create("output/Fusion.singlegene")
library(tidyverse)
library(dplyr)
library(ggplot2)
library(reshape2)
library(survival)
library(survminer)

###draw CTSC-RAB38 oncoplot

genelist<-c("CTSC","RAB38")

pheno<-read_tsv("D:/work/ESCC/runtime/0.rawData/All_156_key_phenotype_20200723.txt")


##==========load Fusion 
Fusion<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957_with_pheno.txt")
Fusion <- Fusion %>% filter(fusionName=="CTSC--RAB38")
fusionsample<-Fusion %>% dplyr::select(sample) %>% unlist() %>% unique()
df_fusion<-tibble(sample=fusionsample, gene="CTSC--RAB38",Fusion="Chimera")

##==========load mutation
####load mutation 
# mutationData<-read_tsv("D:/ESCC/final.data/maf/maf.ESCC/laml_maftools.maf")
mutationData<-read_tsv("D:/sharedata/mutation/155ESCC.maf")

df_mut<-mutationData %>% filter(grepl("T$",Tumor_Sample_Barcode)) %>% 
  filter(Hugo_Symbol %in% genelist) %>% 
  filter(Variant_Classification %in% c( "Missense_Mutation","Nonsense_Mutation",
                                        "In_Frame_Ins","In_Frame_Del","Frame_Shift_Ins",
                                        "Frame_Shift_Del","Splice_Site","Nonstop_Mutation","Silent" )) %>% 
  dplyr::select(Hugo_Symbol,Tumor_Sample_Barcode,Variant_Classification)

colnames(df_mut)<-c("gene","sample","Mutation")

##=======load cnv
cnv <- read_tsv("D:/FusionGene2/7.key_gene/input/CDS-gene_cnv_status.xls")
cnv <- cnv[,c("Sample_id",genelist)]
df_cnv<-gather(cnv,"gene","CNV",-Sample_id)
colnames(df_cnv)<-c("sample","gene","CNV")
df_cnv<-df_cnv %>% filter(CNV !=0)


########Merge 
#left_join
df_cnv_mut <- df_cnv %>% full_join(df_mut,by=c("sample",'gene'))
df_cnv_mut_fusion <- df_cnv_mut %>% full_join(df_fusion,by=c("sample",'gene'))

df_cnv_mut_fusion$CNV<-gsub("-1","Deletion",df_cnv_mut_fusion$CNV)
df_cnv_mut_fusion$CNV<-gsub("1","Amplication",df_cnv_mut_fusion$CNV)
dt<-gather(df_cnv_mut_fusion,"Type","Value",-c(sample,gene))
dt<-dt %>% filter(!is.na(Value))

samples<-pheno %>% dplyr::select(Sample_id)
dt<-left_join(samples,dt,by=c("Sample_id"="sample"))
colnames(dt)<-c("sample","gene","Type","Value")

intersect(df_cnv$sample,df_fusion$sample)
intersect(df_mut$sample,df_fusion$sample)
intersect(df_mut$sample,df_cnv$sample)

# sample<-c(df_fusion$sample,setdiff(df_cnv$sample,df_fusion$sample),df_mut$sample)
# sample_order<-c(sample,setdiff(dt$sample,sample))


#============draw onco (recurrent Fusion gene across patients) ====================
library(pheatmap)
library(ggplot2) 


dt$Value<-factor(dt$Value,levels = c("Chimera","Amplication",
                                     "Deletion","Missense_Mutation","NA"))
# sample_order<-dt %>% dplyr::select(sample,Value) %>% group_by(sample,Value) %>% 
#               summarise(count=n()) %>% dplyr::arrange(Value,desc(count)) %>% 
#               dplyr::select(sample) %>% unlist() %>% unique()

# dt$Value[is.na(dt$Value)]<-"Wild"
dt$gene<-ifelse(is.na(dt$Value),"CTSC--RAB38",dt$gene)
dt$gene<-factor(dt$gene,levels =c("CTSC","RAB38","CTSC--RAB38"))
sample_order<-c(setdiff(fusionsample,"167T"),"167T",
                "039T","208T","029T","138T")


dt$sample<-factor(dt$sample,levels=c(sample_order,setdiff(dt$sample,sample_order)))




onco <-ggplot(dt, aes(sample, gene)) + 
  geom_tile(aes(fill = Value),colour = "gray")+xlab("")+ylab("")

cols <- c("Amplication" = "#FFBE7A", "Deletion" = "#8ECFC9", "Missense_Mutation" = "#BEB8DC", "Chimera" = "#FA7F6F")
onco<-onco + scale_fill_manual(values=cols) + 
      theme_bw()+
      theme(axis.text.y= element_text(size = 18),
        # axis.text.x = element_text(size =18,vjust = 0.9, hjust = 1, angle = 45),
        axis.text.x = element_blank(),
        axis.ticks.x= element_blank(),
        panel.grid.major=element_line(colour=NA),
        plot.margin=unit(c(0,0.5,0.5,0),'lines')
  )+theme(legend.position = "bottom")
# +geom_text(aes(sample, fusionName, label = ReadCount), color = "black", size = 5)


onco



####
##======== Add clinical information  ======================
##############################################################################
#########              phenotype annotation                         ##########
##############################################################################
df<-pheno
df<-df[,c("Sample_id","Status","Survival_time","Overage","DFSurvival_time","Gender","Age" ,"grade","TNM",
          "recurrence_status","Smoking_history","Drinking_history", "Location", "Tumor_size")]

variates<-c("Sample_id","Gender","Age","Status","Survival_time","grade","TNM", "recurrence_status","Smoking_history","Drinking_history","Location")
df<-df[,variates]
df$Age<-factor(ifelse(df$Age>60,">60","<=60"),levels=c("<=60",">60"))
# df$Tumor_size<-cut(as.numeric(df$Tumor_size),3)
df$Gender<-factor(df$Gender,levels=c("Female","Male"))
df$Status<-gsub("0","Alive",df$Status)
df$Status<-gsub("1","Dead",df$Status)
df$Status<-factor(df$Status)
df$grade<-factor(df$grade,levels = c("G1","G2","G3"))
df$TNM<-factor(df$TNM,levels = c("I","II","III","IV"))
df$recurrence_status<-gsub("0","No",df$recurrence_status)
df$recurrence_status<-gsub("1","Yes",df$recurrence_status)
df$recurrence_status<-factor(df$recurrence_status)
df$Smoking_history<-factor(df$Smoking_history,levels=c("never","light","moderate","heavy"))
df$Drinking_history<-factor(df$Drinking_history,levels = c("never","light","moderate","heavy"))
df$Location<-factor(df$Location,c("Lower","Middle","Upper"))
df$Survival_time<-ifelse(df$Survival_time<365,"<1 year",ifelse(df$Survival_time<365*2,"<2 year",ifelse(df$Survival_time<365*3,"<3 year",">=3 year")))
df$Survival_time<-factor(df$Survival_time,levels = c("<1 year","<2 year","<3 year",">=3 year"))


# df<-df %>% filter(Sample_id %in% sample_order)

df$Sample_id<-factor(df$Sample_id,levels = levels(dt$sample))

df2<-gather(df,"pheno","Value",-Sample_id)
df2$pheno<-factor(df2$pheno,levels=c("Survival_time","Status","Gender","Age","grade",
                                    "TNM","recurrence_status",
                                    "Location","Smoking_history",
                                    "Drinking_history"))

# df3<-df2 %>% filter(pheno !="Survival_time")
onco2 <-ggplot(df2, aes(Sample_id, pheno)) + 
  geom_tile(aes(fill = Value),colour = "gray")+xlab("")+ylab("")


cols2 <- c("<1 year" = "#03045e","<2 year"="#0077b6","<3 year"="#00b4d8",">=3 year"="#90e0ef",
           "Male" = "#a26769", "Female" = "#d5b9b2",
           "<=60" ="#cdafe4",">60"="#a06cd5",
           # "never"="#ffd39b","light"="#eec591","moderate"="#cdaa7d","heavy"="#8b7355",
           "never"="#f6efe2","light"="#f7ceb5","moderate"="#ef8f59","heavy"="#e05f10",
           "Yes"="#e69689","No"="#93bde5",
           "Alive"="#e3edd6", "Dead"="#f3c0ad",
           "I"="#f6c0b9", "II"="#ef9f95", "III"="#e78072", "IV"="#e26150",
           "Lower"="#ecd1d5", "Middle"="#cabed0", "Upper"="#8c97ac",
           "G1"="#943c85","G2"="#c282a5","G3"="#e4b7be"
           )
onco2<-onco2 + 
  scale_fill_manual(values=cols2) +
  # scale_fill_gradient(low = "blue", high = "red")+
  theme_bw()+
  theme(axis.text.y= element_text(size = 18),
        # axis.text.x = element_text(size =18,vjust = 0.9, hjust = 1, angle = 45),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major=element_line(colour=NA),
        plot.margin=unit(c(0,0.5,0.5,0),'lines')
  )+theme(legend.position = "bottom")
# +geom_text(aes(sample, fusionName, label = ReadCount), color = "black", size = 5)
onco2




library(cowplot)
pdf("D:/FusionGene2/1.statistics/10.CTSC-RAB38/CTSC_RAB38_CNV_mutation_new.pdf",width = 9,height = 5)
plot_grid(onco2,onco,align=c("v"),ncol=1, rel_heights = c(2.5, 1))
dev.off()

