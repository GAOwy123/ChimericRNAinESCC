rm(list=ls())
setwd("D:/FusionGene2/1.statistics/20.FusionProtein")
library(tidyverse)

##========统计 MS验证的 neoantigen的比例
MS_validated_protein<-read_tsv("MS_data_validation/fusionprotein_LCMS_validation.txt")
#view(head(MS_validated_protein))
MS_validated_protein2<-MS_validated_protein  %>% separate(Protein,c("G5","B5","G3","B3"),sep="[_|-]")
MS_validated_protein2<-MS_validated_protein2 %>% filter(G5!=G3)
MS_validated_protein2$ID<-paste0(MS_validated_protein2$G5,"--",MS_validated_protein2$G3)
MS_validated_protein<-unique(MS_validated_protein2$ID)
MS_validated_protein
length(MS_validated_protein) 
#1717
#2155


all_protein<-read_delim("AGFusion/fusion_protein.combine.txt",delim="\\t",col_names = F)
all_protein2<-all_protein  %>% separate(X1,c("G5","B5","G3","B3"),sep="[_|-]")
all_protein2<-all_protein2 %>% filter(G5!=G3)
all_protein2$ID<-paste0(all_protein2$G5,"--",all_protein2$G3)
all_protein<-unique(all_protein2$ID)
length(all_protein)
# 2204
#2956

length(MS_validated_protein)/length(all_protein) ## 0.7790381



res<-tibble()
res<-rbind(res,c("Protein","validated",length(MS_validated_protein),length(MS_validated_protein)/length(all_protein)))
res<-rbind(res,c("Protein","No validated",length(all_protein)-length(MS_validated_protein),1-length(MS_validated_protein)/length(all_protein)))
colnames(res)<-c("Type","whether MS validation","Number","Ratio")
write_tsv(res,"output/Protein_MS_validated_number_ratio.txt")


###draw pie 
res<-as_tibble(res)
res$Number<-as.numeric(res$Number)
res$Ratio<-as.numeric(res$Ratio)
res$Ratio<-round(res$Ratio*100,2)


# 计算比例，保留一位小数
label <- paste(res$`whether MS validation`, "(", res$Ratio, "%)")
# 定制标签
library(RColorBrewer)
pdf("output/protein_validation_ratio_pie.pdf",width=4,height = 4)
pie(res$Number, border="white", col=c("#f39c80","#8592b5"), label=label,
    main="LC/MS")
dev.off()
