rm(list = ls())
setwd("/Users/gaowy/Desktop/work/keygene_new/")
gene<-read_tsv("output/drivergene_survival_pvalue.txt")
# gene<-read_tsv("key_gene/keygene_survival_pvalue.txt")
gene<-gene %>% filter(survival_pvalue<0.05)
gene<-gene$gene

phenoData<-read_tsv("input/ESCC_155_phenotype_20210512.txt")


ref_data<-read_tsv("input/Homo_sapiens.gene_info")
ref_data<-ref_data %>% select(Symbol,Synonyms)
fusionData<-read_tsv("share_data/ESCC_All_fusion_3957.changed_name.txt")

############============================================================================
############===================drivergene   ======================================
############============================================================================
flag<-0
res<-tibble()
splots <- list()
survival_pvalue<-tibble()
for(i in gene){
  #i<-gene[1]
  # 
  flag<-flag+1
  df<-fusionData %>% dplyr::select(sample,Status,Survival_time) %>% unique()
  
  sample_1<-fusionData %>% filter(gene_5 %in% i) %>% dplyr::select(sample)
  sample_2<-fusionData %>% filter(gene_3 %in% i) %>% dplyr::select(sample)
  fusion_sample<-rbind(sample_1,sample_2)
  fusion_sample<-unique(fusion_sample$sample)
  
  df<-df %>% mutate(Fusion_Status=ifelse(sample %in% fusion_sample,"Chimera","No Chimera"))
  
  df$Survival_time<-as.numeric(df$Survival_time)/30
  df <- filter(df,df$Status!="loss") 
  df$Status <- as.numeric(df$Status)
  
  if(length(unique(df$Fusion_Status))==2 & sum(df$Fusion_Status=="Chimera")>1){
    fit <- survfit(Surv(Survival_time, Status) ~ Fusion_Status, data = df)
    diff<- survdiff(Surv(Survival_time, Status) ~ Fusion_Status, data = df)
    pvalue <-1-pchisq(diff$chisq, length(diff$n)-1)
    
    splots[[flag]]<-ggsurvplot(fit,data = df,pval = T,legend.title = paste0(i),xlab="Time(Months)",
                               legend = c(0.7,0.9) ,# 指定图例位置
                               legend.labs=c("Chimera","No Chimera"),
                               palette = c("pink", "blue"),
                               risk.table = TRUE, # 添加风险表
                               fontsize=6,
                               #surv.median.line = "hv",  # 增加中位生存时间
                               #conf.int = TRUE # 增加置信区间
                               pval.size=10, #pvalue 字的大小
                               font.x=c(14,"plain","black"), 
                               font.y=c(15,"plain","black"),
                               font.tickslab=c(15,"plain","black"),
                               font.main = c(18, "plain", "black"),
                               font.legend=c(15,"plain","black")
    )
    # assign(paste("fig_sur",flag,sep=""),fig_sur)
    
    survival_pvalue<-rbind(survival_pvalue,c(i,pvalue))
  }
}

dir.create("output/keygene")
colnames(survival_pvalue)<-c("gene","survival_pvalue")
write_tsv(survival_pvalue,"output/keygene/survival_pvalue.txt")
figure <- arrange_ggsurvplots(splots, ncol=4,print = FALSE)
ggsave(paste0("output/keygene/","drivergene","_survival",".pdf",sep=""), figure,width=20,height=6)

###+====================最终showed gene

# removed_gene<-c("ZNF680","HSP90AB1","RNF213")
# selected_gene<-setdiff(gene,removed_gene)
selected_gene<-c("MSN","CREBBP","MYH11","ZFHX3","CASP8","EIF3E","PTPRC","HSPG2","BTG1","FN1","ATR")
flag<-0
res<-tibble()
splots <- list()
survival_pvalue<-tibble()
for(i in selected_gene){
  #i<-gene[1]
  # 
  flag<-flag+1
  df<-fusionData %>% dplyr::select(sample,Status,Survival_time) %>% unique()
  
  sample_1<-fusionData %>% filter(gene_5 %in% i) %>% dplyr::select(sample)
  sample_2<-fusionData %>% filter(gene_3 %in% i) %>% dplyr::select(sample)
  fusion_sample<-rbind(sample_1,sample_2)
  fusion_sample<-unique(fusion_sample$sample)
  
  df<-df %>% mutate(Fusion_Status=ifelse(sample %in% fusion_sample,"Chimera","No Chimera"))
  
  df$Survival_time<-as.numeric(df$Survival_time)/30
  df <- filter(df,df$Status!="loss") 
  df$Status <- as.numeric(df$Status)
  
  if(length(unique(df$Fusion_Status))==2 & sum(df$Fusion_Status=="Chimera")>1){
    fit <- survfit(Surv(Survival_time, Status) ~ Fusion_Status, data = df)
    diff<- survdiff(Surv(Survival_time, Status) ~ Fusion_Status, data = df)
    pvalue <-1-pchisq(diff$chisq, length(diff$n)-1)
    
    splots[[flag]]<-ggsurvplot(fit,data = df,pval = T,legend.title = paste0(i),xlab="Time(Months)",
                               legend = c(0.7,0.9) ,# 指定图例位置
                               legend.labs=c("Chimera","No Chimera"),
                               palette = c("pink", "blue"),
                               risk.table = TRUE, # 添加风险表
                               
                               fontsize=6,
                               #surv.median.line = "hv",  # 增加中位生存时间
                               #conf.int = TRUE # 增加置信区间
                               
                               pval.size=10, #pvalue 字的大小
                               # font.x=c(14,"plain","black"), 
                               # font.y=c(15,"plain","black"),
                               # font.tickslab=c(15,"plain","black"),
                               # font.main = c(18, "plain", "black"),
                               # font.risk.table.legend.title=c(15,"italic","black"),
                               font.legend.title=c(15,"plain","black")
                               
                            
    )
    # assign(paste("fig_sur",flag,sep=""),fig_sur)
    
    survival_pvalue<-rbind(survival_pvalue,c(i,pvalue))
  }
}


colnames(survival_pvalue)<-c("gene","survival_pvalue")
write_tsv(survival_pvalue,"output/keygene/survival_pvalue_final.txt")
figure <- arrange_ggsurvplots(splots, ncol=4,print = FALSE)
ggsave(paste0("output/keygene/","drivergene","_survival_final",".pdf",sep=""), figure,width=20,height=5.5)

