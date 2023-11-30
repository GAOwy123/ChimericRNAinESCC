rm(list=ls())
setwd("D:/FusionGene2/1.statistics/10.CTSC-RAB38")
library(tidyverse)
library(reshape2)

library(ggpubr)

pheno<-read_tsv("D:/work/ESCC/runtime/0.rawData/All_156_key_phenotype_20200723.txt")
df<-pheno
df<-df[,c("Sample_id","Status","Survival_time","Overage","DFSurvival_time","Gender","Age" ,"grade","TNM",
          "recurrence_status","Smoking_history","Drinking_history", "Location", "Tumor_size")]

variates<-c("Sample_id","Survival_time","Gender","Age","Status","grade","TNM", "recurrence_status","Smoking_history","Drinking_history","Location","Tumor_size")
df<-df[,variates]
df$Age<-factor(ifelse(df$Age>=60,">=60","<60"),levels = c("<60",">=60"))
df$Tumor_size<-cut(as.numeric(df$Tumor_size),3)
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


Fusion<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957_with_pheno.txt")
Fusion <- Fusion %>% filter(fusionName=="CTSC--RAB38")
fusionsample<-Fusion %>% dplyr::select(sample) %>% unlist() %>% unique()

df$Fusion<-ifelse(df$Sample_id %in% fusionsample,"Chimera","No Chimera")
df$Fusion<-factor(df$Fusion,levels=c("No Chimera","Chimera"))

clinical_state<-colnames(df)[2:11]
dir.create("pheno")
for(state in clinical_state){
  # state<-clinical_state[5]
  sub_df<-df[,c("Fusion",state)]
  sub_df<-sub_df %>% drop_na()
  
  genename<-"Fusion"
  gene_smoke<-table(filter(sub_df[,c(genename,state)],sub_df[,c(genename,state)][,2]!="loss"))
  # nrow<-length(unique(df[,genename]))
  # ncol<-length(unique(df[,state][!is.na(df[,state])]))
  nrow<-nrow(gene_smoke)
  ncol<-ncol(gene_smoke)
  df_smoke = as.data.frame(prop.table(gene_smoke,1))
  pvalue <- fisher.test(matrix(gene_smoke,ncol=ncol,nrow=nrow),simulate.p.value=TRUE) #
  
  write.table(cbind(genename,state,round(pvalue$p.value,5)),"pheno/gene_pheno_pvalue.txt",sep="\t",append = T,quote=F,row.names = F,col.names = F)
  p_smoke =ggplot(data=df_smoke, aes(x=df_smoke[,genename], y=Freq, fill=df_smoke[,state])) +
    geom_bar(stat="identity") +
    # geom_text(aes(label=paste0(round(Freq*100,1),"%")))+
    geom_text(aes(label=paste0(round(Freq*100,1),"%")),
              size=5,position = position_stack(vjust = 0.5))+ 
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0,color="black",size=8,hjust=0.5,lineheight=1))+
    theme(axis.text.y = element_text(color="black",size=8,hjust=1,lineheight=1))+
    labs(title=paste(state,paste0("p_value = ",round(pvalue$p.value,5)),sep="\n"),y="Frequence",x="CSS Class")+
    theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(),axis.line = element_line(colour = "black"),legend.title=element_blank())
  
  ggsave(paste("pheno/",genename,"_",state,".pdf",sep=""), width=6, height=4.5)
  
}


###=========================survival ========================

df_sur<-pheno
df_sur$Fusion<-ifelse(df_sur$Sample_id %in% fusionsample,"Chimera","No Chimera")
df_sur$Fusion<-factor(df_sur$Fusion,levels=c("No Chimera","Chimera"))

df_sur$Status<-as.numeric(df_sur$Status)
df_sur$Survival_time<-as.numeric(df_sur$Survival_time)
df_sur$Survival_time<-df_sur$Survival_time/30


  
fit <- survfit(Surv(Survival_time, Status) ~ Fusion, data = df_sur) #df_sur[,4]
diff<- survdiff(Surv(Survival_time, Status) ~ Fusion, data = df_sur)
pvalue <-1-pchisq(diff$chisq, length(diff$n)-1)


pdf("CTSC-RAB38_survival.pdf", width=7, height=7)
  
print(ggsurvplot(fit,data = df_sur,pval = T,legend.title = "CTSC--RAB38",xlab="Time(Months)",
                   legend = c(0.7,0.9) ,# 指定图例位置
                   # legend.labs=c("High","Low"),
                   palette = c("lightblue","pink"),
                   risk.table = TRUE, # 添加风险表
                   fontsize=6,
                   #surv.median.line = "hv",  # 增加中位生存时间
                   #conf.int = TRUE # 增加置信区间
                   pval.size=10, #pvalue 字的大小
                   font.x=c(14,"plain","black"), 
                   font.y=c(15,"plain","black"),
                   font.tickslab=c(15,"plain","black"),
                   font.main = c(15, "plain", "black"),
                   font.legend=c(15,"plain","black")
  ))
  
dev.off()
  
