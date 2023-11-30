rm(list=ls())
library(ggsci)
library(ggplot2)
library(gridExtra)
library(patchwork)

Fusion<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957_with_pheno.txt")
Fusion$id<-paste0(Fusion$fusionName,"_",Fusion$Breakpoint5,"_",Fusion$Breakpoint3)
Fusion<-Fusion %>% filter(Fusion_Type != "cscRNA")


dt<-Fusion %>% dplyr::select(id,sample) %>% distinct() %>% group_by(sample) %>% dplyr::summarise(FusionNumberPerSample=n())


# RNA_Fusion<-Fusion %>% dplyr::filter(Fusion_Type =="RNA level fusion")
# RNA_Fusion<-RNA_Fusion %>% dplyr::select(id,sample) %>% distinct()
# 
# DNA_Fusion<-Fusion %>% dplyr::filter(Fusion_Type =="DNA level fusion")
# DNA_Fusion<-DNA_Fusion %>% dplyr::select(id,sample) %>% distinct()
# ## draw the density of fusiongene number in per sample 
# dt1<-RNA_Fusion %>% group_by(sample) %>% dplyr::summarise(FusionNumberPerSample=n())
# dt2<-DNA_Fusion %>% group_by(sample) %>% dplyr::summarise(FusionNumberPerSample=n())
# dt1$Type<-"RNA level fusion"
# dt2$Type<-"DNA level fusion"
# dt<-rbind(dt1,dt2)


densFindPeak <- function(x){
  td <- density(x,na.rm = T)
  maxDens <- which.max(td$y)
  list(x=td$x[maxDens],y=td$y[maxDens])
}

limit <- dt$FusionNumberPerSample[dt$FusionNumberPerSample>=0 & dt$FusionNumberPerSample<=max(dt$FusionNumberPerSample)]

# 获得峰值
Peak = densFindPeak(limit)


library(ggsci)#+scale_color_npg()+scale_fill_npg()

# mycolor2<-pal_aaas(alpha = 0.4)(8)

p<-ggplot(dt, aes(x = FusionNumberPerSample))
p<-p+geom_density(alpha=0.4,aes(fill = "#FA7F6F",color="#FA7F6F"))+  #aes(fill = Type,color=Type),
  scale_fill_manual(values=c("#FA7F6F")) +
  scale_color_manual(values=c("#FA7F6F"))+
  geom_vline(xintercept = Peak$x, color="black",linetype="dashed",size=1)+
  # geom_vline(xintercept = normal_Peak$x, color="#8ac6d1",linetype="dashed",size=1)+
  annotate("text",x=Peak$x,y=Peak$y,label=round(Peak$x,3))+
  # annotate("text",x=normal_Peak$x,y=normal_Peak$y,label=round(normal_Peak$x,3))+
  xlab("Number of detected RNA level chimeric RNAs (per sample)")+
  ylab("Density")+theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # axis.title =  element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        legend.title = element_blank()
  )

p
pdf("perSample_fusionNumber_distribution_remove_cscRNA.pdf",width=7,height=5)
print(p)
dev.off()

