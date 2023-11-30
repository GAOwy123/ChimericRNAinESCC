# 展示，kinase fusion event在所有fusion的比例有多少，

rm(list=ls())
library(ggpubr)
library(tidyverse)
setwd("D:/FusionGene2/1.statistics/20.FusionProtein/pfam_domain_prediction")
# dir.create("output/kinase")

allfusionevent<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/allfusioneventID.txt",col_names = F)
allfusionevent<-allfusionevent %>% separate(X1,c("G5","B5","G3","B3"),sep="-|_",remove = F)
allfusionevent<-allfusionevent %>% filter(G5!=G3)

kinasefusionevent<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/kinase_fusion_events_ID.txt",col_names = F)
kinasefusionevent<-kinasefusionevent %>% separate(X1,c("G5","B5","G3","B3"),sep="-|_",remove = F)
kinasefusionevent<-kinasefusionevent %>% filter(G5!=G3)

dir.create("output_remove_cscRNA")
kinasefusionevent$id<-paste0(kinasefusionevent$G5,"--",kinasefusionevent$G3,"_",kinasefusionevent$B5,"_",kinasefusionevent$B3)
write_tsv(kinasefusionevent,"output_remove_cscRNA/kinasefusionevent_remove_cscRNA.txt",col_names = F)
write_tsv(allfusionevent,"output_remove_cscRNA/allfusionevent_remove_cscRNA.txt")


length(setdiff(allfusionevent$X1,kinasefusionevent$X1)) #5929
length(intersect(allfusionevent$X1,kinasefusionevent$X1)) #430

kinase_percentage<-data.frame(kinase=c("kinase","non-kinase"),Number=c(length(intersect(allfusionevent$X1,kinasefusionevent$X1)),length(setdiff(allfusionevent$X1,kinasefusionevent$X1))))

kinase_percentage$Percentage<-round(kinase_percentage$Number/sum(kinase_percentage$Number)*100,digits = 2)
kinase_percentage$kinase<-factor(kinase_percentage$kinase,levels=c("kinase","non-kinase"))



###=======kinse fusion events 中能形成protein的比例是多少
allfusionevent<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/allfusioneventID.txt",col_names = F)
allfusionevent<-allfusionevent %>% separate(X1,c("G5","B5","G3","B3"),sep="-|_",remove = F)
allfusionevent<-allfusionevent %>% filter(G5!=G3)

allprotein<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/fusion_protein.combine2.txt")
allfusionevent<-left_join(allfusionevent,allprotein,by=c("X1"="event"))
kinasefusionevent<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/kinase_fusion_events_ID.txt",col_names = F)
kinasefusionevent<-kinasefusionevent %>% separate(X1,c("G5","B5","G3","B3"),sep="-|_",remove = F)
kinasefusionevent<-kinasefusionevent %>% filter(G5!=G3)

kinasefusionevent<-allfusionevent %>% filter(X1 %in% kinasefusionevent$X1)
kinasefusionevent<-kinasefusionevent %>% mutate(coding=ifelse(is.na(ENSP),"No","Yes"))

coding_percentage<-kinasefusionevent %>% dplyr::select(X1,coding) %>% group_by(coding) %>% summarise(Number=n())
coding_percentage$Percentage<-round(coding_percentage$Number/sum(coding_percentage$Number)*100,digits = 2)
coding_percentage$coding<-factor(coding_percentage$coding,levels=c("Yes","No"))

###============kinase protein 的kinase domain是否完整============================
kinase<-read_tsv("Figure6_new/kinase_chimericevent_pkinasedomain_if_intact.txt")
kinase<-kinase %>% dplyr::select(ID,Pkinase_domain)
colnames(kinase)<-c("event","domain_if_intact")
kinase<-kinase %>% dplyr::select(event,domain_if_intact) %>% distinct()
domain_if_inatact_percentage<-kinase %>% dplyr::select(event,domain_if_intact) %>% group_by(domain_if_intact) %>% summarise(Number=n())
domain_if_inatact_percentage$Percentage<-round(domain_if_inatact_percentage$Number/sum(domain_if_inatact_percentage$Number)*100,digits = 2)

##==================draw figure
kinase_percentage

coding_percentage

domain_if_inatact_percentage

kinase_percentage$Type<-"Kinase"
coding_percentage$Type<-"Protein coding ability"
domain_if_inatact_percentage$Type<-"Kinase domain"

colnames(kinase_percentage)<-c("subtype","Number","Percentage","Type")
colnames(coding_percentage)<-c("subtype","Number","Percentage","Type")
colnames(domain_if_inatact_percentage)<-c("subtype","Number","Percentage","Type")


dt<-rbind(kinase_percentage,coding_percentage)
dt<-rbind(dt,domain_if_inatact_percentage)

dt$Type<-factor(dt$Type,levels=c("Kinase","Protein coding ability","Kinase domain"))

cols<-list("kinase"="#f9999b","non-kinase"="#bebada",
           "Yes"="#fdb462","No"="#619cff",
           "Disrupted"="#ee8265","Intact"="#2fa61d")
p1<-ggplot(data = dt, aes(fill = subtype ,x = Type, y=Percentage)) +
  geom_bar(stat="identity",width = 0.6)+
  geom_text(aes(label=paste0(round(Percentage,1),"%")),
            size=5,position = position_stack(vjust = 0.5))+
  # scale_fill_brewer(palette = 'Set3')+
  scale_fill_manual(values = cols)+
  labs(x = "", y = "Frequency", title = "") +theme_bw()+
  theme(panel.grid=element_blank(),axis.text = element_text(color="black"))

p1

pdf("Figure6_new/kinase_barplot_new.pdf",width=6,height = 5)
print(p1)
dev.off()
tiff("Figure6_new/kinase_barplot_new.tiff",width=600,height = 500)
print(p1)
dev.off()



