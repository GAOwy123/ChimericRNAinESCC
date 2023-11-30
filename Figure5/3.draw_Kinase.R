rm(list=ls())
setwd("D:/FusionGene2/1.statistics/20.FusionProtein/pfam_domain_prediction")
library(tidyverse)
kinase_fusion_domain<-read_tsv("kinase_fusiongene.pfam",comment = "#",col_names = F)
kinase_fusion_domain$X1<-gsub(" +",":",kinase_fusion_domain$X1,perl = TRUE)

# kinase_fusion_domain<-as.data.frame(kinase_fusion_domain)
kinase_fusion_domain2<-kinase_fusion_domain %>% separate(X1,c("seq id","alignment start","alignment end",
                                                              "envelope start", "envelope end",
                                                              "hmm acc", "hmm name", "type",
                                                              "hmm start", "hmm end", "hmm length" ,"bit score", 
                                                              "E-value", "significance" ,"clan"),
                                                         remove = T,sep=":")

kinase_fusion_domain3<-kinase_fusion_domain2 %>% separate(`seq id`,c("G5","B5","G3","B3"),sep="-|_",remove = F)


kinase_domain<-read_tsv("kinase.pfam",comment = "#",col_names = F)
kinase_domain$X1<-gsub(" +",":",kinase_domain$X1,perl = TRUE)
kinase_domain2<-kinase_domain %>% separate(X1,c("seq id","alignment start","alignment end",
                                                "envelope start", "envelope end",
                                                "hmm acc", "hmm name", "type",
                                                "hmm start", "hmm end", "hmm length" ,"bit score", 
                                                "E-value", "significance" ,"clan"),
                                           remove = T,sep=":")

kinase_domain3<-kinase_domain2 %>% separate(`seq id`,c("tmp1","uniprotID","kinase"),remove = F)




tmp<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/output_remove_cscRNA/uniprot-download_true_format_fasta_query__28existence_3A1_29-2023.05.22-14.20.38.81.fasta",col_names = F)
tmp<-tmp %>% filter(grepl("^>",X1))
map_name<-as_tibble(cbind(gene1=str_extract(tmp$X1,"\\|(.+?)\\|(.+?)_HUMAN",group = 2),
                          gene2=str_extract(tmp$X1,"GN=(.+?) ",group = 1)))

kinase_domain4<-inner_join(map_name,kinase_domain3,by=c("gene1"="kinase"))
kinase_domain4$kinase<-kinase_domain4$gene2
kinase_domain3<-kinase_domain4

kinase_list<-read_tsv("D:/FusionGene2/5.kinase/input/human_kinases_list.txt")
kinase_domain3<-kinase_domain3 %>% filter(kinase %in% kinase_list$`HGNC Name`)
length(unique(kinase_domain3$kinase)) ##113

kinase_domain_intact<-tibble()
for(i in unique(kinase_domain3$kinase)){
  # i<-unique(kinase_domain3$kinase)[1]
  sub_kinase_domain<-kinase_domain3 %>% filter(kinase==i)
  sub_kinase_domain<-sub_kinase_domain %>% filter(type=="Domain")
  
  sub_kinase_fusion_domain1<-kinase_fusion_domain3 %>% filter(G5==i)
  sub_kinase_fusion_domain2<-kinase_fusion_domain3 %>% filter(G3==i)
  sub_kinase_fusion_domain<-rbind(sub_kinase_fusion_domain1,sub_kinase_fusion_domain2) %>% distinct()
  sub_kinase_fusion_domain<-sub_kinase_fusion_domain  %>% filter(type=="Domain")
  
  common_domain<-inner_join(sub_kinase_domain,sub_kinase_fusion_domain,by=c("hmm acc"="hmm acc"))
  
  kinase_domain_intact<-rbind(kinase_domain_intact,common_domain)
}
kinase_domain_intact<-kinase_domain_intact %>% filter(`hmm name.x`=="Pkinase") %>% filter(`hmm length.x`==`hmm length.y`)


kinase_chimericevent<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/output_remove_cscRNA/kinasefusionevent_remove_cscRNA.txt")
kinase_chimericevent<-kinase_chimericevent[,1:5]
colnames(kinase_chimericevent)<-c("ID","G5","B5","G3","B3")

kinase_chimericevent<-kinase_chimericevent %>% mutate(Pkinase_domain=ifelse(ID %in% kinase_domain_intact$`seq id.y`,"Intact","Disrupted"))
write_tsv(kinase_chimericevent,"Figure6_new/kinase_chimericevent_pkinasedomain_if_intact.txt")

kinase_summary<-tibble()
for(i in unique(kinase_domain3$kinase)){
  # i<-unique(kinase_domain3$kinase)[1]
  sub1<-kinase_chimericevent %>% filter(G5==i)
  sub2<-kinase_chimericevent %>% filter(G3==i)
  sub<-rbind(sub1,sub2)
  sub<-sub %>% dplyr::select(ID,Pkinase_domain) %>% distinct()
  sub$kinase<-i
  kinase_summary<-rbind(kinase_summary,sub)
  
}
# kinase_domain_count<-kinase_summary %>% distinct() %>% group_by(kinase,Pkinase_domain) %>% summarise(count=n()) %>% dplyr::arrange(desc(count))
# # kinase_domain_count<-kinase_domain_count %>% filter(Pkinase_domain=="Intact")
# dir.create("Figure6_new")
# write_csv(kinase_domain_count,"Figure6_new/kinase_domain_fusionevents_number.csv")

kinase_domain_count<-kinase_summary %>% dplyr::select(ID,kinase) %>% distinct() %>% group_by(kinase) %>% summarise(chimericeventsNumber=dplyr::n())
# colnames(kinase_domain_count)<-c("kinase","count")
###===================================
fusionData<-read_tsv("D:/FusionGene2/share_data/ESCC_All_fusion_3957.changed_name.txt")
fusionData<-fusionData %>% filter(Fusion_Type %in% c("DNA level fusion","RNA level fusion")) 
fusion<-fusionData %>% dplyr::select(G5,Breakpoint5,G3,Breakpoint3,sample) %>% unique()

fusion <- fusion %>% separate(Breakpoint5,into=c("CHR1","POINT1","STRAND1"),sep=":")
fusion <- fusion %>% separate(Breakpoint3,into=c("CHR2","POINT2","STRAND2"),sep=":")
fusion<-fusion %>% select(G5,POINT1,G3,POINT2,sample) %>% unique()

fusion<-fusion %>% mutate(ID=paste0(G5,"-",POINT1,"_",G3,"-",POINT2))
fusion<-fusion %>% dplyr::select(ID,sample) %>% distinct()


# kinase_summary_intact<-kinase_summary %>% %>% distinct() %>% group_by(kinase) %>% summarise(count=n()) %>% dplyr::arrange(desc(count))
kinase_summary_intact<-kinase_summary %>% dplyr::select(ID,kinase) %>% distinct()
kinase_intact_patient_number<-inner_join(kinase_summary_intact,fusion)
kinase_intact_patient_number<-kinase_intact_patient_number %>% dplyr::select(kinase,sample) %>% distinct() 
kinase_intact_patient_number<-kinase_intact_patient_number %>% group_by(kinase) %>% summarise(PatientNumber=n()) %>% dplyr::arrange(desc(PatientNumber))
write_csv(kinase_intact_patient_number,"Figure6_new/kinase_domain_patient_number.csv")





# kinase_fusion_domain3<-kinase_fusion_domain2 %>% dplyr::select(`seq id`,`hmm acc`, `hmm name`,type) 
# kinase_domain_intact3<-kinase_domain_intact %>% dplyr::select(`seq id.x`,`seq id.y`) 
# colnames(kinase_domain_intact3)<-c("kinase_gene","chimeric events")
# colnames(kinase_fusion_domain3)<-c("chimeric events","hmm acc","hmm name","type")
# 
# dt<-left_join(kinase_fusion_domain3,kinase_domain_intact3)



###+=====================draw figure
# dt<-read_tsv("D:/FusionGene2/5.kinase/draw_kinase_domian/input/kinase_fusion_events_finial.txt")
# dt<-dt %>% dplyr::select(event,G5,G3,domain_if_intact,effect,`HGNC Name`,`Kinase Name`,Group, Family,SubFamily, UniprotID,Fusion_Type) %>% distinct()
# dt<-dt %>% dplyr::select(`HGNC Name`,`Kinase Name`,Group,Family,SubFamily,UniprotID)

intact_kinase_summary<-kinase_domain_count

# colnames(intact_kinase_summary)<-c("kinase","type","event_count")
# colnames(kinase_intact_patient_number)<-c("kinase","type","patient_count")
gene_order<-inner_join(intact_kinase_summary,kinase_intact_patient_number,by=c("kinase"="kinase"))
gene_order<-gene_order %>% dplyr::arrange(desc(PatientNumber),desc(chimericeventsNumber))

gene_order<-gene_order %>% filter(PatientNumber>1)
intact_kinase_summary$chimericeventsNumber<-(-intact_kinase_summary$chimericeventsNumber)
intact_kinase_summary<-intact_kinase_summary %>% filter(kinase %in% gene_order$kinase)
intact_kinase_summary$kinase<-factor(intact_kinase_summary$kinase,levels = rev(gene_order$kinase))



p1<-ggplot(intact_kinase_summary)+
  # geom_rect(data= Reces_table, inherit.aes = F,
  #           aes(xmin=Start_gene-0.5, xmax=End_gene+0.5, ymin=-Inf, ymax=+Inf),
  #           fill=Reces_table$Color, alpha=0.3)+
  geom_bar(aes(x=kinase,y=chimericeventsNumber,fill="#B88CD4"),
           stat = 'identity',position="stack",
           width = rep(0.8,nrow(intact_kinase_summary)))+ 
  labs(x='', y = 'The number of chimeric events') +
  coord_flip()+#翻转坐标轴
  scale_fill_manual(values = c("#B88CD4","#2fa61d"))+
  #修改坐标轴刻度和标签
  scale_y_continuous(breaks = seq(-65, 0, 5),
                     labels = as.character(abs(seq(-65, 0, 5))),
                     limits = c(-65, 0),expand = c(0,0)) +
  
  # scale_x_continuous(breaks = tmp1$`HGNC Name2`, 
  #                    labels = tmp1$`HGNC Name`, 
  #                    expand = c(0,0)) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_blank(),
        axis.text = element_text(color="black"),
        axis.ticks.y = element_blank()
  )
p1


kinase_intact_patient_number<-kinase_intact_patient_number %>% filter(kinase %in% gene_order$kinase)
kinase_intact_patient_number$kinase<-factor(kinase_intact_patient_number$kinase,levels = unique(rev(gene_order$kinase)))
p2<-ggplot(kinase_intact_patient_number)+
  # geom_rect(data= Reces_table, inherit.aes = F,
  #           aes(xmin=Start_gene-0.5, xmax=End_gene+0.5, ymin=-Inf, ymax=+Inf),
  #           fill=Reces_table$Color, alpha=0.3)+
  geom_bar(aes(x=kinase,y=PatientNumber,fill="#2fa61d"),
           stat = 'identity',position="stack",
           width = rep(0.8,nrow(kinase_intact_patient_number)))+ 
  labs(x='', y = 'The number of patients with intact kinase domain') +
  coord_flip()+#翻转坐标轴
  scale_fill_manual(values = c("#2fa61d"))+
  #修改坐标轴刻度和标签
  scale_y_continuous(breaks = seq(0, 30, 5),
                     labels = as.character(abs(seq(0, 30, 5))),
                     limits = c(0, 30),expand = c(0,0)) +
  
  # scale_x_continuous(breaks = tmp1$`HGNC Name2`, 
  #                    labels = tmp1$`HGNC Name`, 
  #                    expand = c(0,0)) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.background = element_blank(),
        # axis.text.y = element_blank(),
        axis.text = element_text(color="black"),
        axis.ticks.y = element_blank()
  )
p2


library(cowplot)
p<-plot_grid(p1,p2,ncol =2,rel_widths = c(1,1.3))
p
pdf("Figure6_new/kinases.pdf",width=6,height=7)
print(p)
dev.off()

tiff("Figure6_new/kinases.tiff",width=480,height=520)
print(p)
dev.off()

