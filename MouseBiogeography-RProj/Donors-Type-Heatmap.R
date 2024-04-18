library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)
here::i_am("MouseBiogeography-RProj/Donors-Type-Heatmap.R")


target <- find_features_union_for_type_heatmap(duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Proximal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Distal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv")

#Query the target vector against all_results.tsv for each site
duodenum<-read.delim("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID-DonorID/all_results.tsv", header=TRUE)
duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
duodenum_all$Site<- "Duodenum"
jejunum<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
jejunum_all$Site<- "Jejunum"
ileum<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
ileum_all<-ileum_all[match(target,ileum_all$feature),]
ileum_all$Site<- "Ileum"
cecum<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
cecum_all<-cecum_all[match(target,cecum_all$feature),]
cecum_all$Site<- "Cecum"
pc<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
pc_all<-pc_all[match(target,pc_all$feature),]
pc_all$Site<- "Proximal_Colon"
DC<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
DC_all<-DC_all[match(target,DC_all$feature),]
DC_all$Site<- "Distal_Colon"

duojej<-rbind(duodenum_all,jejunum_all)
ilecec<-rbind(ileum_all, cecum_all)
pcdc<-rbind(pc_all,DC_all)
duojejilecec<-rbind(duojej,ilecec)
duojejilececpcdc<-rbind(duojejilecec,pcdc)