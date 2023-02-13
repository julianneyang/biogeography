###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
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

here::i_am("MouseBiogeography-RProj/CS_Facility_Metabolic_Map_Prep.R")
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", header=TRUE)

### for SI vs Colon ---
input<-read.table("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
data<-filter(input, metadata=="Site_General" & qval<0.05)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data,"CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/MMAP_Luminal-ColonvsSI-GMM.csv")

input<-read.table("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
data<-filter(input, metadata=="Site_General" & qval<0.05)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data,"CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/MMAP_Mucosal-ColonvsSI-GMM.csv")

### for Muc vs Lum ---
input<-read.table("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",header = TRUE)
data<-filter(input, metadata=="Type" & qval<0.05)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data,"CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/MMAP_LumRef-Colon-GMM.csv")

input<-read.table("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",header = TRUE)
data<-filter(input, metadata=="Type" & qval<0.05)
data<- (merge(data, annotation, by = 'feature'))
write.csv(data,"CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/MMAP_LumRef-SI-GMM.csv")

