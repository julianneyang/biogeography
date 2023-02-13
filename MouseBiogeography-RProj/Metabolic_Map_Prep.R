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

###for SITE:DC vs all data
#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER")
input<-read.table("GMM-Maaslin2-SITE/GMM-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
input<-read.table("GMM-Maaslin2-SITE/GMM_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
data<-filter(input, metadata=="Site_General" & qval<0.05)
input<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",header = TRUE)
input<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",header = TRUE)

data<-filter(input, metadata=="Type" & qval<0.05)


site_heatmap<-data

annotation <- read.csv("Revised_Module_Key.csv", header=TRUE)

data<- (merge(site_heatmap, annotation, by = 'feature'))

write.csv(data,"MMAP_Luminal-ColonvsSI-GMM.csv")
write.csv(data,"MMAP_Mucosal-ColonvsSI-GMM.csv")
write.csv(data,"MMAP_LumRef-Colon-GMM.csv")
write.csv(data,"MMAP_LumRef-SI-GMM.csv")
