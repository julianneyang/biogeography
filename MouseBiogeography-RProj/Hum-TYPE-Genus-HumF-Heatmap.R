###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)

###for TYPE:Mucosal vs Luminal Data
#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")
duodenum<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Duodenum-ComBat-SexType-1-MsID/significant_results.tsv", header=TRUE)
duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
a<-duodenum_significant$feature
jejunum<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
b<-jejunum_significant$feature
ileum<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
ileum_significant<-filter(ileum, metadata=="Type" & value=="Mucosal" &qval<0.05)
c<-ileum_significant$feature
cecum<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
d<-cecum_significant$feature  
pc<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
e<-pc_significant$feature  
DC<-read.table("Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef) #totally empty
