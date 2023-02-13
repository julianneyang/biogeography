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
library(plyr)
library(Microbiome.Biogeography)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/HumGBM-SITE-Heatmap-SPF.R")

###for SITE:DC vs all data
find_concordant_features_across_sites <- function(filepath_to_significant_results_tsv) {
  luminal<-read.table(filepath_to_significant_results_tsv, header=TRUE)
  duodenum_significant<-filter(luminal, metadata=="Site" & value=="Duodenum" &qval<0.05)
  a<-duodenum_significant$feature
  jejunum_significant<-filter(luminal, metadata=="Site" & value=="Jejunum" &qval<0.05)
  b<-jejunum_significant$feature
  ileum_significant<-filter(luminal, metadata=="Site" & value=="Ileum" &qval<0.05)
  c<-ileum_significant$feature
  cecum_significant<-filter(luminal, metadata=="Site" & value=="Cecum" &qval<0.05)
  d<-cecum_significant$feature  
  pc_significant<-filter(luminal, metadata=="Site" & value=="Proximal_Colon" &qval<0.05)
  e<-pc_significant$feature  
  DC_significant<-filter(luminal, metadata=="Site" & value=="Distal_Colon" &qval<0.05)
  f<-DC_significant$feature  
  joinab<- union(a,b)
  joincd<- union(c,d)
  joinef<- union(e,f)
  joinabcd <- union(joinab,joincd)
  target<-union(joinabcd,joinef)
  return(target)
}


# Luminal - 0 concordant features 
lum_target <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 0 concordant features 
muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

