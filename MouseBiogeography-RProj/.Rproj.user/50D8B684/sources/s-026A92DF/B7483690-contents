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
here::i_am("MouseBiogeography-RProj/RegionalGBM-SITE-Heatmap.R")

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


# Luminal - 20 concordant features 
lum_target <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 30 concordant features 
muctarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv and make a heatmap---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

muc_GMM_map <- generate_GBM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "UCLA O. SPF Mucosal",
                                             cols,
                                             bk)

dev.new(width=15, height=10)
muc_GMM_map


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
lum_GMM_map <- generate_GBM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "UCLA O. SPF Luminal",
                                             cols,
                                             bk)

dev.new(width=15, height=10)
lum_GMM_map



