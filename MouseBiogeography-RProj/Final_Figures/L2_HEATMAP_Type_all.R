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
library(cowplot)
here::i_am("MouseBiogeography-RProj/Final_Figures/L2_HEATMAP_Type_all.R")

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

### Full heatmap colors- modify this according to max and min coef sizes ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

### UCLA O SPF ---
target <- find_features_union_for_type_heatmap(
  "Maaslin2_L2/UCLA_O_SPF/L2-Duodenum-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Jejunum-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Ileum-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Cecum-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-PC-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-DC-CLR-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Maaslin2_L2/UCLA_O_SPF/L2-Duodenum-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Jejunum-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Ileum-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-Cecum-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-PC-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/UCLA_O_SPF/L2-DC-CLR-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv"
) 

#from here make sure all NA rows are filled with feature name corresponding to NA
#write.csv(df, "Maaslin2_L2/UCLA_O_SPF/L2_type_heatmap.csv")

df <- read.csv("Maaslin2_L2/UCLA_O_SPF/L2_type_heatmap.csv")
heatmap<-df
discard<- heatmap[is.na(heatmap$metadata), ]
offtarget<- discard$feature
offtarget<-unique(offtarget)
heatmap_final<-subset(heatmap,  !heatmap[,3] %in% offtarget )

write.csv(offtarget, "Maaslin2_L2/UCLA_O_SPF/omitted_phyla.csv")

#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
ucla_o_type_L2_heatmap <- generate_taxa_heat_map_by_type(heatmap_final,
                               "UCLA O. SPF",
                               cols,
                               bk)

### CS SPF ---
target <- find_features_union_for_type_heatmap(
  "Maaslin2_L2/CS_SPF/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Maaslin2_L2/CS_SPF/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/CS_SPF/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv"
) 

#from here make sure all NA rows are filled with feature name corresponding to NA
  #and that coef = 0 and qval = 100 for those features 
#write.csv(df, "Maaslin2_L2/CS_SPF/L2_type_heatmap.csv")

df <- read.csv("Maaslin2_L2/CS_SPF/L2_type_heatmap.csv")
heatmap<-df
discard<- heatmap[is.na(heatmap$metadata), ]
offtarget<- discard$feature
offtarget<-unique(offtarget)
heatmap_final<-subset(heatmap,  !heatmap[,3] %in% offtarget )

write.csv(offtarget, "Maaslin2_L2/CS_SPF/omitted_phyla.csv")

#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1)
cs_type_L2_heatmap <- generate_taxa_heat_map_by_type(heatmap_final,
                                                         "CS SPF",
                                                         cols,
                                                         bk)
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_type_L2_heatmap_all <- generate_taxa_heat_map_by_type(heatmap,
                                                     "CS SPF all",
                                                     cols,
                                                     bk)

### SPF Gavage --- no features stand out 
target <- find_features_union_for_type_heatmap(
  "Maaslin2_L2/SPF_Gavage/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/SPF_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/SPF_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/SPF_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/SPF_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/SPF_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

### HUM Gavage ---
target <- find_features_union_for_type_heatmap(
  "Maaslin2_L2/HUM_Gavage/L2-Duodenum-CLR-ComBat-SexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Maaslin2_L2/HUM_Gavage/L2-Duodenum-CLR-ComBat-SexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Maaslin2_L2/HUM_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID/all_results.tsv"
) 

heatmap<-df
discard<- heatmap[is.na(heatmap$metadata), ]
offtarget<- discard$feature
offtarget<-unique(offtarget)
heatmap_final<-subset(heatmap,  !heatmap[,3] %in% offtarget )

#draw heatmap
cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
hum_gavage_type_L2_heatmap <- generate_taxa_heat_map_by_type(heatmap,
                                                     "HUM Gavage",
                                                     cols,
                                                     bk)

plot_grid(ucla_o_type_L2_heatmap, cs_type_L2_heatmap, hum_gavage_type_L2_heatmap, nrow=2, ncol=2,
          labels=c("A","B","C"))
