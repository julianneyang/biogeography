###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(cowplot)
library(here)

here::i_am("MouseBiogeography-RProj/Final_Figures/L2_HEATMAP_Type_all.R")

remove.packages("Microbiome.Biogeography")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

### Full heatmap colors- modify this according to max and min coef sizes ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

### UCLA O SPF ---
target <- find_features_union_for_type_heatmap(
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv")
  
df <- query_type_features_union(
  target,
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv")


#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_o_type_L2_heatmap <- generate_taxa_heat_map_by_type(df,
                               "UCLA O. SPF",
                               cols,
                               bk)

### CS SPF ---
target <- find_features_union_for_type_heatmap(
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")
  
df <- query_type_features_union(
  target,
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS_SPF/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")

#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_type_L2_heatmap <- generate_taxa_heat_map_by_type(df,
                                                         "CS SPF",
                                                         cols,
                                                         bk)

### HUM Gavage --- NO features
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")


### SPF Gavage --- no features
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

### HUM V Gavage ---
target <- find_features_union_for_type_heatmap(
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")
  
df <- query_type_features_union(
  target,
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")

#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
hum_v_gavage_type_L2_heatmap <- generate_taxa_heat_map_by_type(df,
                                                             "HUM MD Gavage",
                                                             cols,
                                                             bk)


heatmap_hum <- plot_grid(NULL,hum_type_L2_heatmap,NULL, nrow=3)                                                      

dev.new()
plot_grid(ucla_o_type_L2_heatmap, cs_type_L2_heatmap, 
          hum_v_gavage_type_L2_heatmap, NULL,
          nrow=2, ncol=2,
          labels=c("A","B","C", "D"))
