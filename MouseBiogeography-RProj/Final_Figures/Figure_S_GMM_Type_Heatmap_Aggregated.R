###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(here)
library(cowplot)
library(ggplot2)
library(dplyr)
library(plyr)

remove.packages("Microbiome.Biogeography")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_GMM_Type_Heatmap_Aggregated.R") #conflicts with plyr so use the :: notation
here::here()

### HUM V Gavage ---
target <- find_features_union_for_type_heatmap(
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GMM_type/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


#construct the heatmap using ggplot

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)

?generate_GMM_heat_map_by_Type()
hum_v_gmm_type <- generate_GMM_heat_map_by_Type(df, 
                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                         Map, 
                                         "Map", 
                                         "HUM V. Gavage", 
                                         colorvector=cols, 
                                         breakvector=bk)
dev.new(width=15, height=10)
hum_v_gmm_type

### UCLA O SPF ---
target <- find_features_union_for_type_heatmap(
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_GMM_type/GMM-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv")


#construct the heatmap using ggplot

cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1, -0.5, 0, 0.5, 1, 1.5, 2)


?generate_GMM_heat_map_by_Type()
ucla_o_gmm_type <- generate_GMM_heat_map_by_Type(df, 
                                                "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                                Map, 
                                                "Map", 
                                                "UCLA O. SPF", 
                                                colorvector=cols, 
                                                breakvector=bk)
dev.new(width=15, height=10)
ucla_o_gmm_type

### SPF Gavage---
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/SPF_GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


#construct the heatmap using ggplot

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF")
bk =c(-1.5, -1, -0.5, 0, 0.5)


?generate_GMM_heat_map_by_Type()
spf_gmm_type <- generate_GMM_heat_map_by_Type(df, 
                                                 "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                                 Map, 
                                                 "Map", 
                                                 "SPF Gavage", 
                                                 colorvector=cols, 
                                                 breakvector=bk)
dev.new(width=15, height=10)
spf_gmm_type

### HUM Gavage---
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_GMM_type/HUM_GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


#Note - these two features were omitted from cecum (due to overfitting) and therefore the map was not included

cols=c("#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-0.5, 0, 0.5, 1, 1.5, 2)


?generate_GMM_heat_map_by_Type()
hum_gmm_type <- generate_GMM_heat_map_by_Type(df, 
                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                              Map, 
                                              "Map", 
                                              "HUM Gavage", 
                                              colorvector=cols, 
                                              breakvector=bk)
dev.new(width=15, height=10)
hum_gmm_type

### CS SPF ---
target <- find_features_union_for_type_heatmap(
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_GMM_type/GMM-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


#construct the heatmap using ggplot

cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1, -0.5, 0, 0.5, 1, 1.5, 2)


?generate_GMM_heat_map_by_Type()
cs_gmm_type <- generate_GMM_heat_map_by_Type(df, 
                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                              Map, 
                                              "Map", 
                                              "CS SPF", 
                                              colorvector=cols, 
                                              breakvector=bk)
dev.new(width=15, height=10)
cs_gmm_type

### Plot Grid ---
top <- cowplot::plot_grid(ucla_o_gmm_type, hum_v_gmm_type,
          labels=c("A","E"))

bottom <- cowplot::plot_grid(cs_gmm_type, spf_gmm_type,hum_gmm_type,
                   labels=c("B","C","D"),nrow=1)

dev.new()
plot_grid(top,bottom,
          nrow=2)
