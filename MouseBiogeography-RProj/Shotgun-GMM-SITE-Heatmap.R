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
here::i_am("MouseBiogeography-RProj/HumGMM-SITE-Heatmap.R")

## Original color vector and breaks- see how to add this to default params 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Shotgun/UCLA_O_SPF/GMM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)

lum <- generate_GMM_heat_map_by_site("Shotgun/UCLA_O_SPF/GMM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Shotgun UCLA",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Shotgun/CS_SPF/GMM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

lum <- generate_GMM_heat_map_by_site("Shotgun/CS_SPF/GMM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Shotgun CS",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum


#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Shotgun/SPF_Gavage/GMM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

lum <- generate_GMM_heat_map_by_site("Shotgun/SPF_Gavage/GMM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Shotgun SPF Gavage",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Shotgun/HUM_Gavage/GMM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

lum <- generate_GMM_heat_map_by_site("Shotgun/HUM_Gavage/GMM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Shotgun HUM Gavage",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum
