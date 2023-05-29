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
here::i_am("MouseBiogeography-RProj/HumGMM-SITE-Heatmap-SPF.R")

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF")
bk =c(-2, -1.5, -1, -0.5, 0)

lum <- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Luminal",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum


lum <- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum


cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

muc<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=Map,
                                    "Map",
                                    "Mucosal",
                                    cols,
                                    bk)
dev.new(width=15, height=10)
muc

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

muc <- read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")
muc<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=metabolic_map,
                                    "metabolic_map",
                                    "Mucosal",
                                    cols,
                                    bk)
dev.new(width=15, height=10)
muc

## Original color vector and breaks- see how to add this to default params 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

### Make heatmaps by hand 
