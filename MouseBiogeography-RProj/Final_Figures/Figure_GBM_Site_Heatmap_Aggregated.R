###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 12.13.2022
### Figure Number: Supplemental GBM 
### Figure Contents: GBM Luminal and Mucosal Intraregional Heatmaps
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)

setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")


setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")

### HUM V Gavage ---
donors_filepath <- "Donors-Analysis/differential_GBM_site/"
lumtarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

muctarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c( -1, -0.5, 0, 0.5, 1)


hum_v_muc_GBM_map <- generate_GBM_heat_map_by_site(paste0(donors_filepath,"GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                                    targetvector = muctarget, 
                                                    path_to_Module_Key = "GBM_Module_Key.csv",
                                                    "HUM V Gavage Mucosal",
                                                    cols,
                                                    bk)


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
hum_v_lum_GBM_map <- generate_GBM_heat_map_by_site(paste0(donors_filepath,"GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                                   targetvector = lumtarget, 
                                                   path_to_Module_Key = "GBM_Module_Key.csv",
                                                   "HUM V Gavage Luminal",
                                                   cols,
                                                   bk)



### UCLA O SPF ---
# Luminal - 20 concordant features 
lum_target <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 30 concordant features 
muctarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv and make a heatmap---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)

ucla_o_muc_GBM_map <- generate_GBM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "UCLA O. SPF Mucosal",
                                             cols,
                                             bk)


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_o_lum_GBM_map <- generate_GBM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "UCLA O. SPF Luminal",
                                             cols,
                                             bk)

### CS SPF ---
# Luminal - 13 concordant features 
lum_target <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 17 concordant features 
muctarget <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv and make a heatmap---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

cs_muc_GBM_map <- generate_GBM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "CS SPF Mucosal",
                                             cols,
                                             bk)


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_lum_GBM_map <- generate_GBM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "CS_SPF_Luminal",
                                             cols,
                                             bk)

### UCLA V SPF ---
# Mucosal - 24 concordant features 
muctarget <- find_concordant_features_across_sites("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WT_Val_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv and make a heatmap---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

ucla_v_muc_GBM_map <- generate_GBM_heat_map_by_site("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WT_Val_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "UCLA V. SPF Mucosal",
                                             cols,
                                             bk)

### HUM Gavage ---
# Luminal - 1 concordant features 
lum_target <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 17 concordant features 
muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv and make a heatmap---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

hum_muc_GBM_map <- generate_GBM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "HUM Gavage Mucosal",
                                             cols,
                                             bk)


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF")
bk =c(-2, -1.5, -1, -0.5, 0)
hum_lum_GBM_map <- generate_GBM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "GBM_Module_Key.csv",
                                             "HUM Gavage Luminal",
                                             cols,
                                             bk)

### SPF Gavage ---

# Luminal - 0 concordant features 
lum_target <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 0 concordant features 
muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Generate multi-panel figure ---

top <- plot_grid(ucla_o_lum_GBM_map, cs_lum_GBM_map, 
          nrow=1, ncol=2, labels=c("A","B"))

middle <- plot_grid(ucla_o_muc_GBM_map, cs_muc_GBM_map, 
          nrow=1, ncol=2, labels=c("C","D"))

bottom_right <- plot_grid(hum_lum_GBM_map,
                          NULL,NULL,NULL, 
                          ncol=1,labels=c("G","H","",""))
bottom <- plot_grid(ucla_v_muc_GBM_map,
          hum_muc_GBM_map, bottom_right,
          nrow=1, ncol=3,
          labels=c("E","F"))
plot_grid(middle, bottom, nrow=2)
