###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(data.table)
library(funrar)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)
library(plyr)
library(cowplot)
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_L2_HEATMAP_Site_all.R")

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

#Full heatmap colors- modify this according to max and min coef sizes
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

### HUM V Gavage ---
lumtarget <- find_concordant_features_across_sites("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

hum_v_lum_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                                                        lumtarget,
                                                        "HUM V. Gavage Luminal",
                                                        cols,
                                                        bk)

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1)

hum_v_muc_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                                                        muctarget,
                                                        "HUM V. Gavage Mucosal",
                                                        cols,
                                                        bk)



### UCLA O SPF ---
lumtarget <- find_concordant_features_across_sites("Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-1, -0.5, 0, 0.5, 1)

ucla_o_lum_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "UCLA O. SPF Luminal",
                                     cols,
                                     bk)

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1)

ucla_o_muc_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                      muctarget,
                                     "UCLA O. SPF Mucosal",
                                     cols,
                                     bk)

### UCLA V SPF ---
muctarget <- find_concordant_features_across_sites("Maaslin2_L2/UCLA_V_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_v_muc_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/UCLA_V_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                      muctarget,
                                      "UCLA V. SPF Mucosal",
                                      cols,
                                      bk)

### CS SPF ---
#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
lumtarget <- find_concordant_features_across_sites("Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cs_lum_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                        lumtarget,
                                                        "CS SPF Luminal",
                                                        cols,
                                                        bk)

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
muctarget <- find_concordant_features_across_sites("Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cs_muc_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                        muctarget,
                                                        "CS SPF Mucosal",
                                                        cols,
                                                        bk)

### SPF Gavage ---
lumtarget <- find_concordant_features_across_sites("Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# No significantly diff abundant phyla

### HUM Gavage ---
lumtarget <- find_concordant_features_across_sites("Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")
# No significantly diff abundant phyla in lum

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
muctarget <- find_concordant_features_across_sites("Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

hum_muc_L2_heatmap <- generate_taxa_heat_map_by_site("Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                    muctarget,
                                                    "HUM Gavage Mucosal",
                                                    cols,
                                                    bk)

dev.new()
plot_grid(ucla_o_muc_L2_heatmap, ucla_v_muc_L2_heatmap, cs_muc_L2_heatmap, hum_muc_L2_heatmap,
          hum_v_muc_L2_heatmap, ucla_o_lum_L2_heatmap,  cs_lum_L2_heatmap, hum_v_lum_L2_heatmap,
          nrow=3, ncol=3, labels=c("A", "B", "C", "D", "E", "F","G","H","I"))


