###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 3.10.2023
### Figure Number: 6
### Figure Contents: GMM Site Heatmaps for all Cohorts 
###### whining ends here ---

library(ggplot2)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
library(plyr)
library(gridExtra)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_GMM_Site_Heatmap_Aggregated.R")

### HUM V Gavage ---
donors_filepath <- "Donors-Analysis/differential_GMM_site/"
lumtarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

muctarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
hum_v_map_lum <- generate_GMM_heat_map_by_site(paste0(donors_filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                                lumtarget,
                                                "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                Y=Map,
                                                "Map",
                                                "Luminal",
                                                cols,
                                                bk) +
  theme_cowplot(20) +
  ggtitle("HUM V. Gavage Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF")
bk =c( -1, -0.5, 0, 0.5)

hum_v_map_muc <- generate_GMM_heat_map_by_site(paste0(donors_filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                               muctarget,
                                               "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                               Y=Map,
                                               "Map",
                                               "Mucosal",
                                               cols,
                                               bk) +
  theme_cowplot(20) +
  ggtitle("HUM V. Gavage Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



### UCLA. O. SPF ---
lumtarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1)
ucla_o_map_lum <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk) +
  theme_cowplot(20) +
                  ggtitle("UCLA O. SPF Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-2,-1.5, -1, -0.5, 0, 0.5, 1)
ucla_o_map_muc <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     muctarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Mucosal",
                                     cols,
                                     bk) +
  theme_cowplot(20)+
                  ggtitle("UCLA O. SPF Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


### CS SPF ---

# Luminal - 36 concordant features 
lum_target <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 37 concordant features 
muctarget <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_muc_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Mucosal", 
                                             cols,
                                             bk) +
  
  theme_cowplot(20) +
                  ggtitle("CS SPF Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
cs_lum_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Luminal",
                                             cols,
                                             bk) +
  theme_cowplot(20) +
                  ggtitle("CS SPF Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


### HUM Gavage ---

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF")
bk =c(-2, -1.5, -1, -0.5, 0)

hum_lum_GMM_map <- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk) +
  theme_cowplot(20) +
  ggtitle("HUM Gavage Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
 



cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

hum_muc_GMM_map<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=Map,
                                    "Map",
                                    "Mucosal",
                                    cols,
                                    bk) +
  theme_cowplot(20) +
  ggtitle("HUM Gavage Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

hum_muc_GMM_map

### SPF Gavage ---
muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

spf_muc_GMM_map<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=Map,
                                    "Map",
                                    "SPF Gavage Muc",
                                    cols,
                                    bk)+
  theme_cowplot(20) +
  ggtitle("SPF Gavage Muc") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

fig6h <- plot_grid(NULL,spf_muc_GMM_map, ncol=1)

### UCLA V. SPF ---
muctarget <- find_concordant_features_across_sites("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_v_muc_GMM_map <- generate_GMM_heat_map_by_site("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             titlestring="UCLA V. SPF Muc",
                                             colorvector = cols,
                                             breakvector = bk) +
  theme_cowplot(20) +
  ggtitle("UCLA V. SPF Muc") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

thelegend <- generate_GMM_heat_map_by_site("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                           targetvector = muctarget, 
                                           path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                           Y=Map,
                                           ystring= "Map",
                                           titlestring="UCLA V. SPF Muc",
                                           colorvector = cols,
                                           breakvector = bk) +
  theme_cowplot(16) +
  ggtitle("UCLA V. SPF Muc") +
  theme(legend.position = "top") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- cowplot::get_legend(thelegend)
grid::grid.newpage()
grid::grid.draw(legend)

### Assemble Figure 6 ---

dev.new(width=10,height=10)
plot_grid(ucla_o_map_lum, cs_lum_GMM_map, hum_lum_GMM_map, nrow=1, labels= c("A", "B", "F"), label_size = 20,greedy = FALSE)
dev.new(width=10,height=10)
plot_grid(ucla_o_map_muc, cs_muc_GMM_map, hum_v_map_lum, nrow=1, labels= c("C", "D", "G"), label_size = 20)
dev.new(width=10,height=10)
plot_grid(ucla_v_muc_GMM_map, hum_muc_GMM_map, hum_v_map_muc, nrow=1, labels = c("H", "I", "J"), label_size = 20)
dev.new(width=10,height=10)
plot_grid(spf_muc_GMM_map, NULL, NULL, nrow=1, labels = c("K", "", "L"), label_size = 20)
