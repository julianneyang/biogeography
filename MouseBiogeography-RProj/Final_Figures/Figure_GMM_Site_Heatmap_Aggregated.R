###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 3.10.2023
### Figure Number: 6
### Figure Contents: GMM Site Heatmaps for all Cohorts 
###### whining ends here ---

library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
library(Microbiome.Biogeography)
library(plyr)
library(gridExtra)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")

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
                  ggtitle("UCLA O. SPF Luminal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")


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
                  ggtitle("UCLA O. SPF Mucosal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")


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
                  ggtitle("CS SPF Mucosal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none") 
cs_muc_GMM_map

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
                  ggtitle("CS SPF Luminal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")

cs_lum_GMM_map

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
  ggtitle("HUM Gavage Luminal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") 
 



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
  ggtitle("HUM Gavage Mucosal") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") 

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
                                    "SPF Gavage Mucosal",
                                    cols,
                                    bk)+
  theme_cowplot(20) +
  ggtitle("SPF Gavage Mucosal") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") 

fig6h <- plot_grid(NULL,spf_muc_GMM_map, ncol=1)

### UCLA V. SPF ---
muctarget <- find_concordant_features_across_sites("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_v_muc_GMM_map <- generate_GMM_heat_map_by_site("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             titlestring="UCLA V. SPF Mucosal",
                                             colorvector = cols,
                                             breakvector = bk) +
  theme_cowplot(20) +
  ggtitle("UCLA V. SPF Mucosal") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") 

thelegend <- generate_GMM_heat_map_by_site("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                           targetvector = muctarget, 
                                           path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                           Y=Map,
                                           ystring= "Map",
                                           titlestring="UCLA V. SPF Mucosal",
                                           colorvector = cols,
                                           breakvector = bk) +
  theme_cowplot(16) +
  ggtitle("UCLA V. SPF Mucosal") +
  theme(legend.position = "top") 

legend <- cowplot::get_legend(thelegend)
grid::grid.newpage()
grid::grid.draw(legend)

### Assemble Figure 6 ---

plot_grid(ucla_o_map_lum, cs_lum_GMM_map, hum_lum_GMM_map, nrow=1, labels= c("A", "B", "F"), label_size = 20)
plot_grid(ucla_o_map_muc, cs_muc_GMM_map, hum_muc_GMM_map, nrow=1, labels= c("C", "D", "G"), label_size = 20)
plot_grid(ucla_v_muc_GMM_map, NULL, fig6h, nrow=1, labels = c("E", "", "H"), label_size = 20)
