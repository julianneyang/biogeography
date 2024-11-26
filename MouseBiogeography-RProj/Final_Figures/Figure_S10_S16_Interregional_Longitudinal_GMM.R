library(dplyr)
library(here)

library("Microbiome.Biogeography")
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S10_S16_Interregional_Longitudinal_GMM.R")

### Luminal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
hum_v_GMM_Map_lum <- generate_interregional_GMM_barplot("Donors-Analysis/differential_GMM_site/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "HUM MD Gavage",
                                                         cols)

ucla_o_GMM_Map_lum <- generate_interregional_GMM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                   "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                   ystring="Map",
                                   titlestring= "UCLA O. SPF",
                                   cols)
cs_spf_GMM_Map_lum <- generate_interregional_GMM_barplot("CS_SPF/OMIXER-RPM Results/CS_GMM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                     ystring="Map",
                                                     titlestring= "CS SPF",
                                                     cols)
spf_gavage_GMM_Map_lum <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                     ystring="Map",
                                                     titlestring= "SPF Gavage",
                                                     cols)
hum_gavage_GMM_Map_lum <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "HUM SD Gavage",
                                                        cols)
dev.new()
cowplot::plot_grid(ucla_o_GMM_Map_lum, cs_spf_GMM_Map_lum, 
                   spf_gavage_GMM_Map_lum, hum_gavage_GMM_Map_lum,hum_v_GMM_Map_lum, nrow=2,
                   labels=c("A","B","C","D","E"), label_size = 16)

### Mucosal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
hum_v_GMM_Map_muc <- generate_interregional_GMM_barplot("Donors-Analysis/differential_GMM_site/GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "HUM MD Gavage",
                                                         cols)


ucla_o_GMM_Map_muc <- generate_interregional_GMM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                     ystring="Map",
                                                     titlestring= "UCLA O. SPF",
                                                     cols)
ucla_v_GMM_Map_muc <- generate_interregional_GMM_barplot("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "UCLA V. SPF",
                                                         cols)
cs_spf_GMM_Map_muc <- generate_interregional_GMM_barplot("CS_SPF/OMIXER-RPM Results/CS_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                     ystring="Map",
                                                     titlestring= "CS SPF",
                                                     cols)
spf_gavage_GMM_Map_muc <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "SPF Gavage",
                                                         cols)
hum_gavage_GMM_Map_muc <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="Map",
                                                         titlestring= "HUM SD Gavage",
                                                         cols)
dev.new()
cowplot::plot_grid(ucla_o_GMM_Map_muc, ucla_v_GMM_Map_muc, cs_spf_GMM_Map_muc, 
                   spf_gavage_GMM_Map_muc, hum_gavage_GMM_Map_muc, hum_v_GMM_Map_muc, 
                   nrow=2, ncol=3,
                   labels=c("A","B","C","D","E","F"),
                   label_size = 12)



### Luminal Barplots: Aggregated by Metabolic_Map  --- 
cols <- viridis::viridis(2)
ucla_o_GMM_MMap_lum <- generate_interregional_GMM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA O. SPF Luminal",
                                                         cols)
cs_spf_GMM_MMap_lum <- generate_interregional_GMM_barplot("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "CS SPF Luminal",
                                                         cols)
spf_gavage_GMM_MMap_lum <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "SPF Gavage Luminal",
                                                             cols)
hum_gavage_GMM_MMap_lum <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/HUM/OMIXER-RPM/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "HUM Gavage Luminal",
                                                             cols)
cowplot::plot_grid(ucla_o_GMM_MMap_lum, cs_spf_GMM_MMap_lum, spf_gavage_GMM_MMap_lum, hum_gavage_GMM_MMap_lum, nrow=2, ncol=2)

### Mucosal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
ucla_o_GMM_Map_muc <- generate_interregional_GMM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA O. SPF Mucosal",
                                                         cols)
ucla_v_GMM_Map_muc <- generate_interregional_GMM_barplot("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA V. SPF Mucosal",
                                                         cols)
cs_spf_GMM_Map_muc <- generate_interregional_GMM_barplot("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "CS SPF Mucosal",
                                                         cols)
spf_gavage_GMM_Map_muc <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "SPF Gavage Mucosal",
                                                             cols)
hum_gavage_GMM_Map_muc <- generate_interregional_GMM_barplot("Humanized-Biogeography-Analysis/Source RPCA/HUM/OMIXER-RPM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "HUM Gavage Mucosal",
                                                             cols)
cowplot::plot_grid(ucla_o_GMM_Map_muc, ucla_v_GMM_Map_muc, cs_spf_GMM_Map_muc, spf_gavage_GMM_Map_muc, hum_gavage_GMM_Map_muc, nrow=2, ncol=3)
