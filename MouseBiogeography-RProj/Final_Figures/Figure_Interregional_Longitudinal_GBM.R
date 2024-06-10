library(dplyr)
library(Microbiome.Biogeography)
library(here)

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Interregional_Longitudinal_GBM.R")

### Luminal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)

hum_v_GBM_Map_lum <- generate_interregional_GBM_barplot("Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "HUM MD Gavage Luminal",
                                                         cols)


ucla_o_GBM_Map_lum <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA O. SPF Luminal",
                                                         cols)
cs_spf_GBM_Map_lum <- generate_interregional_GBM_barplot("CS_SPF/OMIXER-RPM Results/CS_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "CS SPF Luminal",
                                                         cols)
spf_gavage_GBM_Map_lum <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "SPF Gavage Luminal",
                                                             cols)
hum_gavage_GBM_Map_lum <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "HUM SD Gavage Luminal",
                                                             cols)
dev.new()
cowplot::plot_grid(ucla_o_GBM_Map_lum, cs_spf_GBM_Map_lum, 
                   nrow=1, ncol=2,
                   labels=c("A", "B"),
                   label_size = 12)

dev.new()
cowplot::plot_grid(spf_gavage_GBM_Map_lum, hum_gavage_GBM_Map_lum, 
                   hum_v_GBM_Map_lum,
                   nrow=2, ncol=2,
                   labels=c("C","D", "E"),
                   label_size = 12)

### Mucosal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
hum_v_GBM_Map_muc <- generate_interregional_GBM_barplot("Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "HUM MD Gavage Mucosal",
                                                         cols)

ucla_o_GBM_Map_muc <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA O. SPF Mucosal",
                                                         cols)
ucla_v_GBM_Map_muc <- generate_interregional_GBM_barplot("UCLA_V_SPF_Analysis/OMIXER-RPM/WT_Val_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA V. SPF Mucosal",
                                                         cols)
cs_spf_GBM_Map_muc <- generate_interregional_GBM_barplot("CS_SPF/OMIXER-RPM Results/CS_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "CS SPF Mucosal",
                                                         cols)
spf_gavage_GBM_Map_muc <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "SPF Gavage Mucosal",
                                                             cols)
hum_gavage_GBM_Map_muc <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "HUM SD Gavage Mucosal", 
                                                             cols)

dev.new()
cowplot::plot_grid(ucla_o_GBM_Map_muc, ucla_v_GBM_Map_muc, 
                   cs_spf_GBM_Map_muc, spf_gavage_GBM_Map_muc,
                  labels=c("A","B","C", "D"), 
                  label_size = 16,
                  nrow=2, ncol=2)

dev.new()
cowplot::plot_grid(hum_gavage_GBM_Map_muc, hum_v_GBM_Map_muc,
                   nrow=1, ncol=2,
                   labels=c("E", "F"),
                   label_size = 16)



### Luminal Barplots: Aggregated by Metabolic_Map  --- 
cols <- viridis::viridis(2)
ucla_o_GBM_MMap_lum <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                          "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                          ystring="metabolic_map",
                                                          titlestring= "UCLA O. SPF Luminal",
                                                          cols)
cs_spf_GBM_MMap_lum <- generate_interregional_GBM_barplot("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                          ystring="metabolic_map",
                                                          titlestring= "CS SPF Luminal",
                                                          cols)
spf_gavage_GBM_MMap_lum <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                              ystring="metabolic_map",
                                                              titlestring= "SPF Gavage Luminal",
                                                              cols)
hum_gavage_GBM_MMap_lum <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/HUM/OMIXER-RPM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                              ystring="metabolic_map",
                                                              titlestring= "HUM Gavage Luminal",
                                                              cols)
cowplot::plot_grid(ucla_o_GBM_MMap_lum, cs_spf_GBM_MMap_lum, spf_gavage_GBM_MMap_lum, hum_gavage_GBM_MMap_lum, nrow=2, ncol=2)

### Mucosal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
ucla_o_GBM_Map_muc <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA O. SPF Mucosal",
                                                         cols)
ucla_v_GBM_Map_muc <- generate_interregional_GBM_barplot("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA V. SPF Mucosal",
                                                         cols)
cs_spf_GBM_Map_muc <- generate_interregional_GBM_barplot("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "CS SPF Mucosal",
                                                         cols)
spf_gavage_GBM_Map_muc <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "SPF Gavage Mucosal",
                                                             cols)
hum_gavage_GBM_Map_muc <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/HUM/OMIXER-RPM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "HUM Gavage Mucosal",
                                                             cols)
cowplot::plot_grid(ucla_o_GBM_Map_muc, ucla_v_GBM_Map_muc, cs_spf_GBM_Map_muc, spf_gavage_GBM_Map_muc, hum_gavage_GBM_Map_muc, nrow=2, ncol=3)
