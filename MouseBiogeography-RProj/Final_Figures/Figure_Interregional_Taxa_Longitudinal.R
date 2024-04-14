library(Microbiome.Biogeography)
library(dplyr)
library(here)
library(ggplot2)

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Interregional_Taxa_Longitudinal.R")

### Luminal --- 
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
cols <- cols_general
ucla_o_DAT_lum <- generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                    "UCLA O. SPF Luminal",
                                    cols)
hum_v_DAT_lum <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                      "HUM V. Gavage Luminal",
                                                      cols)

cs_DAT_lum <- generate_interregional_taxa_barplot_SITE("CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "CS SPF Luminal",
                                                      cols)

spf_gavage_DAT_lum <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "SPF Gavage Luminal",
                                                  cols)

hum_gavage_DAT_lum <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                    "HUM Gavage Luminal",
                                    cols)

dev.new()
cowplot::plot_grid(ucla_o_DAT_lum, cs_DAT_lum,spf_gavage_DAT_lum, ncol=3,
                   labels=c("A","B","C"))
dev.new()
cowplot::plot_grid(hum_gavage_DAT_lum,hum_v_DAT_lum,ncol=2,
                   labels=c("D","E"))

 ### Mucosal --- 
ucla_o_DAT_muc <- generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                      "UCLA O. SPF Mucosal",
                                                      cols)

ucla_v_DAT_muc <- generate_interregional_taxa_barplot_SITE("UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "UCLA V. SPF Mucosal",
                                                      cols)

cs_DAT_muc<- generate_interregional_taxa_barplot_SITE("CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "CS SPF Mucosal",
                                                  cols)


hum_gavage_DAT_muc <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "HUM Gavage Mucosal",
                                                          cols)

hum_v_gavage_DAT_muc <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                               "HUM V. Gavage Mucosal",
                                                               cols)


spf_gavage_DAT_muc <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "SPF Gavage Mucosal",
                                                          cols)

dev.new()
cowplot::plot_grid(ucla_o_DAT_muc, ucla_v_DAT_muc, cs_DAT_muc,ncol=3,
                   labels=c("A","B","C"),
                   label_size = 16)

dev.new()
cowplot::plot_grid(spf_gavage_DAT_muc, hum_gavage_DAT_muc, hum_v_gavage_DAT_muc, 
                   ncol=3, labels=c("D","E","F"),
                   label_size = 16)
