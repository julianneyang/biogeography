library(Microbiome.Biogeography)
library(dplyr)

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_Interregional_Taxa_Transverse.R")

### Colon --- 

Type_cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
seecolor::print_color(c("#481567FF"))

hum_v_gavage_DAT_colon <- generate_interregional_taxa_barplot_TYPE("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                             "HUM V. Gavage Colon",
                                                             Type_cols)

ucla_o_spf_DAT_colon <- generate_interregional_taxa_barplot_TYPE("Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                               "UCLA O. SPF Colon",
                                                               Type_cols)

cs_DAT_colon <- generate_interregional_taxa_barplot_TYPE("CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                             "CS SPF Colon",
                                                             Type_cols)

spf_gavage_DAT_colon <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "SPF Gavage Colon",
                                                          Type_cols)

hum_gavage_DAT_colon <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "HUM Gavage Colon",
                                                          Type_cols)

dev.new()
cowplot::plot_grid(ucla_o_spf_DAT_colon, hum_v_gavage_DAT_colon,ncol=2,
                   labels = c("A","E"))

dev.new()
cowplot::plot_grid(cs_DAT_colon,spf_gavage_DAT_colon, hum_gavage_DAT_colon,ncol=3,
                   labels=c("B","C","D"))

### Small Intestine --- 
Type_cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
seecolor::print_color(c("#481567FF"))

hum_v_gavage_DAT_si <- generate_interregional_taxa_barplot_TYPE("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                   "HUM V. Gavage SI",
                                                                   Type_cols)

ucla_o_spf_DAT_si <- generate_interregional_taxa_barplot_TYPE("Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                                 "UCLA O. SPF SI",
                                                                 Type_cols)

cs_DAT_si <- generate_interregional_taxa_barplot_TYPE("CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                         "CS SPF SI",
                                                         Type_cols)

spf_gavage_DAT_si <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "SPF Gavage SI",
                                                                 Type_cols)

hum_gavage_DAT_si <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "HUM Gavage SI",
                                                                 Type_cols)

dev.new()
cowplot::plot_grid(ucla_o_spf_DAT_si, hum_v_gavage_DAT_si,ncol=2,
                   labels = c("A","E"))

dev.new()
cowplot::plot_grid(cs_DAT_si,spf_gavage_DAT_si, hum_gavage_DAT_si,ncol=3,
                   labels=c("B","C","D"))

