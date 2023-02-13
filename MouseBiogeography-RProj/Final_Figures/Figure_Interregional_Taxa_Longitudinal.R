library(Microbiome.Biogeography)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/RegionalGBM-SITE-Heatmap.R")

### Luminal --- 
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
cols <- cols_general
ucla_o_DAT_lum <- generate_interregional_taxa_dotplot("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                    "ucla_original_lum",
                                    "UCLA O. SPF Luminal",
                                    cols)

cs_DAT_lum <- generate_interregional_taxa_dotplot("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "cs",
                                                      "CS SPF Luminal",
                                                      cols)

spf_gavage_DAT_lum <- generate_interregional_taxa_dotplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "spf_gavage",
                                                  "SPF Gavage Luminal",
                                                  cols)

hum_gavage_DAT_lum <- generate_interregional_taxa_dotplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                    "hum_gavage",
                                    "HUM Gavage Luminal",
                                    cols)

cowplot::plot_grid(ucla_o_DAT_lum, cs_DAT_lum,ncol=2)
cowplot::plot_grid(spf_gavage_DAT_lum, hum_gavage_DAT_lum,ncol=2)

 ### Mucosal --- 
ucla_o_DAT_muc <- generate_interregional_taxa_dotplot("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                      "ucla_original_muc",
                                                      "UCLA O. SPF Luminal",
                                                      cols)

ucla_v_DAT_muc <- generate_interregional_taxa_dotplot("ImmDef-Mouse-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "ucla_validation",
                                                      "UCLA V. SPF Mucosal",
                                                      cols)

cs_DAT_muc<- generate_interregional_taxa_dotplot("CS-Facility-Analysis/Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "cs",
                                                  "CS SPF Mucosal",
                                                  cols)


spf_gavage_DAT_muc <- generate_interregional_taxa_dotplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "spf_gavage",
                                                          "SPF Gavage Mucosal",
                                                          cols)

hum_gavage_DAT_muc <- generate_interregional_taxa_dotplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "hum_gavage",
                                                          "HUM Gavage Mucosal",
                                                          cols)

cowplot::plot_grid(ucla_o_DAT_muc, ucla_v_DAT_muc, cs_DAT_muc,ncol=3)
cowplot::plot_grid(spf_gavage_DAT_muc, hum_gavage_DAT_muc,ncol=2)
