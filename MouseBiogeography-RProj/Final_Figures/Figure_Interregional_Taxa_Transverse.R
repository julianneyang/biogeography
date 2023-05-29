library(Microbiome.Biogeography)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/RegionalGBM-SITE-Heatmap.R")

### Luminal --- 

Type_cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
seecolor::print_color(c("#481567FF"))
ucla_o_DAT_colon <- generate_interregional_taxa_dotplot_TYPE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                      "ucla_original",
                                                      "UCLA O. SPF Colon",
                                                      Type_cols)

cs_DAT_colon <- generate_interregional_taxa_dotplot_TYPE("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                             "cs",
                                                             "CS. SPF Colon",
                                                             Type_cols)

spf_gavage_DAT_colon <- generate_interregional_taxa_dotplot_TYPE("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Type_L6/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "spf_gavage",
                                                          "SPF Gavage Colon",
                                                          Type_cols)

hum_gavage_DAT_colon <- generate_interregional_taxa_dotplot_TYPE("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "hum_gavage",
                                                          "HUM Gavage Colon",
                                                          Type_cols)

cowplot::plot_grid(ucla_o_DAT_colon, cs_DAT_colon,ncol=2)
cowplot::plot_grid(spf_gavage_DAT_colon, hum_gavage_DAT_colon,ncol=2)

### Mucosal --- 
ucla_o_DAT_si <- generate_interregional_taxa_dotplot_TYPE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/L6-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                             "ucla_original",
                                                             "UCLA O. SPF SI",
                                                             Type_cols)
# Dropped feature: 861fb4882b4c65b68e5998bb89ad59e8
cs_DAT_si <- generate_interregional_taxa_dotplot_TYPE("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                         "cs",
                                                         "CS. SPF SI",
                                                         Type_cols)

spf_gavage_DAT_si <- generate_interregional_taxa_dotplot_TYPE("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Type_L6/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "spf_gavage",
                                                                 "SPF Gavage SI",
                                                                 Type_cols)

hum_gavage_DAT_si <- generate_interregional_taxa_dotplot_TYPE("Humanized-Biogeography-Analysis/Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "hum_gavage",
                                                                 "HUM Gavage SI",
                                                                 Type_cols)

cowplot::plot_grid(ucla_o_DAT_si, cs_DAT_si,ncol=2)
cowplot::plot_grid(spf_gavage_DAT_si, hum_gavage_DAT_si,ncol=2)
