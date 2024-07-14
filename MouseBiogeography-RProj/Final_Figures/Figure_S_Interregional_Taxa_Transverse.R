library(Microbiome.Biogeography)
library(dplyr)
library(ggplot2)
library(here)

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_Interregional_Taxa_Transverse.R")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

### Colon --- 

Type_cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
seecolor::print_color(c("#481567FF"))

hum_v_gavage_DAT_colon_result <- generate_interregional_taxa_barplot_TYPE("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                             "HUM MD Gavage Colon",
                                                             Type_cols)

hum_md_genera <-  (hum_v_gavage_DAT_colon_result$dataframe)$annotation
hum_v_gavage_DAT_colon <- hum_v_gavage_DAT_colon_result$plot

ucla_o_spf_DAT_colon_result <- generate_interregional_taxa_barplot_TYPE("Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                               "UCLA O. SPF Colon",
                                                               Type_cols)
ucla_o_genera <-  (ucla_o_spf_DAT_colon_result$dataframe)$annotation
ucla_o_spf_DAT_colon <- ucla_o_spf_DAT_colon_result$plot

intersect(hum_md_genera, ucla_o_genera)
cs_DAT_colon_result <- generate_interregional_taxa_barplot_TYPE("CS_SPF/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                             "CS SPF Colon",
                                                             Type_cols)
cs_genera <-  (cs_DAT_colon_result$dataframe)$annotation
cs_DAT_colon <- cs_DAT_colon_result$plot
intersect(ucla_o_genera,cs_genera)

cs_DAT_colon <- cs_DAT_colon_result$plot

spf_gavage_DAT_colon_result <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "SPF Gavage Colon",
                                                          Type_cols)

spf_genera <-  (spf_gavage_DAT_colon_result$dataframe)$annotation
spf_gavage_DAT_colon <- spf_gavage_DAT_colon_result$plot

hum_gavage_DAT_colon_result <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                          "HUM SD Gavage Colon",
                                                          Type_cols)
hum_sd_genera <- (hum_gavage_DAT_colon_result$dataframe)$annotation
hum_gavage_DAT_colon<- hum_gavage_DAT_colon_result$plot

# intersections 
length(hum_md_genera)
length(hum_sd_genera)
length(ucla_o_genera)
length(cs_genera)
length(spf_genera)
c1 <- intersect(hum_md_genera, ucla_o_genera)
c2 <- intersect(cs_genera,ucla_o_genera)
c3 <- intersect(spf_genera,ucla_o_genera)
c4 <- intersect(spf_genera,cs_genera)
c5 <- intersect(spf_genera,hum_md_genera)
c6 <- intersect(cs_genera, hum_md_genera)
c7 <- intersect(spf_genera, hum_sd_genera)
c8 <- intersect(hum_md_genera, hum_sd_genera)

intersect(c1,c2)

dev.new()
cowplot::plot_grid(ucla_o_spf_DAT_colon, hum_v_gavage_DAT_colon,ncol=2,
                   labels = c("A","E"))

dev.new()
cowplot::plot_grid(cs_DAT_colon,spf_gavage_DAT_colon, hum_gavage_DAT_colon,ncol=3,
                   labels=c("B","C","D"))

### Small Intestine --- 
Type_cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
seecolor::print_color(c("#481567FF"))

hum_v_gavage_DAT_si_result <- generate_interregional_taxa_barplot_TYPE("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                   "HUM MD Gavage SI",
                                                                   Type_cols)
hum_md_genera <-  (hum_v_gavage_DAT_si_result$dataframe)$annotation
hum_v_gavage_DAT_si <- hum_v_gavage_DAT_si_result$plot

ucla_o_spf_DAT_si_result <- generate_interregional_taxa_barplot_TYPE("Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv",
                                                                 "UCLA O. SPF SI",
                                                                 Type_cols)
ucla_o_genera <-  (ucla_o_spf_DAT_si_result$dataframe)$annotation
ucla_o_spf_DAT_si <- ucla_o_spf_DAT_si_result$plot

cs_DAT_si_result <- generate_interregional_taxa_barplot_TYPE("CS_SPF/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                         "CS SPF SI",
                                                         Type_cols)

cs_genera <-  (cs_DAT_si_result$dataframe)$annotation
cs_DAT_si <- cs_DAT_si_result$plot

spf_gavage_DAT_si_result <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "SPF Gavage SI",
                                                                 Type_cols)

spf_genera <-  (spf_gavage_DAT_si_result$dataframe)$annotation
spf_gavage_DAT_si <- spf_gavage_DAT_si_result$plot

hum_gavage_DAT_si <- generate_interregional_taxa_barplot_TYPE("Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID/significant_results.tsv",
                                                                 "HUM SD Gavage SI",
                                                                 Type_cols)

# intersections 
length(hum_md_genera)
length(ucla_o_genera)
length(cs_genera)
length(spf_genera)
c1 <- intersect(hum_md_genera, ucla_o_genera)
c2 <- intersect(cs_genera,ucla_o_genera)
c3 <- intersect(spf_genera,ucla_o_genera)
c4 <- intersect(spf_genera,cs_genera)
c5 <- intersect(spf_genera,hum_md_genera)
c6 <- intersect(cs_genera, hum_md_genera)

dev.new()
cowplot::plot_grid(ucla_o_spf_DAT_si, hum_v_gavage_DAT_si,ncol=2,
                   labels = c("A","E"))

dev.new()
cowplot::plot_grid(cs_DAT_si,spf_gavage_DAT_si, hum_gavage_DAT_si,ncol=3,
                   labels=c("B","C","D"))

