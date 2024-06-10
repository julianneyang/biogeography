library(Microbiome.Biogeography)
library(dplyr)
library(here)
library(tidyr)
library(ggplot2)

remove.packages("Microbiome.Biogeography")
setwd("../microbiome.biogeography/")
devtools::document()
setwd("../biogeography/")
#devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Interregional_Taxa_Longitudinal.R")


### Function ---
data_all <- NULL
process_results_for_upset_plot_interregional <- function(file_paths, cohort_prefixes) {

  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.delim(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site_General" & qval<0.05)
    
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    data <- data %>% mutate(Cohort = cohort)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  
  return(data_all)
}

### Identify overlapping features ---
lum_file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/all_results.tsv",  
  "CS_SPF/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv")

lum_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "SPF_Gavage",
                         "HUM_SD_Gavage")

all_taxa <- process_results_for_upset_plot_interregional(file_paths = lum_file_paths,
                                           cohort_prefixes = lum_cohort_prefixes)

id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)

id_df_wide$count_ones <- rowSums(id_df_wide[, c(3:7)])
df_filtered <- id_df_wide[id_df_wide$count_ones >= 3, ]
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
shared_genera <- df_filtered$feature


# Create a named vector of colors using the phylum color vector

select_cols <- c("Firmicutes"="#aa0000ff", "Bacteroidetes"="#800080ff","Actinobacteria"="#008000ff",
                 "Bacteria_unclassified"="black", "Candidatus_Saccharibacteria"="#808000ff","Proteobacteria"="#00ffffff")
seecolor::print_color(select_cols)
phylum_colors <- select_cols
names(select_cols)


color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

ucla_o_shotgun_species <- result2$plot+
  theme(axis.text.y = element_text(colour = color_mapping))+
  theme(legend.position = "right")
ucla_o_shotgun_species

### Luminal --- 
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
cols <- cols_general
ucla_o_DAT_lum <- generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                    "UCLA O. SPF Luminal",
                                    cols)
hum_v_DAT_lum <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                      "HUM MD Gavage Luminal",
                                                      cols)

cs_DAT_lum <- generate_interregional_taxa_barplot_SITE("CS_SPF/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "CS SPF Luminal",
                                                      cols)

spf_gavage_DAT_lum <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "SPF Gavage Luminal",
                                                  cols)

hum_gavage_DAT_lum <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                    "HUM SD Gavage Luminal",
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

cs_DAT_muc<- generate_interregional_taxa_barplot_SITE("CS_SPF/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "CS SPF Mucosal",
                                                  cols)


hum_gavage_DAT_muc <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "HUM SD Gavage Mucosal",
                                                          cols)

hum_v_gavage_DAT_muc <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                               "HUM MD Gavage Mucosal",
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
