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

data_all <- data.frame(metadata=character(),
                       feature=character(), 
                       value=character(), 
                       coef=character(),
                       stderr=character(),
                       N=character(),
                       N.not.0=character(),
                       pval=character(),
                       qval=character(),
                       Cohort=character(),
                       stringsAsFactors=FALSE) 
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

df <- as.data.frame(shared_genera)
df$feature <- df$shared_genera
df$Phylum <- gsub(".*\\.p__", "", df$feature)
df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
df$Order <- gsub(".*\\.o__", "", df$feature)
df$Order <- gsub("\\.f__.*", "", df$Order)
df$Order <- paste0(df$Order, " (o)")
df$Family <- gsub(".*\\.f__", "", df$feature)
df$Family <- gsub("\\.g__.*", "", df$Family)
df$Family<- paste0(df$Family, " (f)")
df$Genus <- gsub(".*\\.g__", "", df$feature)
df$Genus <- gsub("\\.g__.*", "", df$Genus)
df$Species <- gsub(".*\\.s__", "", df$feature)

df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))

df <- df %>%
  mutate(level1 = ifelse(nchar(Genus) != 0, Genus, Family))
df <- df %>%
  mutate(annotation = ifelse(level1!= " (f)", level1, Order))
shared_genera <-df$annotation


# Create a named vector of colors using the phylum color vector

ucla_o_spf_result<-generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                             "UCLA O. SPF Luminal",
                                                             cols)


highlight_genera <- function (genera_to_highlight, specific_genera){

firebrick_vector <- rep("firebrick", length(intersect(specific_genera,genera_to_highlight)))
names(firebrick_vector) <- intersect(specific_genera,genera_to_highlight)
black_vector <- rep("black", length(setdiff(specific_genera,genera_to_highlight)))
names(black_vector) <- setdiff(specific_genera,genera_to_highlight)
combined <- c(firebrick_vector,black_vector)

reordered_vector <- combined[match(specific_genera, names(combined))]
return(reordered_vector)

}



### Luminal --- 
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
cols <- cols_general

ucla_o_spf_result<-generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                            "UCLA O. SPF Luminal",
                                                            cols)
region_specific_genus <- ucla_o_spf_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

ucla_o_DAT_lum <- ucla_o_spf_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

hum_v_DAT_result <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                      "HUM MD Gavage Luminal",
                                                      cols)

region_specific_genus <- hum_v_DAT_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

hum_v_DAT_lum <- hum_v_DAT_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))
cs_DAT_result<- generate_interregional_taxa_barplot_SITE("CS_SPF/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "CS SPF Luminal",
                                                      cols)
region_specific_genus <- cs_DAT_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

cs_DAT_lum <- cs_DAT_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

spf_gavage_DAT_result <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "SPF Gavage Luminal",
                                                  cols)
region_specific_genus <- spf_gavage_DAT_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

spf_gavage_DAT_lum <- spf_gavage_DAT_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))
hum_gavage_DAT_result <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                    "HUM SD Gavage Luminal",
                                    cols)

region_specific_genus <- hum_gavage_DAT_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

hum_gavage_DAT_lum <- hum_gavage_DAT_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))
dev.new()
cowplot::plot_grid(ucla_o_DAT_lum, cs_DAT_lum,spf_gavage_DAT_lum, ncol=3,
                   labels=c("A","B","C"))
dev.new()
cowplot::plot_grid(hum_gavage_DAT_lum,hum_v_DAT_lum,ncol=2,
                   labels=c("D","E"))


### Identify overlapping features ---
muc_file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/all_results.tsv",  
  "CS_SPF/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv",
  "UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite_General-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv")

muc_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "UCLA_V_SPF",
                         "HUM_MD_Gavage",
                         "SPF_Gavage",
                         "HUM_SD_Gavage")

all_taxa <- process_results_for_upset_plot_interregional(file_paths = muc_file_paths,
                                                         cohort_prefixes = muc_cohort_prefixes)

id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)

id_df_wide$count_ones <- rowSums(id_df_wide[, c(3:8)])
df_filtered <- id_df_wide[id_df_wide$count_ones >= 3, ]
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
shared_genera <- df_filtered$feature

df <- as.data.frame(shared_genera)
df$feature <- df$shared_genera
df$Phylum <- gsub(".*\\.p__", "", df$feature)
df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
df$Order <- gsub(".*\\.o__", "", df$feature)
df$Order <- gsub("\\.f__.*", "", df$Order)
df$Order <- paste0(df$Order, " (o)")
df$Family <- gsub(".*\\.f__", "", df$feature)
df$Family <- gsub("\\.g__.*", "", df$Family)
df$Family<- paste0(df$Family, " (f)")
df$Genus <- gsub(".*\\.g__", "", df$feature)
df$Genus <- gsub("\\.g__.*", "", df$Genus)
df$Species <- gsub(".*\\.s__", "", df$feature)

df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))

df <- df %>%
  mutate(level1 = ifelse(nchar(Genus) != 0, Genus, Family))
df <- df %>%
  mutate(annotation = ifelse(level1!= " (f)", level1, Order))
shared_genera <-df$annotation


 ### Mucosal --- 
ucla_o_DAT_muc_result<- generate_interregional_taxa_barplot_SITE("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID/significant_results.tsv",
                                                      "UCLA O. SPF Mucosal",
                                                      cols)

region_specific_genus <- ucla_o_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

ucla_o_DAT_muc <- ucla_o_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


ucla_v_DAT_muc_result <- generate_interregional_taxa_barplot_SITE("UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                      "UCLA V. SPF Mucosal",
                                                      cols)

region_specific_genus <- ucla_v_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

ucla_v_DAT_muc <- ucla_v_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

cs_DAT_muc_result<- generate_interregional_taxa_barplot_SITE("CS_SPF/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                  "CS SPF Mucosal",
                                                  cols)
region_specific_genus <- cs_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

cs_DAT_muc <- cs_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


hum_gavage_DAT_muc_result <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "HUM SD Gavage Mucosal",
                                                          cols)

region_specific_genus <- hum_gavage_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

hum_gavage_DAT_muc <- hum_gavage_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


hum_v_gavage_DAT_muc_result <- generate_interregional_taxa_barplot_SITE("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                               "HUM MD Gavage Mucosal",
                                                               cols)
region_specific_genus <- hum_v_gavage_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

hum_v_gavage_DAT_muc <- hum_v_gavage_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


spf_gavage_DAT_muc_result <- generate_interregional_taxa_barplot_SITE("Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                          "SPF Gavage Mucosal",
                                                          cols)

region_specific_genus <- spf_gavage_DAT_muc_result$dataframe
region_specific_genus <- region_specific_genus$annotation
reordered_vector <- highlight_genera(genera_to_highlight = shared_genera,specific_genera = region_specific_genus)

spf_gavage_DAT_muc <- spf_gavage_DAT_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


dev.new()
cowplot::plot_grid(ucla_o_DAT_muc, ucla_v_DAT_muc, cs_DAT_muc,ncol=3,
                   labels=c("A","B","C"),
                   label_size = 16)

dev.new()
cowplot::plot_grid(spf_gavage_DAT_muc, hum_gavage_DAT_muc, hum_v_gavage_DAT_muc, 
                   ncol=3, labels=c("D","E","F"),
                   label_size = 16)
