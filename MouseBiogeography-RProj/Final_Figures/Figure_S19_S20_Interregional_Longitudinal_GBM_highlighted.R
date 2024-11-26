library(dplyr)
library(here)


library("Microbiome.Biogeography")


here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S19_S20_Interregional_Longitudinal_GBM.R")

### Find shared GBMs ---
lum_file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",  
  "CS_SPF/OMIXER-RPM Results/CS_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")

lum_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "SPF_Gavage",
                         "HUM_SD_Gavage")

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
annotation <- read.csv("GBM_Module_Key.csv")
annotated_df_filtered <- merge(df_filtered, annotation, by="feature")
shared_luminal_GBMs <- annotated_df_filtered$annotation

### Identify shared mucosal GBMs ---
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
muc_file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",  
  "CS_SPF/OMIXER-RPM Results/CS_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
  "Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
  "UCLA_V_SPF_Analysis/OMIXER-RPM/WT_Val_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv")

muc_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "SPF_Gavage",
                         "HUM_SD_Gavage",
                         "UCLA_V_SPF")

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
annotation <- read.csv("GBM_Module_Key.csv")
annotated_df_filtered <- merge(df_filtered, annotation, by="feature")
shared_mucosal_GBMs <- annotated_df_filtered$annotation

### Highlight shared GBMs
highlight_GBM <- function (genera_to_highlight, specific_genera){
  
  firebrick_vector <- rep("firebrick", length(intersect(specific_genera,genera_to_highlight)))
  names(firebrick_vector) <- intersect(specific_genera,genera_to_highlight)
  black_vector <- rep("black", length(setdiff(specific_genera,genera_to_highlight)))
  names(black_vector) <- setdiff(specific_genera,genera_to_highlight)
  combined <- c(firebrick_vector,black_vector)
  
  reordered_vector <- combined[match(specific_genera, names(combined))]
  return(reordered_vector)
  
}

### Luminal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)

hum_v_GBM_Map_result <- generate_interregional_GBM_barplot("Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "HUM MD Gavage Luminal",
                                                         cols)

region_specific_GBM <- hum_v_GBM_Map_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_luminal_GBMs,specific_genera = region_specific_GBM)

hum_v_GBM_Map_lum <- hum_v_GBM_Map_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


ucla_o_GBM_Map_result <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA O. SPF Luminal",
                                                         cols)
region_specific_GBM <- ucla_o_GBM_Map_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_luminal_GBMs,specific_genera = region_specific_GBM)

ucla_o_GBM_Map_DAT_lum <- ucla_o_GBM_Map_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

cs_spf_GBM_Map_result <- generate_interregional_GBM_barplot("CS_SPF/OMIXER-RPM Results/CS_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "CS SPF Luminal",
                                                         cols)

region_specific_GBM <- cs_spf_GBM_Map_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_luminal_GBMs,specific_genera = region_specific_GBM)

cs_spf_GBM_Map_lum <- cs_spf_GBM_Map_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

spf_gavage_GBM_Map_result <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "SPF Gavage Luminal",
                                                             cols)

region_specific_GBM <- spf_gavage_GBM_Map_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_luminal_GBMs,specific_genera = region_specific_GBM)

spf_gavage_GBM_Map_lum <- spf_gavage_GBM_Map_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

hum_gavage_GBM_Map_result <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "HUM SD Gavage Luminal",
                                                             cols)
region_specific_GBM <- hum_gavage_GBM_Map_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_luminal_GBMs,specific_genera = region_specific_GBM)

hum_gavage_GBM_Map_lum <- hum_gavage_GBM_Map_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


dev.new()
cowplot::plot_grid(ucla_o_GBM_Map_DAT_lum, cs_spf_GBM_Map_lum, 
                   nrow=1, ncol=2,
                   labels=c("A", "B"),
                   label_size = 20)

dev.new()
cowplot::plot_grid(spf_gavage_GBM_Map_lum, hum_gavage_GBM_Map_lum, 
                   hum_v_GBM_Map_lum,
                   nrow=2, ncol=2,
                   labels=c("C","D", "E"),
                   label_size = 20)

### Mucosal Barplots: Aggregated by Map (median coef) --- 
cols <- viridis::viridis(2)
hum_v_GBM_Map_muc_result <- generate_interregional_GBM_barplot("Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "HUM MD Gavage Mucosal",
                                                         cols)
region_specific_GBM <- hum_v_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
hum_v_GBM_Map_muc <- hum_v_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


ucla_o_GBM_Map_muc_result <- generate_interregional_GBM_barplot("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA O. SPF Mucosal",
                                                         cols)
region_specific_GBM <-ucla_o_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
ucla_o_GBM_Map_muc <- ucla_o_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


ucla_v_GBM_Map_muc_result <- generate_interregional_GBM_barplot("UCLA_V_SPF_Analysis/OMIXER-RPM/WT_Val_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "UCLA V. SPF Mucosal",
                                                         cols)
region_specific_GBM <-ucla_v_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
ucla_v_GBM_Map_muc <- ucla_v_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

cs_spf_GBM_Map_muc_result <- generate_interregional_GBM_barplot("CS_SPF/OMIXER-RPM Results/CS_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                         "GBM_Module_Key.csv",
                                                         titlestring= "CS SPF Mucosal",
                                                         cols)
region_specific_GBM <-cs_spf_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
cs_spf_GBM_Map_muc<- cs_spf_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


spf_gavage_GBM_Map_muc_result <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "SPF Gavage Mucosal",
                                                             cols)
region_specific_GBM <- spf_gavage_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
spf_gavage_GBM_Map_muc<- spf_gavage_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))

hum_gavage_GBM_Map_muc_result <- generate_interregional_GBM_barplot("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             titlestring= "HUM SD Gavage Mucosal", 
                                                             cols)
region_specific_GBM <- hum_gavage_GBM_Map_muc_result$dataframe
region_specific_GBM <- region_specific_GBM$annotation
reordered_vector <- highlight_GBM(genera_to_highlight = shared_mucosal_GBMs,specific_genera = region_specific_GBM)
hum_gavage_GBM_Map_muc<- hum_gavage_GBM_Map_muc_result$plot + theme(axis.text.y = element_text(colour = reordered_vector))


dev.new()
cowplot::plot_grid(ucla_o_GBM_Map_muc, ucla_v_GBM_Map_muc, 
                   cs_spf_GBM_Map_muc, spf_gavage_GBM_Map_muc,
                  labels=c("A","B","C", "D"), 
                  label_size = 20,
                  nrow=2, ncol=2)

dev.new()
cowplot::plot_grid(hum_gavage_GBM_Map_muc, hum_v_GBM_Map_muc,
                   nrow=1, ncol=2,
                   labels=c("E", "F"),
                   label_size = 20)



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
