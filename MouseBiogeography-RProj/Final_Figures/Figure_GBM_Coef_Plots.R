###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 5.29.23
### Figure Number: Main Figure GBM  
### Figure Contents: GBM Luminal and Mucosal Intraregional Coef plots
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(Microbiome.Biogeography)
library(paletteer)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")

### Wrangle GBM data for target features ---
process_results_files <- function(file_paths, feature_value, new_value, new_coef, cohort_prefixes) {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(file_path, header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site" & feature == feature_value)
    
    # Add a new row to the data frame
    new_row <- list(metadata = "Site",
                    feature = feature_value,
                    value = new_value,
                    coef = new_coef,
                    stderr = 0,
                    N = NA,
                    N.not.0 = NA,
                    pval = NA,
                    qval = 100)
    data <- rbind(data, new_row)
    
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

# Define the file paths, cohort prefixes, and other parameters
file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WT_Val_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"
)

new_value <- "Distal_Colon"
new_coef <- 0
cohort_prefixes <- c("UCLA_O_SPF_Lum", 
                     "UCLA_O_SPF_Muc",
                     "CS_SPF_Muc",
                     "CS_SPF_Lum",
                     "UCLA_V_SPF_Muc",
                     "HUM_Gavage_Muc",
                     "HUM_Gavage_Lum",
                     "SPF_Gavage_Lum",
                     "SPF_Gavage_Muc")

file_paths_3 <- c(
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"
)

new_value <- "Distal_Colon"
new_coef <- 0
cohort_prefixes_3 <- c("HUM_Gavage_Muc",
                     "HUM_Gavage_Lum",
                     "SPF_Gavage_Lum",
                     "SPF_Gavage_Muc")


file_paths_2 <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WT_Val_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"
  )

new_value <- "Distal_Colon"
new_coef <- 0
cohort_prefixes_2 <- c("UCLA_O_SPF_Lum", 
                     "UCLA_O_SPF_Muc",
                     "CS_SPF_Muc",
                     "CS_SPF_Lum",
                     "UCLA_V_SPF_Muc")

# Process the results files --
# GABA degradation
feature_value <- "MGB019"
data_all_2 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_2 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_2 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# 17 B- estradiol degradation
feature_value <- "MGB031"
data_all_3 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_3 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_3 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# NO degradation I
feature_value <- "MGB027"
data_all_4 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_4 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_4 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Tryptophan degradation
feature_value <- "MGB049"
data_all_5 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_5 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_5 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Butyrate Synthesis I 
feature_value <- "MGB052"
data_all_6 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_6 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_6 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Plot the combined data frame
paletteer::palettes_d_names
my_palette <- paletteer_d("basetheme::brutal",9)
data_all_2$Cohort <- factor(data_all_2$Cohort)
names(my_palette) <-levels(data_all_2$Cohort)

plot_data <- function(data, title) {
  # Update 'value' column with labels
  data$value <- revalue(data$value,c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value <- factor(data$value, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  # Plotting
  plot <- ggplot(data, aes(x = value, y = coef, group = Cohort, color = Cohort)) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.1) +
    labs(x = "Value", y = "Coefficient") +
    scale_color_manual(values = my_palette, name="") +
    theme_cowplot(16) +
    ggtitle(title) +
    geom_point(size = 3, shape = 21, fill = "black") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}

gaba <- plot_data(data_all_2, "GABA degradation") 

gaba_legend <- plot_data(data_all_2, "GABA degradation") + 
  theme(legend.position = "right")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(1, 1, 10, 10)) 
legend <- cowplot::get_legend(gaba_legend)
grid::grid.newpage()
grid::grid.draw(legend)

estradiol <- plot_data(data_all_3, "17-B-Estradiol degradation" )

no_degradation <- plot_data(data_all_4, "Nitric oxide degradation I (NO dioxygenase)" )

trp_degradation <- plot_data(data_all_5, "Tryptophan degradation" )

butyrate_synthesis <- plot_data(data_all_6, "Butyrate synthesis" )

plot_grid(gaba, estradiol, no_degradation,
          trp_degradation, butyrate_synthesis, 
          labels=c("A","B","C","D","E"),ncol = 2,label_size = 20)

