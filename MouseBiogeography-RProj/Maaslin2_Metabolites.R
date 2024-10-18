library(Maaslin2)
library(dplyr)
library(tidyr)
library(here)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Maaslin2")

here::i_am("MouseBiogeography-RProj/Maaslin2_L2_SITE_all.R")
here::here()

### Maaslin2 Function ---
run_maaslin <- function(input_data_path, metadata_path, output_path, 
                        fixed_effects, random_effects, metadata_factors,
                        min_prevalence = 0.15, normalization = "none", transform = "log",
                        reference = NULL) {
  # Load and prepare the input data
  input_data <- read.delim(here(input_data_path), header = TRUE, row.names = 1) 
  df_input_data <- as.data.frame(t(input_data))
  names(df_input_data) <- gsub("X", "", names(df_input_data))
  
  # Load and prepare the metadata
  input_metadata <- read.delim(here(metadata_path), header = TRUE, row.names = 1)
  target <- colnames(df_input_data)
  input_metadata <- input_metadata[match(target, row.names(input_metadata)),]
  
  # Convert necessary columns to factors
  for (col in metadata_factors) {
    input_metadata[[col]] <- factor(input_metadata[[col]])
  }
  
  # Run Maaslin2
  fit_data <- Maaslin2(
    input_data = df_input_data, 
    input_metadata = as.data.frame(input_metadata), 
    output = here(output_path), 
    fixed_effects = fixed_effects, 
    random_effects = random_effects,
    normalization = normalization, 
    transform = transform,
    min_prevalence = min_prevalence,
    plot_heatmap = FALSE, 
    plot_scatter = FALSE,
    reference = reference
  )
  
  return(fit_data)
}

### HUM MD Gavage Dataset ---
run_maaslin(
  input_data_path = "Donors-Analysis/melonnpan/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Donors-Analysis/starting_files/Donors_Metadata.tsv", 
  output_path = "Donors-Analysis/melonnpan/Mets-SeqRunSexSite_General-1-MsID-DonorID",
  fixed_effects = c("Sequencing_Run", "Sex", "Site_General"),
  random_effects = c("MouseID", "DonorID"),
  metadata_factors = c("Sequencing_Run", "Donor_ID", "MouseID", "Sex", "Site", "Site_General"),
  reference = c('Sequencing_Run,Jan_2017', 'Site_General,Colon')
)

### UCLA O SPF Dataset ---
run_maaslin(
  input_data_path = "Regional-Mouse-Biogeography-Analysis/melonnpan/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv", 
  output_path = "Regional-Mouse-Biogeography-Analysis/melonnpan/Mets-SeqRunLineSexSite_General_1-MouseID/",
  fixed_effects = c("Sequencing_Run", "Line", "Sex", "Site_General"),
  random_effects = c("MouseID_Line"),
  metadata_factors = c("Sequencing_Run", "Line", "MouseID_Line", "Sex", "Type", "Site_General"),
  reference =  c("Sequencing_Run,Hiseq_April_Nineteen","Line,ItgCre")
)

### HUM SD Gavage Dataset  ---
run_maaslin(
  input_data_path = "Humanized-Biogeography-Analysis/melonnpan/HUM_SD_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv", 
  output_path = "Humanized-Biogeography-Analysis/melonnpan/HUM_SD_Gavage/Mets-SeqRunSexSite_General-1-MsID/",
  fixed_effects = c("Sequencing_Run", "Sex", "Site_General"),
  random_effects = c("MouseID"),
  metadata_factors = c("Sequencing_Run", "MouseID", "Sex", "Site_General"),
  reference = NULL
)

### SPF Gavage Dataset ---
run_maaslin(
  input_data_path = "Humanized-Biogeography-Analysis/melonnpan/SPF_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv", 
  output_path = "Humanized-Biogeography-Analysis/melonnpan/SPF_Gavage/Mets-SeqRunSexSite_General-1-MsID/",
  fixed_effects = c("Sequencing_Run", "Sex", "Site_General"),
  random_effects = c("MouseID"),
  metadata_factors = c("Sequencing_Run", "MouseID", "Sex", "Site_General"),
  reference = NULL
)

### CS SPF Dataset ---
run_maaslin(
  input_data_path = "CS_SPF/melonnpan/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "CS_SPF/starting_files/CS_Facility_Metadata.tsv", 
  output_path = "CS_SPF/melonnpan/Mets-SeqRunSexSite_General-1-MsID/",
  fixed_effects = c("Sequencing_Run", "Sex", "Site_General"),
  random_effects = c("MouseID"),
  metadata_factors = c("Sequencing_Run", "MouseID", "Sex", "Site_General"),
  reference = NULL
)

### Shotgun Dataset ---
run_maaslin(
  input_data_path = "Shotgun/melonnpan/UCLA_O_SPF/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/UCLA_O_SPF/Mets-LineSexSite_General-1-MsID/",
  fixed_effects = c("Line","Sex", "Site"), 
  random_effects = c("MouseID"),
  reference= c('Line,JJWT','Site,Distal_Colon'),
  metadata_factors = c("Line", "MouseID", "Sex", "Site"),
)

run_maaslin(
  input_data_path = "Shotgun/melonnpan/CS_SPF/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/CS_SPF/Mets-SexSite-1-MsID/",
  fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
  reference= c('Site,Distal_Colon'),
  metadata_factors = c("MouseID", "Sex", "Site"),
)

run_maaslin(
  input_data_path = "Shotgun/melonnpan/SPF_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/",
  fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
  reference= c('Site,Distal_Colon'),
  metadata_factors = c("MouseID", "Sex", "Site"),
)

run_maaslin(
  input_data_path = "Shotgun/melonnpan/SPF_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/SPF_Gavage/Mets-Site-1-MsID/",
  fixed_effects = c("Site"), random_effects = c("MouseID"),
  reference= c('Site,Distal_Colon'),
  metadata_factors = c("MouseID", "Sex", "Site"),
)

run_maaslin(
  input_data_path = "Shotgun/melonnpan/HUM_SD_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/HUM_SD_Gavage/Mets-SexSite-1-MsID/",
  fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
  reference= c('Site,Distal_Colon'),
  metadata_factors = c("MouseID", "Sex", "Site"),
)

run_maaslin(
  input_data_path = "Shotgun/melonnpan/HUM_SD_Gavage/MelonnPan_Predicted_Metabolites.txt", 
  metadata_path = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv", 
  output_path = "Shotgun/melonnpan/HUM_SD_Gavage/Mets-Site-1-MsID/",
  fixed_effects = c("Site"), random_effects = c("MouseID"),
  reference= c('Site,Distal_Colon'),
  metadata_factors = c("MouseID", "Sex", "Site"),
)


