library(Maaslin2)
library(funrar)
library(dplyr)
library(here)

here::i_am("MouseBiogeography-RProj/Donors-Maaslin2-TYPE-Genus.R")

perform_maaslin_analysis <- function(tissue_type, file_path) {
  # Construct file path for the input data
  input_data_path <- here::here(file_path, paste0("export_L6_", tissue_type, "_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"))
  
  # Read input data
  input_data <- read.delim(input_data_path, header=TRUE, row.names=1)
  df_input_data <- as.data.frame(input_data)
  df_input_data <- dplyr::select(df_input_data, -c("taxonomy"))
  names(df_input_data) <- gsub("X", "", names(df_input_data))
  
  # Read input metadata
  input_metadata <- read.delim(here::here("Donors-Analysis/starting_files/Donors_Metadata.tsv"), header=TRUE)
  input_metadata$SampleID <- gsub("-", ".", input_metadata$SampleID)
  row.names(input_metadata) <- input_metadata$SampleID
  
  # Match metadata with target
  target <- colnames(df_input_data)
  input_metadata <- input_metadata[match(target, row.names(input_metadata)), ]
  
  # Prepare metadata
  df_input_metadata <- input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
  df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
  df_input_metadata$Sex <- factor(df_input_metadata$Sex)
  df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
  
  # Perform Maaslin2 analysis
  fit_data <- Maaslin2(
    input_data = df_input_data,
    input_metadata = df_input_metadata,
    output = here::here(paste0("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-", tissue_type, "-ComBat-SeqRunSexType-1-MsID-DonorID")),
    fixed_effects = c("Sequencing_Run", "Sex", "Type"),
    random_effects = c("MouseID", "DonorID"),
    normalization = "clr",
    transform = "none",
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    min_prevalence = 0.15,
    reference = c("Sequencing_Run,Jan_2017", "Type,Luminal")
  )
}

tissue_types <- c("Cecum", "Proximal_Colon", "Distal_Colon", "Duodenum", "Jejunum", "Ileum","Colon","SI")
filepath <- "Donors-Analysis/type_subsets/collapsed_taxa/genus_level/"


for (tissue_type in tissue_types) {
  perform_maaslin_analysis(tissue_type, filepath)
}




