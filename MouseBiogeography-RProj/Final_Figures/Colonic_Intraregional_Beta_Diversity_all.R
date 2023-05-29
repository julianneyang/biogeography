library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
library(Microbiome.Biogeography)

setwd("C:/Users/Jacobs Laboratory/Desktop/biogeography_github/biogeography/")

### UCLA O. SPF ---
permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Line","Sex", "MouseID")

# Cecum vs PC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "2021-8-Microbiome-Batch-Correction-Analysis/intraregional/dm_rpca_Cec_PC_Mucosal-Combat.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Cecum vs DC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "2021-8-Microbiome-Batch-Correction-Analysis/intraregional/dm_rpca_Cec_DC_Mucosal-Combat.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# PC vs DC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "2021-8-Microbiome-Batch-Correction-Analysis/intraregional/dm_rpca_PC_DC_Mucosal-Combat.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
### UCLA V SPF (no combat) ---

permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Cecum vs PC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/intraregional_nocombat/dm_rpca_Cec_PC_Mucosal_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Cecum vs DC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/intraregional_nocombat/dm_rpca_Cec_DC_Mucosal_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# PC vs DC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/intraregional_nocombat/dm_rpca_PC_DC_Mucosal_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

### CS SPF ---
permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Cecum vs PC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/intraregional/dm_rpca_Cec_PC_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Cecum vs DC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/intraregional/dm_rpca_Cec_DC_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# PC vs DC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/intraregional/dm_rpca_PC_DC_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

### HUM Gavage ---

permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Cecum vs PC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Hum/intraregional/dm_rpca_Cec_PC_Mucosal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Hum/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Cecum vs DC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Hum/intraregional/dm_rpca_Cec_DC_Mucosal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Hum/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# PC vs DC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Hum/intraregional/dm_rpca_PC_DC_Mucosal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Hum/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

## SPF Gavage --- 
permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Cecum vs PC 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "SPF/intraregional/dm_rpca_Cec_PC_Mucosal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "SPF/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Cecum vs DC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "SPF/intraregional/dm_rpca_Cec_DC_Mucosal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "SPF/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# PC vs DC
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "SPF/intraregional/dm_rpca_PC_DC_Mucosal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "SPF/Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
