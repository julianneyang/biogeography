library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
library(here)

here::i_am("MouseBiogeography-RProj/Donors-BetaDiversity.R")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

run_repeated_PERMANOVA_donors <- function(path_to_distance_matrix_tsv,path_to_metadata_csv,permute_columns_vector, subject_metadata_vector){
  #data<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
  #metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", header=TRUE, row.names=1)
  #data<- read.table( "Donors-Analysis/site_rpca/dm_rpca_SI_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv")
  #metadata <- read.delim("Donors-Analysis/starting_files/Donors_Metadata.tsv",header=T,row.names=1)
 
   # Read in files ---
  data<-read.table(path_to_distance_matrix_tsv)
  row.names(data) <- gsub("-",".",row.names(data))
  metadata <- read.delim(here(path_to_metadata_csv), header=T, row.names=1)
  row.names(metadata) <- gsub("-",".",row.names(metadata))
                              
  # Ensure metadata matches sample order in distance matrix
  data.dist <- as.dist(as(data, "matrix"))
  target <- row.names(data)
  metadata <- metadata[match(target, row.names(metadata)),]
  target == row.names(metadata)
  
  # Read in relevant metadata, where permute_within (Timepoint, SampleType) and subject data (Age,Sex)
  permute_within <- c(permute_columns_vector)
  subject_data <- c(subject_metadata_vector)
  #permute_within <- c("Site")
  #subject_data <- c("Sequencing_Run", "Sex", "MouseID")
  
  # Wrangle metadata into appropriate formats 
   general_metadata<- dplyr::select(metadata, c(permute_within))
  metadata_subj <- dplyr::select(metadata, c(subject_data))
  metadata_subj <-as.data.frame(metadata_subj[!duplicated(metadata$MouseID),]) #one of these columns is your SubjectID
  row.names(metadata_subj) <- metadata_subj$MouseID
  metadata_subj <- dplyr::select(metadata_subj, -MouseID)
  
  subjectvector <- c(metadata$MouseID)
  order_vector <- head(subject_data,-1)
  order_vector <- c(order_vector, permute_within)
  
  # Run repeat-measrues aware PERMANOVA (Lloyd-Price et al., 2019)

  data.adonis <- PERMANOVA_repeat_measures(D = data.dist, permutations=10000,
                                           permute_within= general_metadata, 
                                           blocks= subjectvector, 
                                           block_data=metadata_subj,
                                           metadata_order=order_vector)
  data.adonis
}



permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Mucosal SI
run_repeated_PERMANOVA_donors(path_to_distance_matrix_tsv = "Donors-Analysis/site_rpca/dm_rpca_SI_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Donors-Analysis/starting_files/Donors_Metadata.tsv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Mucosal Colon
run_repeated_PERMANOVA_donors(path_to_distance_matrix_tsv = "Donors-Analysis/site_rpca/dm_rpca_Colon_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Donors-Analysis/starting_files/Donors_Metadata.tsv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Luminal SI 
run_repeated_PERMANOVA_donors(path_to_distance_matrix_tsv = "Donors-Analysis/site_rpca/dm_rpca_SI_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Donors-Analysis/starting_files/Donors_Metadata.tsv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal Colon
data<- read.table( "Donors-Analysis/site_rpca/dm_rpca_Colon_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv")
metadata <- read.delim("Donors-Analysis/starting_files/Donors_Metadata.tsv",header=T,row.names=1)
row.names(data) <- gsub("-",".",row.names(data))
row.names(metadata) <- gsub("-",".",row.names(metadata))

data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Donor_ID+Site, data=metadata, permutations=10000)
data.adonis

permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Luminal 
data<- read.table( "Donors-Analysis/site_rpca/dm_rpca_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv")
metadata <- read.delim("Donors-Analysis/starting_files/Donors_Metadata.tsv",header=T,row.names=1)
row.names(data) <- gsub("-",".",row.names(data))
row.names(metadata) <- gsub("-",".",row.names(metadata))

data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Site_General, data=metadata, permutations=10000)
data.adonis

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Donor_ID + Site_General, data=metadata, permutations=10000)
data.adonis

run_repeated_PERMANOVA_donors(path_to_distance_matrix_tsv = "Donors-Analysis/site_rpca/dm_rpca_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Donors-Analysis/starting_files/Donors_Metadata.tsv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal 
data<- read.table( "Donors-Analysis/site_rpca/dm_rpca_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv")
metadata <- read.delim("Donors-Analysis/starting_files/Donors_Metadata.tsv",header=T,row.names=1)
row.names(data) <- gsub("-",".",row.names(data))
row.names(metadata) <- gsub("-",".",row.names(metadata))

data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

metadata$Donor_ID <- factor(metadata$Donor_ID)

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Site_General, data=metadata, permutations=10000)
data.adonis

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Donor_ID + Site_General, data=metadata, permutations=10000)
data.adonis


run_repeated_PERMANOVA_donors(path_to_distance_matrix_tsv = "Donors-Analysis/site_rpca/dm_rpca_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Donors-Analysis/starting_files/Donors_Metadata.tsv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

### Run Adonis on Type subset ---
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")

metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", header=TRUE, row.names=1)

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Type, data=metadata, permutations=10000)
data.adonis$aov.tab

# Repeated Measures
permute_within <- c("Site", "Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Repeated Measures
permute_within <- c("Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Duodenum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Duodenum_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Jejunum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Jejunum_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Ileum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Ileum_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Cecum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Cecum_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)


# Proximal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Proximal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Distal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Distal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
