library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

metadata <- read.table("Humanized Metadata.tsv.txt", sep="\t", header=TRUE,row.names = 1)
names(metadata)
metadata <- write.csv(metadata, "Humanized_Metadata.csv")

Colon_cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
SI_cols<- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
all_cols <-  c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green","Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Microbiota_cols <-c("Humanized"="purple", "Cedars_SPF" = "turquoise")
Type_cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")

lum<- generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - Luminal.csv", "Humanized Metadata.tsv.txt", "Humanized Luminal", Site_General,cols_general)
muc <- generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - Mucosal.csv", "Humanized Metadata.tsv.txt", "Humanized Mucosal", Site_General,cols_general)
lc <- generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - LC.csv", "Humanized Metadata.tsv.txt", "Humanized LumCol", Site,Colon_cols)
lsi <- generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - LSI.csv", "Humanized Metadata.tsv.txt", "Humanized LumSI", Site,SI_cols)
mc <- generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - MC.csv", "Humanized Metadata.tsv.txt", "Humanized MucCol", Site,Colon_cols)
msi <-generate_pcoA_plots("Source RPCA/Hum/Site/Source RPCA -Humanized - MSI.csv", "Humanized Metadata.tsv.txt", "Humanized MucSI", Site,SI_cols)

plot_grid(lum,muc, lsi,lc,msi,mc, align="hv", ncol=2)

colon <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Colon.csv", "Humanized Metadata.tsv.txt", "Humanized Colon", Type, Type_cols)
si <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - SI.csv", "Humanized Metadata.tsv.txt", "Humanized SI", Type, Type_cols)
plot_grid(colon,si, align ="hv")

DC <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - DC.csv", "Humanized Metadata.tsv.txt", "Humanized DC", Type, Type_cols)
PC <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - PC.csv", "Humanized Metadata.tsv.txt", "Humanized PC", Type, Type_cols)
cec <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Cec.csv", "Humanized Metadata.tsv.txt", "Humanized Cec", Type, Type_cols)
ile <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Ile.csv", "Humanized Metadata.tsv.txt", "Humanized Ile", Type, Type_cols)
jej <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Jej.csv", "Humanized Metadata.tsv.txt", "Humanized Jej", Type, Type_cols)
duo <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Duo.csv", "Humanized Metadata.tsv.txt", "Humanized Duo", Type, Type_cols)
plot_grid(duo,jej,ile,cec,PC,DC, ncol=2, align ="hv")

##Run Adonis on Site subset
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Luminal_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Mucosal_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Luminal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Luminal_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Mucosal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Mucosal_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")

metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)


target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex*Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis$aov.tab

##Run Adonis on Type subset
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/Hum/dm_rpca_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")

metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Type, data=metadata, permutations=10000)
data.adonis$aov.tab

## Run Repeated Measures PERMANOVA - SITE 
permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Mucosal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Mucosal_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Mucosal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Luminal SI 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Luminal_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Luminal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Luminal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Luminal_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Mucosal_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)


run_repeated_PERMANOVA <- function(path_to_distance_matrix_tsv,path_to_metadata_csv,permute_columns_vector, subject_metadata_vector){
  #data<-read.table(file ="Source RPCA/Hum/dm_rpca_Luminal_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
  #metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)
  
  # Read in files ---
  data<-read.table(path_to_distance_matrix_tsv)
  metadata <- read.csv(path_to_metadata_csv, header=T, row.names=1)
  
  # Ensure metadata matches sample order in distance matrix
  data.dist <- as.dist(as(data, "matrix"))
  target <- row.names(data)
  metadata <- metadata[match(target, row.names(metadata)),]
  target == row.names(metadata)
  
  # Fix metadata columns 
  if("MouseID_Line" %in% names(metadata)){
    metadata$MouseID <- metadata$MouseID_Line
  }
  
  if("Site.1" %in% names(metadata)){
    metadata$Site <- factor(metadata$Site.1)
  }
  if(class(metadata$MouseID)=="integer"){
    metadata$MouseID <- paste("Mouse_",metadata$MouseID)
  }
  
  # Read in relevant metadata, where permute_within (Timepoint, SampleType) and subject data (Age,Sex)
  permute_within <- c(permute_columns_vector)
  subject_data <- c(subject_metadata_vector)
  
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
                                           metadata_order = order_vector)
  print(data.adonis$aov.tab)
}

## Run Repeated Measures PERMANOVA - TYPE 

permute_within <- c("Site", "Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site", "Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_SI_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Duodenum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Duodenum_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Jejunum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Jejunum_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Ileum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Ileum_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Cecum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Cecum_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Proximal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Proximal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Distal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/Hum/dm_rpca_Distal_Colon_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
