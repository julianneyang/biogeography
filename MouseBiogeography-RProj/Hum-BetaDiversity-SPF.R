library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

metadata <- read.table("Humanized Metadata.tsv.txt", sep="\t", header=TRUE)
names(metadata)
Colon_cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
SI_cols<- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
all_cols <-  c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green","Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Microbiota_cols <-c("Humanized"="purple", "Cedars_SPF" = "turquoise")
Type_cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")

lum<- generate_pcoA_plots("Source RPCA/SPF/Site/Source RPCA - SPF - Lum.csv", "Humanized Metadata.tsv.txt", "SPF Luminal", Site_General,cols_general)
muc <- generate_pcoA_plots("Source RPCA/SPF/Site/Source RPCA - SPF - Muc.csv", "Humanized Metadata.tsv.txt", "SPF Mucosal", Site_General,cols_general)
lc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LC.csv", "Humanized Metadata.tsv.txt", "SPF LumCol", Site,Colon_cols)
lsi <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LSI.csv", "Humanized Metadata.tsv.txt", "SPF LumSI", Site,SI_cols)
mc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MC.csv", "Humanized Metadata.tsv.txt", "SPF MucCol", Site,Colon_cols)
msi <-generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MSI.csv", "Humanized Metadata.tsv.txt", "SPF MucSI", Site,SI_cols)

plot_grid(lum,muc, lsi,lc,msi,mc, align="hv", ncol=2)

colon <- generate_pcoA_plots("Source RPCA/SPF/Type/Source RPCA - SPF - Colon.csv", "Humanized Metadata.tsv.txt", "SPF Colon", Type, Type_cols)
si <- generate_pcoA_plots("Source RPCA/SPF/Type/Source RPCA - SPF - SI.csv", "Humanized Metadata.tsv.txt", "SPF SI", Type, Type_cols)
plot_grid(colon,si, align ="hv")

DC <- generate_pcoA_plots("Source RPCA/", "Humanized Metadata.tsv.txt", "Humanized DC", Type, Type_cols)
PC <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - PC.csv", "Humanized Metadata.tsv.txt", "Humanized PC", Type, Type_cols)
cec <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Cec.csv", "Humanized Metadata.tsv.txt", "Humanized Cec", Type, Type_cols)
ile <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Ile.csv", "Humanized Metadata.tsv.txt", "Humanized Ile", Type, Type_cols)
jej <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Jej.csv", "Humanized Metadata.tsv.txt", "Humanized Jej", Type, Type_cols)
duo <- generate_pcoA_plots("Source RPCA/Hum/Type/Source RPCA -Humanized - Duo.csv", "Humanized Metadata.tsv.txt", "Humanized Duo", Type, Type_cols)
plot_grid(duo,jej,ile,cec,PC,DC, ncol=2, align ="hv")

##Run Adonis on Site subset
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Luminal_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Mucosal_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Luminal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Luminal_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Mucosal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Mucosal_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")

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
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/SPF/dm_rpca_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")

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
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Mucosal_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Mucosal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Luminal SI 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Luminal_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Luminal_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Luminal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Luminal_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Mucosal_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

## Run Repeated Measures PERMANOVA - TYPE 

permute_within <- c("Site", "Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)


# SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/dm_rpca_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Type")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Duodenum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Duodenum_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Jejunum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Jejunum_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Ileum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Ileum_SI_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# cecum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Cecum_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Proximal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Proximal_Colon_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Distal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Source RPCA/SPF/Type_RPCA/dm_rpca_Distal_Colon_Colon_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Humanized Metadata - All-Humanized-Metadata (1).csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
