library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)


setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here()

generate_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
  data<-read.csv(ordination_file, header=FALSE)
  metadata <- read.csv(metadata, header=TRUE)
  metadata$SampleID <- metadata$X.SampleID
  
  #store PC1 and Pc2
  PC1<-data[5,1]
  PC1 <-round(as.numeric(PC1)*100, digits=1)
  PC2<-data[5,2]
  PC2 <-round(as.numeric(PC2)*100, digits=1)
  PC1 <-as.character(PC1)
  str_PC1<-paste("PC 1 (", PC1,"%)")
  str_PC2<-paste("PC 2 (", PC2, "%)")
  
  #edit dataframe
  data<-data[,1:4]
  data <- slice(data, 1:(n() - 4))     # Apply slice & n functions
  data<-data[-c(1,2,3,4,5,6,7,8,9),]
  
  #rename columns
  names(data)[names(data) == "V1"] <- "SampleID" 
  names(data)[names(data)=="V2"] <- "PC1" 
  names(data)[names(data)=="V3"] <- "PC2"
  names(data)[names(data)=="V4"] <- "PC3"
  # data$SampleID<-gsub(".","",data$SampleID)
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #declare factors
  data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))

  p<- ggplot(data, aes(x=PC1, y=PC2, colour={{colorvariable}})) + 
    geom_point(size=3) + 
    scale_colour_manual(values={{colorvector}}) +
    #scale_color_viridis_d()+
    xlab(str_PC1) +
    ylab(str_PC2) +
    theme_cowplot(12)+
    theme(legend.position="top") +
    #coord_fixed(ratio=1/2)+
    labs(title= paste0({{title}}, " RPCA")) 
  p
}
metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", header=TRUE)
names(metadata)
Colon_cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
SI_cols<- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
all_cols <-  c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green","Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Type_cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")

lum<- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility_Luminal", Site_General,cols_general)
muc <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility_Mucosal", Site_General,cols_general)
lc <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumCol", Site,Colon_cols)
lsi <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumSI", Site,SI_cols)
mc <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility MucCol", Site,Colon_cols)
msi <-generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility MucSI", Site,SI_cols)

dev.new(width=15, height=10) 
plot_grid(lum,muc, lsi,lc,msi,mc, align="hv", ncol=2)

colon <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Colon", Type, Type_cols)
si <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility SI", Type, Type_cols)
dev.new(width=15, height=10) 
plot_grid(colon,si, align ="hv")

DC <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Distal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility DC", Type, Type_cols)
PC <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Proximal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility PC", Type, Type_cols)
cec <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Cecum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Cec", Type, Type_cols)
ile <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Ileum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Ile", Type, Type_cols)
jej <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Jejunum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Jej", Type, Type_cols)
duo <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Duodenum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Duo", Type, Type_cols)
dev.new(width=15, height=10) 
plot_grid(duo,jej,ile,cec,PC,DC, ncol=2, align ="hv")

### Run Adonis on Site subset ---
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")

metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", header=TRUE, row.names=1)


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

data.adonis=adonis2(data.dist ~ Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
print(data.adonis)

run_repeated_PERMANOVA <- function(path_to_distance_matrix_tsv,path_to_metadata_csv,permute_columns_vector, subject_metadata_vector){
  #data<-read.table(file ="CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
  #metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", header=TRUE, row.names=1)
  
  # Read in files ---
  data<-read.table(path_to_distance_matrix_tsv)
  metadata <- read.csv(path_to_metadata_csv, header=T, row.names=1)
  
  # Ensure metadata matches sample order in distance matrix
  data.dist <- as.dist(as(data, "matrix"))
  target <- row.names(data)
  metadata <- metadata[match(target, row.names(metadata)),]
  target == row.names(metadata)
  
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
                            metadata_order=order_vector)
  print(data.adonis$aov.tab)
}



permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Mucosal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Luminal SI 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_SI_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_Colon_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

# Luminal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "CS-Facility-Analysis/RPCA/rpca_dm/dm_rpca_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "CS-Facility-Analysis/CS_Facility_Metadata.csv",
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
