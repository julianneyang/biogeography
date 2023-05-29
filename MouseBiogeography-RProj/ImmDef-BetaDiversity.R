library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)

here::i_am("MouseBiogeography-RProj/ImmDef-BetaDiversity.R")
generate_pcoA_plots <- function(ordination_file, metadata, title, colorvector, colorvariable){
  data<-read.csv(ordination_file, header=FALSE)
  metadata <- read.table(metadata, sep="\t", header=TRUE)

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
  
  
  p<- ggplot(data, aes(x=PC1, y=PC2, colour={{colorvariable}})) + 
    mytheme + geom_point(size=3) + 
    scale_colour_manual(name="Site", values={{colorvector}}) +
    #scale_color_viridis_d()+
     xlab(str_PC1) +
    ylab(str_PC2) +
    theme(legend.position="top") 
  p+labs(title= paste0(title, "RPCA"), envir=parent.frame())
  
  
}

names(metadata)
Colon_cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")

SI_cols<- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
all_cols <-  c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green","Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
dev.new(width=15, height=10)
generate_pcoA_plots("WTCohort Site RPCA Post Combat Seq - Mucosal.csv", "WTCohort-Metadata.tsv.txt", "Mucosal", cols_general, Site_General)

dev.new(width=15, height=10)
generate_pcoA_plots("WTCohort Site RPCA Post Combat Seq - Mucosal_Colon.csv", "WTCohort-Metadata.tsv.txt", "Mucosal Colon", Colon_cols, Site)
dev.new(width=15, height=10)
generate_pcoA_plots("WTCohort Site RPCA Post Combat Seq - Mucosal_SI.csv", "WTCohort-Metadata.tsv.txt", "Mucosal SI", SI_cols, Site)


#Run Adonis --- Post-Combat
data.dist<-read.table(file ="ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_Colon_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_SI_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv")

data.dist<-read.table(file ="ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/dm_rpca_Mucosal_SI_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv")

metadata <- read.table("ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.tsv.txt", sep="\t", header=TRUE, row.names=1)
metadata$MouseID <- metadata$MouseID_Original
write.csv(metadata, "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv")

target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex*Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex*Site, data=metadata, permutations=10000)
data.adonis$aov.tab

## Use Repeated Measures PERMANOVA --- Pre-Combat 
permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

#Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/dm_rpca_Mucosal_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_Colon_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_SI_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

## Use Repeated Measures PERMANOVA --- Post-Combat \
#Mucosal
permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

#Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/dm_rpca_Mucosal_min10000_WTCohort-ImmDef-ComBat-Adjusted-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")
# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/dm_rpca_Mucosal_Colon_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/dm_rpca_Mucosal_SI_min10000_WTCohort-ImmDef-ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
