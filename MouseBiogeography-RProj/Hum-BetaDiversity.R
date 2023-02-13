library(ggplot2)
library(vegan)
library(dplyr)
library(rlang)
library(cowplot)
library(viridis)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

generate_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
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
  
  #declare factors
  data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))
  data$Microbiota <-factor(data$Microbiota, levels=c("Humanized", "Cedars_SPF"))
  
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
metadata <- read.table("Humanized Metadata.tsv.txt", sep="\t", header=TRUE)
metadata$Microbiota
names(metadata)
Colon_cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
SI_cols<- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
all_cols <-  c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green","Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Microbiota_cols <-c("Humanized"="purple", "Cedars_SPF" = "turquoise")

##full dataset RPCA
fulldataset <- generate_pcoA_plots("Humanized RPCA Full - ordination.csv", "Humanized Metadata.tsv.txt", "Humanized vs SPF", Site_General,cols_general)
a<- fulldataset +facet_grid(~Microbiota)
d <- generate_pcoA_plots("Humanized RPCA Full - ordination.csv", "Humanized Metadata.tsv.txt", "Humanized vs SPF", Site,all_cols)
b<- d + facet_grid(~Microbiota)

dev.new(width=40, height=20)
plot_grid(fulldataset, a, d,b,cols=2)

fulldataset <- generate_pcoA_plots("Humanized RPCA Full - ordination.csv", "Humanized Metadata.tsv.txt", "Humanized vs SPF", Type, Type_cols)
a<- fulldataset +facet_grid(~Microbiota)

dev.new(width=40, height=20)
plot_grid(fulldataset, a, cols=2)

dev.new(width=40, height=20)
generate_pcoA_plots("Humanized RPCA Full - ordination.csv", "Humanized Metadata.tsv.txt", "Humanized vs SPF", Microbiota,Microbiota_cols)


##Site RPCA 

lumcol <- generate_pcoA_plots("Site RPCA/Site RPCA - LumCol.csv", "Humanized Metadata.tsv.txt", "Colonic Luminal", Site,Colon_cols)
lumsi <- generate_pcoA_plots("Site RPCA/Site RPCA - LumSI.csv", "Humanized Metadata.tsv.txt", "SI Luminal", Site,SI_cols)
muccol <- generate_pcoA_plots("Site RPCA/Site RPCA - MucCol.csv", "Humanized Metadata.tsv.txt", "Colonic Mucosal", Site,Colon_cols)
mucsi <- generate_pcoA_plots("Site RPCA/Site RPCA - MucSI.csv", "Humanized Metadata.tsv.txt", "SI Mucosal", Site,SI_cols)

lumcol <- generate_pcoA_plots("Site RPCA/Site RPCA - LumCol.csv", "Humanized Metadata.tsv.txt", "Colonic Luminal", Microbiota,Microbiota_cols)
lumsi <- generate_pcoA_plots("Site RPCA/Site RPCA - LumSI.csv", "Humanized Metadata.tsv.txt", "SI Luminal", Microbiota,Microbiota_cols)
muccol <- generate_pcoA_plots("Site RPCA/Site RPCA - MucCol.csv", "Humanized Metadata.tsv.txt", "Colonic Mucosal", Microbiota,Microbiota_cols)
mucsi <- generate_pcoA_plots("Site RPCA/Site RPCA - MucSI.csv", "Humanized Metadata.tsv.txt", "SI Mucosal", Microbiota,Microbiota_cols)

dev.new(width=40, height=20)
plot_grid(lumcol,lumsi,muccol,mucsi, cols=2)

##Source RPCA
humanized <- generate_pcoA_plots("Source RPCA/Source RPCA - Humanized.csv", "Humanized Metadata.tsv.txt", "Humanized Mice", Site_General,cols_general)
spf <- generate_pcoA_plots("Source RPCA/Source RPCA - Cedars_SPF.csv", "Humanized Metadata.tsv.txt", "SPF Mice", Site_General,cols_general)

dev.new(width=40, height=20)
plot_grid(humanized,spf, align="hv")

humanized <- generate_pcoA_plots("Source RPCA/Source RPCA - Humanized.csv", "Humanized Metadata.tsv.txt", "Humanized Mice", Site,all_cols)
spf <- generate_pcoA_plots("Source RPCA/Source RPCA - Cedars_SPF.csv", "Humanized Metadata.tsv.txt", "SPF Mice", Site,all_cols)

dev.new(width=40, height=20)
plot_grid(humanized,spf, align="hv")

humanized <- generate_pcoA_plots("Source RPCA/Source RPCA - Humanized.csv", "Humanized Metadata.tsv.txt", "Humanized Mice", Type,Type_cols)
spf <- generate_pcoA_plots("Source RPCA/Source RPCA - Cedars_SPF.csv", "Humanized Metadata.tsv.txt", "SPF Mice", Type,Type_cols)

dev.new(width=40, height=20)
plot_grid(humanized,spf, align="hv")

##Type RPCA
Type_cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
Colon <- generate_pcoA_plots("Type RPCA/Type RPCA - Colon.csv", "Humanized Metadata.tsv.txt", "Colon : Mucosal vs Luminal", Type, Type_cols)
SI <- generate_pcoA_plots("Type RPCA/Type RPCA - SI.csv", "Humanized Metadata.tsv.txt", "SI : Mucosal vs Luminal", Type, Type_cols)

dev.new(width=40, height=20)
plot_grid(SI, Colon)

Colon <- generate_pcoA_plots("Type RPCA/Type RPCA - Colon.csv", "Humanized Metadata.tsv.txt", "Colon : Mucosal vs Luminal", Microbiota,Microbiota_cols)
SI <- generate_pcoA_plots("Type RPCA/Type RPCA - SI.csv", "Humanized Metadata.tsv.txt", "SI : Mucosal vs Luminal", Microbiota,Microbiota_cols)

dev.new(width=40, height=20)
plot_grid(SI, Colon)

duo <- generate_pcoA_plots("Type RPCA/Type RPCA - Duodenum.csv", "Humanized Metadata.tsv.txt", "Duodenum : Mucosal vs Luminal", Type, Type_cols)
jej <- generate_pcoA_plots("Type RPCA/Type RPCA - Jejunum.csv", "Humanized Metadata.tsv.txt", "Jejunum : Mucosal vs Luminal", Type, Type_cols)
ile <- generate_pcoA_plots("Type RPCA/Type RPCA - Ileum.csv", "Humanized Metadata.tsv.txt", "Ileum : Mucosal vs Luminal", Type, Type_cols)
cec <- generate_pcoA_plots("Type RPCA/Type RPCA - Cecum.csv", "Humanized Metadata.tsv.txt", "Cecum : Mucosal vs Luminal", Type, Type_cols)
PC <- generate_pcoA_plots("Type RPCA/Type RPCA - PC.csv", "Humanized Metadata.tsv.txt", "PC : Mucosal vs Luminal", Type, Type_cols)
DC <- generate_pcoA_plots("Type RPCA/Type RPCA - DC.csv", "Humanized Metadata.tsv.txt", "DC : Mucosal vs Luminal", Type, Type_cols)

dev.new(width=40, height=20)
plot_grid(duo,cec,jej,PC,ile,DC,ncol=2, nrow=3)

duo <- generate_pcoA_plots("Type RPCA/Type RPCA - Duodenum.csv", "Humanized Metadata.tsv.txt", "Duodenum : Mucosal vs Luminal", Microbiota,Microbiota_cols)
jej <- generate_pcoA_plots("Type RPCA/Type RPCA - Jejunum.csv", "Humanized Metadata.tsv.txt", "Jejunum : Mucosal vs Luminal", Microbiota,Microbiota_cols)
ile <- generate_pcoA_plots("Type RPCA/Type RPCA - Ileum.csv", "Humanized Metadata.tsv.txt", "Ileum : Mucosal vs Luminal", Microbiota,Microbiota_cols)
cec <- generate_pcoA_plots("Type RPCA/Type RPCA - Cecum.csv", "Humanized Metadata.tsv.txt", "Cecum : Mucosal vs Luminal", Microbiota,Microbiota_cols)
PC <- generate_pcoA_plots("Type RPCA/Type RPCA - PC.csv", "Humanized Metadata.tsv.txt", "PC : Mucosal vs Luminal", Microbiota,Microbiota_cols)
DC <- generate_pcoA_plots("Type RPCA/Type RPCA - DC.csv", "Humanized Metadata.tsv.txt", "DC : Mucosal vs Luminal", Microbiota,Microbiota_cols)

dev.new(width=40, height=20)
plot_grid(duo,cec,jej,PC,ile,DC,ncol=2, nrow=3)


##Run Adonis on full dataset
data.dist<-read.table(file ="dm_rpca_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)


target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Type + Microbiota*Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Type + Microbiota*Site, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site + Type*Microbiota, data=metadata, permutations=10000)
data.adonis$aov.tab

##Run Adonis on Site subset
data.dist<-read.table(file ="Site RPCA/dm_rpca_LumCol_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Site RPCA/dm_rpca_LumSI_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Site RPCA/dm_rpca_MucCol_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Site RPCA/dm_rpca_MucSI_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")

metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)


target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Microbiota*Site, data=metadata, permutations=10000)
data.adonis$aov.tab

##Run Adonis on Type subset
data.dist<-read.table(file ="Type RPCA/dm_rpca_Colon_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_SI_min10k_Humanized_Combat.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Distal_Colon_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Proximal_Colon_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Cecum_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Ileum_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Jejunum_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Type RPCA/dm_rpca_Duodenum_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")

metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)


target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site+ Microbiota*Type, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex +  Microbiota*Type, data=metadata, permutations=10000)
data.adonis$aov.tab

##Run Adonis on Source subset
data.dist<-read.table(file ="Source RPCA/dm_rpca_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")
data.dist<-read.table(file ="Source RPCA/dm_rpca_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza.txt/distance-matrix.tsv")


metadata <- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv", header=TRUE, row.names=1)


target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site+ Type, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Type + Site, data=metadata, permutations=10000)
data.adonis$aov.tab
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Type + Site_General, data=metadata, permutations=10000)
data.adonis$aov.tab
