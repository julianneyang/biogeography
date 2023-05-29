library(GGally)
library(ggplot2)
library(car)
library(vegan)
here::i_am("MouseBiogeography-RProj/Regional-BetaDiversity-RPCA.R")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis")
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_LumCol_PcoA.csv", row.names=1)
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv", row.names=1)
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_MucSI_PcoA.csv", row.names=1)
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_LumSI_PcoA.csv", row.names=1)
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_Luminal_PcoA.csv", row.names=1)
data <- read.csv("RPCA-PCoA/RPCA for all Sites - RPCA_Mucosal_PcoA.csv", row.names=1)

mytheme <- theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                 panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                 panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                 #panel.border=element_blank(), #gets rid of square going around the entire graph
                 axis.line.x = element_line(colour = 'black', size = 0.5),#sets the axis line size
                 axis.line.y = element_line(colour = 'black', size = 0.5),#sets the axis line size
                 axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
                 axis.title.x = element_text(size=16, color="black"), #size of x-axis title
                 axis.title.y = element_text(size=16, color="black"), #size of y-axis title
                 axis.text.x = element_text(size=16, color="black"), #size of x-axis text
                 axis.text.y = element_text(size=16, color="black"))#size of y-axis text

#double check levels
data$Sex <- factor(data$Sex)
data$Sequencing_Run <- factor(data$Sequencing_Run)
data$Site <- factor(data$Site.1)
data$Line <- factor(data$Line)
data$MouseID_Line <- factor(data$MouseID_Line)
data$Genotype <- factor(data$Genotype)
data$Site_General<- factor(data$Site_General, levels =c("SI", "Colon"))
sapply(data,levels)

###Luminal
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  #scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Luminal Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="RPCA_DM/Luminal-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site_General, data=metadata, permutations=10000)
data.adonis

# dbRDA 
cols <- c("SI" = "#F8766D","Colon" ="#00BFC4")
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site_General,data)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site_General)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Mucosal
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  #scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Mucosal Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="RPCA_DM/Mucosal-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site_General, data=metadata, permutations=10000)
data.adonis

# dbRDA 
cols <- c("SI" = "#F8766D","Colon" ="#00BFC4")
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site_General,data)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site_General)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Luminal Colon
cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site.1)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  theme(legend.position="none") 
p + labs(title="LuminalColon Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="RPCA-Distance-Matrices/RPCA_DM_LumCol/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 16S 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
plot(cap_16s, type = "n")
col.list=cols
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19)
ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
plot(ef_16s)


###Mucosal Colon
cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site.1)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  theme(legend.position="none") 
p + labs(title="MucosalColon Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data<-read.table(file ="RPCA_DM/MucCol-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv")
data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 16S Mucosal colon
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)


###Mucosal SI
cols <- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site.1)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="none") 
p + labs(title="MucosalSI Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data<-read.table(file ="RPCA_DM/MucSI-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv")
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)
data.dist <- as.dist(as(data, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Luminal SI
cols <- c("Duodenum" = "red", "Jejunum" = "gold", "Ileum" = "green")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  theme(legend.position="none") 
p + labs(title="LuminalSI Site Differences RPCA")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="RPCA_DM/LumSI-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 16S 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

### Repeated Measures PERMANOVA 
permute_within <- c("Site_General")
subject_data <- c("Sequencing_Run", "Line","Sex", "MouseID")


#Luminal
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/Luminal-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_Luminal_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

#Mucosal 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/Mucosal-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_Mucosal_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Line", "Sex", "MouseID")

# Mucosal Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/MucCol-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Mucosal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/MucSI-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_MucSI_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal Colon 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/LumCol-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_LumCol_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Luminal SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_DM/LumSI-Combat-RPCA-distance-matrix.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA-PCoA/RPCA for all Sites - RPCA_LumSI_PcoA.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

