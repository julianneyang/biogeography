library(ggplot2)
library(vegan)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis")
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Cecum PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Colon PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Duodenum PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Jejunum PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Ileum PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Proximal_Colon PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/Distal_Colon PcoA.csv", header = TRUE)
data <- read.csv("RPCA_LuminalvsMucosal_PcoA/SI PcoA.csv", header = TRUE)

metadata <- read.csv("RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All.csv", header=TRUE)

intermediate<- (merge(data, metadata, by = 'Description'))
data<- intermediate

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
data$Site <- factor(data$Site)
data$Line <- factor(data$Line)
data$MouseID_Line <- factor(data$MouseID_Line)
data$Genotype <- factor(data$Genotype)
data$Type <- factor(data$Type)
data$Site_General<- factor(data$Site_General, levels =c("SI", "Colon"))
sapply(data,levels)

###Duodenum

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Duodenum Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Duodenum_Combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Jejunum

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Jejunum Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Jejunum_Combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type, metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Ileum

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Ileum Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Ileum_Combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Proximal_Colon

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Proximal Colon Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Proximal_Colon_Combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Distal_Colon

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Distal Colon Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Distal_Colon_combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)


###Cecum

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Cecum Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Cecum_Combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

###Colon

cols<-c("Luminal"="#481567FF", Mucosal = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="Colon Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/Colon-combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site + Type, data=metadata, permutations=10000)
data.adonis
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type*Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)


###SI

cols<-c("Luminal"="#481567FF", "Mucosal" = "#3CBB75FF")
dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Type", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  #theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  #theme(axis.text = element_text(colour = "black", size="20"))
  theme(legend.position="top") 
p + labs(title="SI Luminal vs Mucosal RPCA")

#Adonis
metadata<- data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
metadata$Type <- factor(metadata$Type)
row.names(metadata) <- metadata$Description
data.dist<-read.table(file ="RPCA_LuminalvsMucosal_DM/SI-combat_distmatrix_export/distance-matrix.tsv")
target <- row.names(data.dist)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data.dist) == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site + Type, data=metadata, permutations=10000)
data.adonis
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Type*Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 
bray_dist_16s <- data.dist
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site + Type,metadata)
col.list=cols
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(metadata$Type)],pch=19) 

ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site_General")], permu = 10000)  # only show significant covariates
ef_16s
dev.new(width=12, height=10)
plot(cap_16s, type = "n")
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19) 
plot(ef_16s)

### Repeated Measures PERMANOVA- Mucosal vs Luminal ---

metadata <- read.csv("RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All.csv", header=TRUE)
colnames <- names(metadata)

colnames <- colnames[c(2,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17)]
metadata <- metadata %>% select(colnames)
row.names(metadata)<-metadata$Description
#metadata <- write.csv(metadata, "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv")
metadata <-read.csv("RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",row.names=1)

permute_within <- c("Site","Type")
subject_data <- c("Sequencing_Run", "Line","Sex", "MouseID")

# SI
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/SI-combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

# Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Colon-combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

permute_within <- c("Type")
subject_data <- c("Sequencing_Run", "Line","Sex", "MouseID")
# Duodenum 
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Duodenum_Combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Jejunum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Jejunum_Combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Ileum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Ileum_Combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Cecum
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Cecum_Combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Proximal_Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Proximal_Colon_Combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
# Distal_Colon
run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "RPCA_LuminalvsMucosal_DM/Distal_Colon_combat_distmatrix_export/distance-matrix.tsv",
                       path_to_metadata_csv = "RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All-wrangled.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)
