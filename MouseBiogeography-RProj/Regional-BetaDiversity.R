library(GGally)
library(ggplot2)
library(car)
library(vegan)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis")
data <- read.csv("Bray Curtis PCoA/LumCol.csv", row.names=1)
data <- read.csv("Bray Curtis PCoA/MucCol.csv", row.names=1)
data <- read.csv("Bray Curtis PCoA/LumSI.csv", row.names=1)
data <- read.csv("Bray Curtis PCoA/MucSI.csv", row.names=1)
data <- read.csv("Bray Curtis PCoA/Mucosal.csv", row.names=1)
data <- read.csv("Bray Curtis PCoA/Luminal.csv", row.names=1)

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

sapply(data,levels)
data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))

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
p + labs(title="Luminal Site Differences Bray-Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/Luminal-dist/distance-matrix.tsv")
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
p + labs(title="Mucosal Site Differences Bray-Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/Mucosal-dist/distance-matrix.tsv")
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
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site.1)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  theme(axis.text = element_text(colour = "black", size="20"))
p + labs(title="LuminalColon Site Differences Bray Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/LumCol-dist/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 16S 
#https://sites.ualberta.ca/~ahamann/teaching/graphics/LabRDA.pdf
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/LumCol-dist/distance-matrix.tsv", row.names=1)
bray_dist_16s <- as.dist(as(data.dist, "matrix"))
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
plot(cap_16s, type = "n")
col.list=cols
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19)
ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
plot(ef_16s)


###Mucosal Colon
cols <- c("Cecum" = "cyan", "Proximal_Colon" = "blue", "Distal_Colon" = "magenta")
p <- ggplot(data, aes(x=PC1, y=PC2, colour=Site.1)) + 
  mytheme + geom_point(size=3) + xlab("PC1") + ylab("PC2") +
  #scale_shape_manual(name="Genotype", values=c(16,17,10)) +
  #stat_ellipse()+
  #geom_label(label=MouseID) +
  #scale_fill_manual(name="Genotype", values=cols) +
  scale_colour_manual(name="Site", values=cols) +
  #geom_text(aes(label=MouseID, hjust=0.40, vjust=1.50)) +
  theme(legend.text = element_text(size = 09), legend.title = element_text(size= 9)) + 
  theme(axis.text = element_text(colour = "black", size="20"))
p + labs(title="MucosalColon Site Differences Bray Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data<-read.table(file ="Bray-Curtis-Distance-Matrices/MucCol-dist/distance-matrix.tsv")
data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

# dbRDA 16S Mucosal colon
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/MucCol-dist/distance-matrix.tsv", row.names=1)
bray_dist_16s <- as.dist(as(data.dist, "matrix"))
cap_16s <- capscale(bray_dist_16s ~ Sequencing_Run + Line + Sex + Site.1,data)
plot(cap_16s, type = "n")
col.list=cols
points(cap_16s, display = "sites", col = col.list[paste(data$Site.1)],pch=19)
ef_16s <- envfit(cap_16s, data[, c("Sequencing_Run", "Line", "Sex", "Site.1")], permu = 10000)  # only show significant covariates
ef_16s
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
p + labs(title="MucosalSI Site Differences Bray-Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/MucSI-dist/distance-matrix.tsv")
data.dist <- as.dist(as(data.dist, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

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
p + labs(title="LuminalSI Site Differences Bray-Curtis")

#Adonis
metadata<-data
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site.1)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)
data.dist<-read.table(file ="Bray-Curtis-Distance-Matrices/LumSI-dist/distance-matrix.tsv")
data.dist <- as.dist(as(data, "matrix"))
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
row.names(data) == row.names(metadata)

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



##################
#Adonis with the Bray-Curtis Matrix, using overlapping ASVs corrected file 

data<-read.table(file ="", row.names=1)
data.dist <- as.dist(as(data, "matrix"))

metadata<-read.table(file="Mucosal-Overlap-Metadata.tsv",header=T,row.names=1) #mapping file
metadata<-read.csv(file="Luminal Regional ASV and Metadata - Luminal-Regional-metadata.csv",header=T,row.names=1) #mapping file
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis

#Adonis with the Bray-Curtis Matrix, using overlapping ASVs and CombatSeq corrected file 

data<-read.table(file ="ComBat-Mucosal-distance-matrix.tsv", row.names=1)
data.dist <- as.dist(as(data, "matrix"))

metadata<-read.table(file="Mucosal-Overlap-Metadata.tsv",header=T,row.names=1) #mapping file
metadata<-read.csv(file="Luminal Regional ASV and Metadata - Luminal-Regional-metadata.csv",header=T,row.names=1) #mapping file
metadata$Sex <- factor(metadata$Sex)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Site <- factor(metadata$Site)
metadata$Line <- factor(metadata$Line)
metadata$MouseID_Line <- factor(metadata$MouseID_Line)

sapply(metadata,levels)
data.adonis=adonis(data.dist ~ Sequencing_Run + Line + Sex + Site, data=metadata, permutations=10000)
data.adonis
