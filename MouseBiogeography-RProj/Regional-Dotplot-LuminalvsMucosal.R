library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis")

data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_Cecum.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_DistalColon.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_ProximalColon.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_Colon.csv", header=TRUE) 
relA <- read.csv("Maaslin2_LuminalvsMucosal/Relative_Abundance-Colon-ComBat.csv", header=TRUE)

data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_Duodenum.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_Jejunum.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_Ileum.csv", header=TRUE) 
data <- read.csv("Maaslin2_LuminalvsMucosal/clr Maaslin2 Luminal vs Mucosal - LumRef_SI.csv", header=TRUE) 
relA <- read.csv("Maaslin2_LuminalvsMucosal/Relative_Abundance-SI-ComBat.csv", header=TRUE)

relA$feature<-relA$X #make sure your search term (ASV sequence) is named the same across whatever you're trying to merge, here it is "feature"
taxonomy <- read.csv("Maaslin2_LuminalvsMucosal/taxonomy.csv", header=TRUE)
taxonomy$feature<-taxonomy$X.OTU.ID
intermediate<- (merge(data, relA, by = 'feature'))
data<- (merge(intermediate, taxonomy, by = 'feature'))
data$Relative_Abundance<-data$V1

mytheme <- theme(panel.grid.minor=element_blank(), #gets rid of grey and lines in the middle
                 panel.grid.major=element_blank(), #gets rid of grey and lines in the middle
                 panel.background=element_rect(fill="white"),#gets rid of grey and lines in the middle
                 panel.border=element_blank(), #gets rid of square going around the entire graph
                 axis.line.x = element_line(colour = 'black', size = 0.5),#sets the axis line size
                 axis.line.y = element_line(colour = 'black', size = 0.5),#sets the axis line size
                 axis.ticks=element_line(colour = 'black', size = 0.5), #sets the tick lines
                 axis.title.x = element_text(size=16, color="black"), #size of x-axis title
                 axis.title.y = element_text(size=16, color="black"), #size of y-axis title
                 axis.text.x = element_text(size=12, color="black"), #size of x-axis text
                 axis.text.y = element_text(size=12, color="black"))#size of y-axis text

###Duodenum
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Jejunum
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)


###Ileum
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.25)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)


###SI
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.25)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Colon
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Cecum
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Proximal Colon
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Distal Colon
data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme

data<-filter(data, metadata == "Type" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Luminal Colon
#Cecum as reference
data<-filter(data, metadata == "Site" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = value)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#DC as reference
data<-filter(data, metadata == "Site" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme

###Mucosal Colon
#DC as reference
data<-filter(data, metadata == "Site" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#Cecum as reference
data<-filter(data, metadata == "Site" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = value)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme

###Luminal 
#Colon as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=20)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#SI as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=20)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme

###Mucosal
#Colon as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=20)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + #switch to Value for Site
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#SI as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
dev.new(width=15, height=20)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Taxon, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme