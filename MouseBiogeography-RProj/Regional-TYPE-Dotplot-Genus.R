library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/")

data <- read.csv("clr Maaslin2 TYPE Genus Level - Colon.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Cecum.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Proximal_Colon.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Distal_Colon.csv", header=TRUE) 
relA <- read.csv("Relative_Abundance-Colon-L6.csv", header=TRUE)


data <- read.csv("clr Maaslin2 TYPE Genus Level - SI.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Duodenum.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Jejunum.csv", header=TRUE) 
data <- read.csv("clr Maaslin2 TYPE Genus Level - Ileum.csv", header=TRUE) 
relA <- read.csv("Relative_Abundance-SI-L6.csv", header=TRUE)

relA$feature<-relA$X #make sure your search term (ASV sequence) is named the same across whatever you're trying to merge, here it is "feature"
taxonomy <- read.csv("genus_taxonomy.csv", header=TRUE)
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

###Mucosal vs Luminal, SI Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.53)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("SI: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal vs Luminal, Colon Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.34)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Colon: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
?tapply()

###Mucosal vs Luminal, Cecum Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.34)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Cecum: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal vs Luminal, Proximal Colon Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.34)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Proximal Colon: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal vs Luminal, Distal Colon Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.34)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Distal_Colon: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)


###Mucosal vs Luminal, Duodenum Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.34)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Duodenum: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal vs Luminal, Jejunum Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.53)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Jejunum: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal vs Luminal, Ileum Dataset
data<-filter(data, metadata=="Type" & value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color=Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.53)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Ileum: Mucosal vs. Luminal") + 
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
