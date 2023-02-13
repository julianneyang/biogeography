library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/[Archive] clr Maaslin2 Pathway/")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction")

data <- read.csv("Pathways_Dotplot/clr Maaslin2 Pathway Luminal Colon vs SI - Colonref_Heatmap.csv", header=TRUE) 
relA <- read.csv("Pathways_Dotplot/Relative_Abundance-Luminal-Pathways.csv", header=TRUE)
relA$X <- gsub("-",".", relA$X)

data <- read.csv("Pathways_Dotplot/[TYPE] clr Maaslin2 Pathways Mucosal vs Luminal - Colon_MucosalvsLuminal.csv", header=TRUE) 
relA <- read.csv("Pathways_Dotplot/Relative_Abundance-Colon-Pathways.csv", header=TRUE)
relA$X <- gsub("-",".", relA$X)

data <- read.csv("Pathways_Dotplot/[TYPE] clr Maaslin2 Pathways Mucosal vs Luminal - SI_MucosalvsLuminal.csv", header=TRUE) 
relA <- read.csv("Pathways_Dotplot/Relative_Abundance-SI-Pathways.csv", header=TRUE)
relA$X <- gsub("-",".", relA$X)


relA$feature<-relA$X #make sure your search term (ASV sequence) is named the same across whatever you're trying to merge, here it is "feature"
taxonomy <- read.csv("Pathways_Dotplot/MetaCyc_pathwaynames_Key.csv", header=TRUE)
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

###SI vs Colon, Luminal Dataset
data<-filter(data, value=="SI", qval < 0.05)
y = tapply(data$coef, data$L2A, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$L2A= factor(as.character(data$L2A), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = L2A, color = L4A)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
?tapply()

#heatmap in ggplot
y = tapply(data$coef, data$L3A, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$L3A= factor(as.character(data$L3A), levels = names(y))
ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = value, y=L3A)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
# As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "aquamarine", high = "violetred")
g1

###Mucosa vs Lumen, Colon Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$L1...Terminal.Pathway, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$L1...Terminal.Pathway= factor(as.character(data$L1...Terminal.Pathway), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = L1...Terminal.Pathway, color = L4A)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
?tapply()

ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = value, y=L1...Terminal.Pathway)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "aquamarine", high = "violetred")
g1

###Mucosa vs Lumen, SI Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$L1...Terminal.Pathway, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$L1...Terminal.Pathway= factor(as.character(data$L1...Terminal.Pathway), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = L1...Terminal.Pathway, color = L4A)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
?tapply()

ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = value, y=L1...Terminal.Pathway)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_continuous(low = "aquamarine", high = "violetred")
g1
