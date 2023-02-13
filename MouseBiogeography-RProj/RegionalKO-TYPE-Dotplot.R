library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/")

data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotated_significant_colon.csv", header=TRUE) 
data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotated_Cecum_Significant.csv", header=TRUE) 
data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotated_ProximalColon_Significant.csv", header=TRUE) 
data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotation_DistalColon_Significant.csv", header=TRUE) 
relA <- read.csv("Relative_Abundance-Colon-KO.csv", header=TRUE)


data <- read.csv("clr Maaslin2 Mucosal vs Luminal KO  - Annotated_significant_SI.csv", header=TRUE) 
data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotated_Duodenum_significant.csv", header=TRUE) 
data <- read.csv("Annotated clr Maaslin2 Mucosal vs Luminal KO  - Annotated_Jejunum_significant.csv", header=TRUE) 

relA <- read.csv("Relative_Abundance-SI-KO.csv", header=TRUE)

relA$feature<-relA$X #make sure your search term (ASV sequence) is named the same across whatever you're trying to merge, here it is "feature"
taxonomy <- read.csv("KO_names_Key.csv", header=TRUE)
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
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000002),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme


data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$feature, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$feature= factor(as.character(data$feature), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = feature, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000002),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Colon Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$L2..Superpathway.= factor(as.character(data$L2..Superpathway.), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = L2..Superpathway.)) + 
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

data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme


data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$feature, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$feature= factor(as.character(data$feature), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = feature, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Cecum Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Proximal Colon Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Distal Colon Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme


###Mucosal vs Luminal, Duodenum Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Jejunum Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme

###Mucosal vs Luminal, Ileum Dataset
data<-filter(data, value == "Mucosal" & qval < 0.05)
y = tapply(data$coef, data$annotation, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$annotation= factor(as.character(data$annotation), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = annotation, color=annotation)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.22)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="none") +
  mytheme
