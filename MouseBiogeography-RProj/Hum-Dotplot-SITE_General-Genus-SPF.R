library(ggplot2)
library(ggpubr)
library(dplyr)
library(funrar)
library(dplyr)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

data<-read.table("Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
data<- data%>% filter(metadata=="Site_General" & qval<0.05)
luminal<-read.table("Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/significant_results.tsv", header=TRUE)
data<- luminal%>% filter(metadata=="Site_General" & qval<0.05)

#Calculate Relative Abundance 
df_input_data <- read.csv("Source RPCA/SPF/Maaslin2_Site_L6/Luminal_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
df_input_data <- select(df_input_data, -c("taxonomy"))
transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
relA <- as.data.frame(t(Relative_Abundance))

relA$feature<-row.names(relA) #make sure your search term (ASV sequence) is named the same across whatever you're trying to merge, here it is "feature"
relA$feature<-gsub(";",".",relA$feature)
relA$feature<-gsub("-",".",relA$feature)
relA$feature<-gsub("/",".",relA$feature)

data<- (merge(data, relA, by = 'feature'))
data$Genus<-gsub(".*g__","",data$feature)
data$Phylum<-gsub(".*p__","",data$feature)
data$Phylum <- gsub("\\.c.*","",data$Phylum)
#write.csv(data,"Source RPCA/SPF/Maaslin2_Site_L6/Significant_Luminal_Site_General.csv")

data<-read.csv("Source RPCA/SPF/Maaslin2_Site_L6/Significant_Luminal_Site_General.csv",header=TRUE)
data$Relative_Abundance<-data$V1


###Luminal Site General
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme_cowplot(16)+
  theme(legend.position="top") +
  ggtitle("Luminal: SI vs colon(ref)")+
  guides(fill=guide_legend(nrow=2, byrow=TRUE))
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#Ileum as reference
data<-filter(data, metadata == "Site" & value=="Jejunum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal SI: Jejunum vs Ileum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

data<-filter(data, metadata == "Site" & value=="Duodenum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal SI: Duodenum vs Ileum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)


###Mucosal SI
#Duodenum as reference
data<-filter(data, metadata == "Site" & value=="Jejunum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal SI: Jejunum vs Duodenum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#Ileum as reference
data<-filter(data, metadata == "Site" & value=="Jejunum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal SI: Jejunum vs Ileum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

data<-filter(data, metadata == "Site" & value=="Duodenum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal SI: Duodenum vs Ileum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Luminal Colon
#Cecum as reference
data<-filter(data, metadata == "Site" & value=="Proximal_Colon" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal Colon: Proximal Colon vs Cecum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#DC as reference
data<-filter(data, metadata == "Site" & value=="Proximal_Colon" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal Colon: Proximal Colon vs Distal Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

data<-filter(data, metadata == "Site" & value=="Cecum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal Colon: Cecum vs Distal Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal Colon
#DC as reference
data<-filter(data, metadata == "Site" & value=="Proximal_Colon" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal Colon: Proximal Colon vs Distal Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

data<-filter(data, metadata == "Site" & value=="Cecum" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal Colon: Cecum vs Distal Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)
#Cecum as reference
data<-filter(data, metadata == "Site" & value=="Proximal_Colon" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal Colon: Proximal Colon vs Cecum (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)


###Luminal 
#Colon as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal: SI vs Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#SI as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Luminal: Colon vs SI (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

###Mucosal
#Colon as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal: SI vs Colon (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)

#SI as reference
data<-filter(data, metadata == "Site_General" & qval < 0.05)
y = tapply(data$coef, data$Genus, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data, aes(x = coef, y = Genus, color = Phylum)) + 
  geom_point(aes(size = sqrt(Relative_Abundance))) + 
  scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.000004),sqrt(0.5)),breaks=c(sqrt(0.00001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.00001","0.001","0.01","0.1")) + 
  geom_vline(xintercept = 0) + 
  xlab(label="Log2 Fold Change")+
  ylab(label=NULL)+
  theme(legend.position="top") +
  ggtitle("Mucosal: Colon vs SI (ref)")+
  mytheme
max(data$Relative_Abundance)
min(data$Relative_Abundance)