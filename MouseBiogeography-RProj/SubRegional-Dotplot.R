library(ggplot2)
library(ggpubr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2/")

data<-read.csv("Maaslin2 for All Sites - LuminalColon-Maaslin2.csv",header=T)  
data<-read.csv("Maaslin2 for All Sites - LuminalSI-Maaslin2.csv",header=T)  
data<-read.csv("Maaslin2 for All Sites - MucosalColon-Maaslin2.csv",header=T)  
data<-read.csv("Maaslin2 for All Sites - MucosalSI-Maaslin2.csv",header=T)  

data<-as.data.frame(data)  
#data<-data[data$Line=="SLCA393T"|data$Line=="Vil1Cre_Pos",]

min(data$Relative_Abundance)
max(data$Relative_Abundance)
#data$Line<- factor(data$Line,levels=c("Vil1Cre_Pos","SLCA393T"))


y = tapply(data$coef, data$Taxon, function(y) max(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Taxon= factor(as.character(data$Taxon), levels = names(y))
data$Relative_Abundance
dev.new(width=15, height=15)  # can adjust window size of the plot output this way
p<- ggplot(data, aes(x = coef, y = Taxon, color=value)) + geom_point(aes(size = sqrt(Relative_Abundance))) + scale_size_continuous(name="Abundance",range = c(0.5,8),limits=c(sqrt(0.0000295),sqrt(0.09)),breaks=c(sqrt(0.0001),sqrt(0.001),sqrt(0.01),sqrt(0.1)),labels=c("0.0001","0.001","0.01","0.1")) + geom_vline(xintercept = 0) + theme_pubr()+xlab(label="Log2 Fold Change")+ylab(label=NULL)+theme(legend.position="right",legend.direction="vertical")
p + ggtitle("Site Differences in Small Intestinal Mucosal Samples")
p
