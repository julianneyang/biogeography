library(ggplot2)
library(ggpubr)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Type_Modules_Results/Heatmap/")

data <- read.csv("OMIXER TYPE Module Test Results - Combined_Colon_Significant_OMIXER_MucosalvsLuminal.csv", header=TRUE) 
data <- read.csv("OMIXER TYPE Module Test Results - Combined_SI_significant_OMIXER_MucosalvsLuminal.csv", header=TRUE) 
data <- read.csv("OMIXER TYPE Module Test Results - Combined_SIandColon_OMIXER_MucosalvsLuminal.csv", header=TRUE) 

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


#Combined Colon
data$Site<-factor(data$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
y = tapply(data$Rank_color, data$Feature, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Feature= factor(as.character(data$Feature), levels = names(y))
ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = Site, y=Feature)) + geom_tile(aes(fill = Rank_color),colour="white",size=0.25) +
# As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  theme(axis.text.x = element_text(angle = 90))
g1

#Combined SI
data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
y = tapply(data$Rank_color, data$Feature, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Feature= factor(as.character(data$Feature), levels = names(y))
ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = Site, y=Feature)) + geom_tile(aes(fill = Rank_color),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  theme(axis.text.x = element_text(angle = 90))

#Combined all
data$Site <-data$Taxon
data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
y = tapply(data$Rank_color, data$Feature, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Feature= factor(as.character(data$Feature), levels = names(y))
ggplotdata<-data
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = Site, y=Feature)) + geom_tile(aes(fill = Rank_color),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") 
g1
g1