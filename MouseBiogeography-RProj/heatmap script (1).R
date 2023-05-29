setwd(choose.dir())
library(gplots)
library(ggplot2)
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "Blues"))(5)#, the number (25) represents the number of shades of the gradient

#coul<- colorRampPalette(c("red", "white", "blue"))(256)

data<-read.csv("MucosalSI-Overlap.csv",header=T, row.names=1)
row.names(data)
names=c(as.character(row.names(data)))
formattednames= gsub("^.*?.o_",".o", names)
formattednames
rownames(data)=formattednames
data <- data[,-c(58:60)]
matrix.data<-as.matrix(data)

heatmap.2(matrix.data, dendrogram = NULL, scale="row",density.info="none", trace="none", cexCol = 1, cexRow = 0.55, margins=c(9,20), col=heat.colors, keysize=1.5, Rowv=FALSE)

?heatmap.2()
  #heatmap.2(matrix.data,density.info="none",col = coul, Colv = NA, Rowv = NA, trace="none", cexCol = 1, cexRow = 0.55, margins=c(5,15))

#wrangling some other genera names 
morenames <- read.csv("Avg_SixFT_only.csv")
names=c(as.character(morenames$Euth_treatment))
formattednames= gsub("^.*?.g",".g", names)
write.csv(formattednames, "Feb14_SixFTonly.g_names.csv")

########Heatmap -Luminal vs Mucosal-Colon
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/Pathway Heatmaps/LuminalvMucosal/")

data<-read.csv("Figure-LuminalvsMucosal-Colon.csv",header=T,row.names=1)
row.names(data)
names=c(as.character(row.names(data)))
formattednames= gsub("^.*?.o_",".o", names)
formattednames
rownames(data)=formattednames
data <- data[,-c(58:60)]
matrix.data<-as.matrix(data)

ggplotdata<-read.csv("ggplot-Colon-LuminalvsMucosal.csv")
g1 <- ggplot(ggplotdata, aes(x = Site, Pathway)) + geom_tile(aes(fill = coef),colour="white",size=0.25) 
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #+ scale_fill_continuous(low = "violetred", high = "aquamarine")
g1

Site <- ggplotdata$Site
Site <- factor(ggplotdata$Site)
Site 
plot1 <- heatmap(matrix.data, Colv = NA, Rowv = NA, scale = "column")
plot1
heatmap.2(matrix.data, colv= NA, rowv =NA, dendrogram ="none", scale="row",density.info="none", trace="none", cexCol = 1, cexRow = 1, margins=c(9,20), col=coul, keysize=1.5, Rowv=FALSE)

?heatmap.2()