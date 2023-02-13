###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)
library(plyr)
here::i_am("MouseBiogeography-RProj/RegionalGMM-SITE-Heatmap.R")

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")
?generate_adiv_plots()
#Feed in the significant results and generate a target vector with the union of all features 
#setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER")

lumtarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

lum <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                                              lumtarget,
                                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                              Y=metabolic_map,
                                                              "metabolic_map",
                                                              "Luminal",
                                                              cols,
                                                              bk)
dev.new(width=15, height=10)
lum

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1)
lum <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
lum

muc <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     muctarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=metabolic_map,
                                     "metabolic_map",
                                     "Mucosal",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
muc

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-2,-1.5, -1, -0.5, 0, 0.5, 1)
muc <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     muctarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Mucosal",
                                     cols,
                                     bk)
dev.new(width=15, height=10)
muc

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

## Do it the old way 
luminal<-read.table("GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
luminal<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)

target <- lumtarget
luminal_all<-filter(luminal, metadata=="Site")
#length(luminal_all$value[luminal_all$value=="Duodenum"])
data<-luminal_all[luminal_all$feature %in% target, ]


#make an empty dataframe to store the reference variable 
y <- data.frame(matrix(NA,nrow=length(target),ncol=9))
#Assign x, a string vector, to y as its column names:
x <- c(colnames(data))
colnames(y) <- x
y$feature<-target
y$coef <- 0
y$value <- "Distal_Colon"
y$metadata <-"Site"
y$qval<-100

site_heatmap<-rbind(data,y)

site_heatmap$feature <- gsub("X","",as.character(site_heatmap$feature))

#write.csv(site_heatmap,"Mucosal-DCvsall-ggplot-Heatmap.csv")

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", header=TRUE)

data<- (merge(site_heatmap, annotation, by = 'feature'))
data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
data$metabolic_map<-paste(data$Hierarchy_L1,data$annotation,sep=" : ")

#write.csv(data,"Luminal-DCvsall-ggplot-Heatmap.csv")
qval<-data$qval
asterisk<-c("")
for (item in qval){
  if (item < 0.05){
    asterisk<-c(asterisk,"*")
  }
  else {
    asterisk<-c(asterisk,"")
  }
}
asterisk<-asterisk[-1]
data$asterisk<-asterisk
data$value<-factor(data$value, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
ggplotdata<-data
cols=viridis(8)
max(ggplotdata$coef)
min(ggplotdata$coef)
#bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = value, y=hierachy_annotations)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  #geom_text(aes(label=asterisk)) +
  #scale_fill_stepsn(low = muted("aquamarine"), mid = "white", high = muted("violetred"), midpoint = 0) 
  scale_fill_stepsn(n.breaks=10, values = NULL, colors = cols)
g1

dev.new(width=15, height=10)
g1 <- ggplot(ggplotdata, aes(x = value, y=metabolic_map)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  geom_text(aes(label=asterisk)) +
  #scale_fill_stepsn(low = muted("aquamarine"), mid = "white", high = muted("violetred"), midpoint = 0) 
  scale_fill_stepsn(n.breaks=8, values = NULL, colors = cols)
g1
data$hierachy_annotations= factor(as.character(data$hierachy_annotations))
ggplotdata<-data
dev.new(width=15, height=10)
g1 <- ggplot(ggplotdata, aes(x = value, y=hierachy_annotations)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  geom_text(aes(label=asterisk)) +
  #scale_fill_stepsn(low = muted("aquamarine"), mid = "white", high = muted("violetred"), midpoint = 0) 
  scale_fill_stepsn(n.breaks=8, values = NULL, colors = cols)
g1

#construct heatmap using heatmap2 with dendrogram
?pivot_wider
data_long<-pivot_wider(data, id_cols=hierachy_annotations, names_from = value, values_from =coef)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,Duodenum,Jejunum, Ileum,Cecum,Proximal_Colon,Distal_Colon)
row.names(data_long_final)= data_long$hierachy_annotations
matrix.data<- as.matrix.data.frame(data_long_final)
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "Blues"))(5)#, the number (25) represents the number of shades of the gradient
dev.new(width=15, height=10)
heatmap.2(matrix.data, colv= NA, rowv =TRUE, dendrogram ="row", scale="row",density.info="none", trace="none", cexCol = 1, cexRow = 1, margins=c(9,20), col=coul, keysize=1.5, )
#looks like 5 clusters could work

#construct heatmap using heatmap2 with hierarchical clustering
library(dendextend) #Creating a color palette & color breaks

coul <- colorRampPalette(c("aquamarine", "white", "violetred"))(n = 299)
bk = c(seq(-2,-0.25,length=100),  # aquamarine
       seq(-0.24,0.24,length=100), # white
       seq(0.25,2,length=100))    # violetred

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")
?dist
?hclust
cols_branches <- c("red", "green", "orange", "purple", "pink", "grey", "navy") # Set the colors of branches
cols_branches <-c("purple", "pink", "grey", "navy")
cols_branches <-c("cyan","purple", "pink", "grey", "navy")
dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 5, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 5)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
nrow(matrix.data)

data_long_qval<-pivot_wider(data, id_cols=hierachy_annotations, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]


for(i in 1:ncol(data_long_qval)){       # for-loop over columns
  v<-data_long_qval %>% pull(i)
  for(ctr in 1:length(v)){
    if (data_long_qval[ctr,i]<0.05){
      data_long_qval[ctr,i]<- 100
    }
    else {
      data_long_qval[ctr,i]<- 0
    }
  }
}

for(i in 1:ncol(data_long_qval)){ 
  v<-data_long_qval %>% pull(i)
  data_long_qval[,i]<-as.character(v)
}


for(i in 1:ncol(data_long_qval)){       # for-loop over columns
  v<-data_long_qval %>% pull(i)
  for(ctr in 1:length(v)){
    if (data_long_qval[ctr,i]=="100"){
      data_long_qval[ctr,i]<- "*"
    }
    else {
      data_long_qval[ctr,i]<- ""
    }
  }
}

asterisk_matrix<-as.matrix.data.frame(data_long_qval)
dev.new(width=15, height=10)
heatmap.2(matrix.data,
          Colv= FALSE,
          breaks=bk,
          Rowv = as.dendrogram(hcluster),
          symkey=FALSE,
          dendrogram="row",
          scale="none",
          key.xlab="coef",
          density.info="none",
          trace="none",
          cexCol = 1,
          cexRow = 1,
          margins=c(5,25),
          col=coul,
          keysize=0.5,
          RowSideColors = col_labels,
          #cellnote=asterisk_matrix,
          colRow=col_labels)
