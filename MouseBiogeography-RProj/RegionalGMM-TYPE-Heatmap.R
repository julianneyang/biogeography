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

here::i_am("MouseBiogeography-RProj/RegionalGMM-TYPE-Heatmap.R")

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/")
duodenum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  a<-duodenum_significant$feature
jejunum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  b<-jejunum_significant$feature
ileum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  ileum_significant<-filter(ileum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  c<-ileum_significant$feature
cecum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
  d<-cecum_significant$feature  
pc<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
  e<-pc_significant$feature  
DC<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
  DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
  f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

#Query the target vector against all_results.tsv for each site
duodenum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
  duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
  duodenum_all$Site<- "Duodenum"
jejunum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
  jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
  jejunum_all$Site<- "Jejunum"
ileum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
  ileum_all<-ileum_all[match(target,ileum_all$feature),]
  ileum_all$Site<- "Ileum"
cecum<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
  cecum_all<-cecum_all[match(target,cecum_all$feature),]
  cecum_all$Site<- "Cecum"
pc<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
  pc_all<-pc_all[match(target,pc_all$feature),]
  pc_all$Site<- "Proximal_Colon"
DC<-read.table("GMM-Maaslin2-TYPE/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
  DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
  DC_all<-DC_all[match(target,DC_all$feature),]
  DC_all$Site<- "Distal_Colon"

duojej<-rbind(duodenum_all,jejunum_all)
ilecec<-rbind(ileum_all, cecum_all)
pcdc<-rbind(pc_all,DC_all)
duojejilecec<-rbind(duojej,ilecec)
duojejilececpcdc<-rbind(duojejilecec,pcdc)

#write.csv(duojejilececpcdc, "GMM-TYPE-Heatmap.csv") #from here make sure all NA rows are filled with feature name corresponding to NA

#remove across six sites all GMM that failed to converge
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
duojejilececpcdc<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

# Use functions to construct the heatmap
library(viridis)
cols<-viridis(8)
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)



cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF")

bk =c(-1, -0.5, 0, 0.5)
mucvlum <- generate_GMM_heat_map_by_Type(gmm_heatmap_final, 
                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                         Map, 
                                         "Map", 
                                         "Mucosal vs Luminal",
                                         cols,
                                         bk)
dev.new(width=15, height=10)
mucvlum

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
mucvlum <- generate_GMM_heat_map_by_Type(gmm_heatmap_final, 
                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                         metabolic_map, 
                                         "metabolic_map", 
                                         "Mucosal vs Luminal",
                                         cols,
                                         bk)
dev.new(width=15, height=10)
mucvlum

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("Revised_Module_Key.csv", header=TRUE)
data<- (merge(gmm_heatmap_final, annotation, by = 'feature'))

data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")

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
data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
ggplotdata<-data
cols=viridis(5)
max(ggplotdata$coef)
min(ggplotdata$coef)
#bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = Site, y=hierachy_annotations)) + geom_tile(aes(fill = coef),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  #geom_text(aes(label=asterisk)) +
  #scale_fill_stepsn(low = muted("aquamarine"), mid = "white", high = muted("violetred"), midpoint = 0) 
  scale_fill_stepsn(n.breaks=5, values = NULL, colors = cols)
g1

#construct heatmap using heatmap2 with dendrogram
?pivot_wider
data_long<-pivot_wider(data, id_cols=hierachy_annotations, names_from = Site, values_from =coef)
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

#coul <- colorRampPalette(c("aquamarine", "white", "violetred"))(n = 299)
coul<-viridis(8)
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="mcquitty")
?hclust
?dist

cols_branches <- rainbow(6)
dend1 <- color_branches(dend1, k = 6, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 6)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
nrow(matrix.data)

data_long_qval<-pivot_wider(data, id_cols=hierachy_annotations, names_from = Site, values_from =qval)
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
