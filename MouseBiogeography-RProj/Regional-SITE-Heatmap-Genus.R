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
library(here)
library(dplyr)
library(plyr)

here::i_am("MouseBiogeography-RProj/Regional-Site-Heatmap-Genus.R")

###for SITE:DC vs all data
#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
luminal<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv", header=TRUE)
luminal<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv", header=TRUE)

duodenum_significant<-filter(luminal, metadata=="Site" & value=="Duodenum" &qval<0.05)
a<-duodenum_significant$feature
jejunum_significant<-filter(luminal, metadata=="Site" & value=="Jejunum" &qval<0.05)
b<-jejunum_significant$feature
ileum_significant<-filter(luminal, metadata=="Site" & value=="Ileum" &qval<0.05)
c<-ileum_significant$feature
cecum_significant<-filter(luminal, metadata=="Site" & value=="Cecum" &qval<0.05)
d<-cecum_significant$feature  
pc_significant<-filter(luminal, metadata=="Site" & value=="Proximal_Colon" &qval<0.05)
e<-pc_significant$feature  
DC_significant<-filter(luminal, metadata=="Site" & value=="Distal_Colon" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

## Grab the Genera names for the venn diagram ---
regionalgenera<-target
regionalgenera <- gsub("X", "", regionalgenera)
regionalgenera<-as.data.frame(regionalgenera)
regionalgenera$feature <- regionalgenera[,1]

#Lum
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
#Muc
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Mucosal_taxonomy.csv", header=TRUE)

annotation$feature<-annotation$X
tempdf<-merge(regionalgenera,annotation, by= "feature")
regionalgenera<-as.data.frame(tempdf$Genus)
regionalgenera <- regionalgenera$`tempdf$Genus`

here()
readr::write_rds(regionalgenera, file=here("regionalluminalgenera.RDS"))
here()
readr::write_rds(regionalgenera, file=here("regionalmucosalgenera.RDS"))

## Query the target vector against all_results.tsv ---
luminal<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
luminal<-read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/L6-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
#luminal<-read.table("L6-DuodvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
#luminal<-read.table("L6-DuodvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)

luminal_all<-filter(luminal, metadata=="Site")
#length(luminal_all$value[luminal_all$value=="Duodenum"])
data<-luminal_all[luminal_all$feature %in% target, ]

length(target)
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
#write.csv(site_heatmap,"SITE Genus Heatmap.csv")

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Luminal_taxonomy.csv", header=TRUE)
annotation <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/genus_Mucosal_taxonomy.csv", header=TRUE)

annotation$feature<-annotation$X
data<- (merge(site_heatmap, annotation, by = 'feature'))
data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")

qval<-data$qval
asterisk<-c("")
for (item in qval){
  if (item < 0.05){
    asterisk<-c(asterisk,"*")
  }
  else if (item=="NA"){}
  else {
    asterisk<-c(asterisk,"")
  }
}
asterisk<-asterisk[-1]
data$asterisk<-asterisk


data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 
y = tapply(data$coef_d, data$Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
ggplotdata<-data


#construct heatmap using heatmap2 with dendrogram
?pivot_wider
data_long<-pivot_wider(data, id_cols=Phylum_Genus, names_from = value, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,Duo,Jej, Ile,Cec,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)

#construct heatmap using heatmap2 with hierarchical clustering
library(dendextend) #Creating a color palette & color breaks

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

# Luminal Cluster Colors
cols_branches <- c("purple", "forestgreen", "royalblue", "firebrick") # Set the colors of branches
# Mucosal Cluster Colors
cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple") # Set the colors of branches


dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
nrow(matrix.data)

data_long_qval<-pivot_wider(data, id_cols=feature, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]
data_long_qval <- select(data_long_qval,c("Duo","Ile", "Jej","Cec","PC","DC"))

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
row.names(matrix.data) <- c()
dev.new(width=15, height=10)
heatmap.2(matrix.data,
          Colv= FALSE,
          breaks=bk,
          Rowv = as.dendrogram(hcluster),
          symkey=FALSE,
          dendrogram="none",
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
          cellnote=asterisk_matrix,
          notecol="black",
          labCol = c("Duo", "Jej","Ile","Cec","PC","DC"),
          colRow=col_labels,
          srtCol=0)
?heatmap.2
