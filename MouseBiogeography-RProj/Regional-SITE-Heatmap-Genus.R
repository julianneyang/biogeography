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
library(gplots)
library(dendextend) 

here::i_am("MouseBiogeography-RProj/Regional-Site-Heatmap-Genus.R")

## Function to clean up the script --
generate_matrix_for_heatmap_clustering <- function(path_to_significant_results,
                                                   path_to_annotation_file,
                                                   path_to_all_results){
  
## First, find all the features that are significant in at least one comparison --
luminal<-readr::read_delim(here({{path_to_significant_results}}), delim="\t")

luminal <- as.data.frame(luminal)

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
print(unique(target))

## Second, Match target taxa to cleaned names in annotation file ---


#Lum
annotation <- readr::read_delim(here(path_to_annotation_file))

regionalgenera<-target
#regionalgenera <- gsub("X", "", regionalgenera)
regionalgenera<-as.data.frame(regionalgenera)
regionalgenera$feature <- regionalgenera[,1]
tempdf<-merge(regionalgenera,annotation, by= "feature")
regionalgenera<-as.data.frame(tempdf$annotation)
regionalgenera <- regionalgenera$`tempdf$annotation`
print(unique(regionalgenera))
here()
#readr::write_rds(regionalgenera, file=here(paste0(filepath,"regionalluminalgenera.RDS")))

## Query the target vector against all_results.tsv ---
luminal<-readr::read_delim(here({{path_to_all_results}}), delim="\t")
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

#site_heatmap$feature <- gsub("X","",as.character(site_heatmap$feature))
#write.csv(site_heatmap,"SITE Genus Heatmap.csv")

## Third, wrangle the dataframe into something that can be used for heatma
library(viridis)
annotation <- readr::read_delim(here({{path_to_annotation_file}}))

annotation$Phylum <- gsub(".*\\.p__", "", annotation$feature)
annotation$Phylum <- gsub("\\.c__.*", "", annotation$Phylum)
data<- (merge(site_heatmap, annotation, by = 'feature'))
data <- unique(data)
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
data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
data$value = factor(data$value, levels=c("D", "J", "I", "C", "PC", "DC"))
ggplot_data <- unique(data)
ggplot_data$Phylum_Genus<-paste(ggplot_data$Phylum,ggplot_data$annotation,sep=" : ")


#construct heatmap using heatmap2 with dendrogram
data_long<-pivot_wider(ggplot_data, id_cols=Phylum_Genus, names_from = value, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,D,J, I,C,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)


#add asterisks for qval <0.05

data_long_qval<-pivot_wider(ggplot_data, id_cols=annotation, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval <- select(data_long_qval,c("D","I", "J","C","PC","DC"))

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

newList <- list("matrix" = matrix.data, "dataframe" = ggplot_data,
                "asterisks"=asterisk_matrix)
return(newList)
}


## UCLA O. SPF Luminal --
filepath <- "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/"
path_to_significant_results <- paste0(filepath,"L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")
path_to_annotation_file <- paste0(filepath,"Genus_Luminal_taxonomy.csv")
path_to_all_results <- paste0(filepath,"L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv")

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                      path_to_annotation_file = path_to_annotation_file,
                                                      path_to_all_results = path_to_all_results)
matrix.data <- results$matrix
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("purple", "royalblue", "forestgreen", "firebrick") # Set the colors of branches


dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2

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
          cellnote=results$asterisks,
          notecol="black",
          labCol = c("D", "J","I","C","PC","DC"),
          colRow=col_labels,
          srtCol=0)


## UCLA O. SPF Mucosal --
filepath <- "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/"
path_to_significant_results <- paste0(filepath,"L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")
path_to_annotation_file <- paste0(filepath,"Genus_Mucosal_taxonomy.csv")
path_to_all_results <- paste0(filepath,"L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv")

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
matrix.data <- results$matrix
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple") # Set the colors of branches


dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2

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
          cellnote=results$asterisks,
          notecol="black",
          labCol = c("D", "J","I","C","PC","DC"),
          colRow=col_labels,
          srtCol=0)

