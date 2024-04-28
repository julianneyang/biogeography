###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(ggplot2)
library(tidyr)
library(plyr)
library(dplyr)
library(cowplot)
library(here)
library(dendextend) 

here::i_am("MouseBiogeography-RProj/Final_Figures/L2_HEATMAP_Type_all.R")

remove.packages("Microbiome.Biogeography")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

### Full heatmap colors- modify this according to max and min coef sizes ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)


### CS SPF ---
target <- find_features_union_for_type_heatmap(
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")

#draw heatmap
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_type_L6_heatmap <- generate_taxa_heat_map_by_type_L6(df,
                                                     "CS SPF",
                                                     cols,
                                                     bk)

### HUM Gavage --- NO features
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")


### SPF Gavage --- no features
target <- find_features_union_for_type_heatmap(
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  "Humanized-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

### HUM V Gavage ---
target <- find_features_union_for_type_heatmap(duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
                                               jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
                                               ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
                                               cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
                                               pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
                                               dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")
df <- query_type_features_union(
  target,
  duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


df$Phylum <- gsub(".*\\.p__", "", df$feature)
df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
df$Order <- gsub(".*\\.o__", "", df$feature)
df$Order <- gsub("\\.f__.*", "", df$Order)
df$Order <- paste0(df$Order, " (o)")
df$Family <- gsub(".*\\.f__", "", df$feature)
df$Family <- gsub("\\.g__.*", "", df$Family)
df$Family<- paste0(df$Family, " (f)")
df$Genus <- gsub(".*\\.g__", "", df$feature)
df$Genus <- gsub("\\.g__.*", "", df$Genus)
df$Species <- gsub(".*\\.s__", "", df$feature)

df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))

df <- df %>%
  mutate(level1 = ifelse(nchar(Genus) != 0, Genus, Family))
df <- df %>%
  mutate(annotation = ifelse(level1!= " (f)", level1, Order))

data <-df 
data$Phylum_Genus<-paste(data$Phylum,data$annotation,sep=" : ")

qval<-data$qval
print(qval)
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

# Make sure max and min coef are 2 and -2 so the scale stays the same
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 

# orders the genera by the mean fold change within row 
y = tapply(data$coef_d, data$Genus, function(y) mean(y))  
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
data$Site = factor(data$Site, levels=c("D", "J", "I", "C", "PC", "DC"))

#Heatmap with clustering
data <- unique(data)
data <- data %>% select(-c("value", "metadata"))
data <- unique(data)
data$coef_d <- as.numeric(data$coef_d)
data_long<-pivot_wider(data, id_cols=Phylum_Genus, names_from = Site, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,D,J, I,C,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)

#construct heatmap using heatmap2 with hierarchical clustering


coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple", "navy","black") # Set the colors of branches
#cols_branches <-c("purple", "pink", "grey", "navy")
#cols_branches <-c("cyan","purple", "pink", "grey", "navy")
dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from =Site, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from = Site, values_from =qval)
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
          labCol = c("D", "J","I","C","PC","DC"),
          colRow=col_labels,
          srtCol=0)

### UCLA O SPF ---
target <- find_features_union_for_type_heatmap(
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv")

df <- query_type_features_union(
  target,
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv")

df$Phylum <- gsub(".*\\.p__", "", df$feature)
df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
df$Order <- gsub(".*\\.o__", "", df$feature)
df$Order <- gsub("\\.f__.*", "", df$Order)
df$Order <- paste0(df$Order, " (o)")
df$Family <- gsub(".*\\.f__", "", df$feature)
df$Family <- gsub("\\.g__.*", "", df$Family)
df$Family<- paste0(df$Family, " (f)")
df$Genus <- gsub(".*\\.g__", "", df$feature)
df$Genus <- gsub("\\.g__.*", "", df$Genus)
df$Species <- gsub(".*\\.s__", "", df$feature)

df$Family_Species <- paste(df$Family,  gsub("^.*_","",df$Species))
df$Order_Species <- paste(df$Order,  gsub("^.*_","",df$Species))

df <- df %>%
  mutate(level1 = ifelse(nchar(Genus) != 0, Genus, Family))
df <- df %>%
  mutate(annotation = ifelse(level1!= " (f)", level1, Order))

data <-df 
data$Phylum_Genus<-paste(data$Phylum,data$annotation,sep=" : ")

qval<-data$qval
print(qval)
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

# Make sure max and min coef are 2 and -2 so the scale stays the same
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 

# orders the genera by the mean fold change within row 
y = tapply(data$coef_d, data$Genus, function(y) mean(y))  
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
data$Site = factor(data$Site, levels=c("D", "J", "I", "C", "PC", "DC"))

#Heatmap with clustering
data <- unique(data)
data <- data %>% select(-c("value", "metadata"))
data <- unique(data)
data$coef_d <- as.numeric(data$coef_d)
data_long<-pivot_wider(data, id_cols=Phylum_Genus, names_from = Site, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,D,J, I,C,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)

#construct heatmap using heatmap2 with hierarchical clustering


coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple", "navy","black") # Set the colors of branches
#cols_branches <-c("purple", "pink", "grey", "navy")
#cols_branches <-c("cyan","purple", "pink", "grey", "navy")
dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from =Site, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from = Site, values_from =qval)
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
          labCol = c("D", "J","I","C","PC","DC"),
          colRow=col_labels,
          srtCol=0)

### Plot the figure ---

heatmap_hum <- plot_grid(NULL,hum_type_L2_heatmap,NULL, nrow=3)                                                      

dev.new()
plot_grid(ucla_o_type_L2_heatmap, cs_type_L2_heatmap, 
          hum_v_gavage_type_L2_heatmap, NULL,
          nrow=2, ncol=2,
          labels=c("A","B","C", "D"))

