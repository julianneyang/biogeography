library(data.table)
library(stringi)
library(stringr)
library(funrar)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)

here::i_am("MouseBiogeography-RProj/Donors-Type-Heatmap.R")


target <- find_features_union_for_type_heatmap(duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv for each site
df <- query_type_features_union(
  target,
  duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv")


#from here make sure all NA rows are filled with feature name corresponding to NA
# from here make sure all NA rows are filled with feature name corresponding to NA via copy paste
write.csv(df, here("Donors-Analysis/type_subsets/L6_type_heatmap.csv"))

heatmap<-df
discard<- heatmap[is.na(heatmap$metadata), ]
offtarget<- discard$feature
offtarget<-unique(offtarget)
heatmap_final<-subset(heatmap,  !heatmap[,3] %in% offtarget )

write.csv(offtarget, here("Donors-Analysis/type_subsets/L6_omitted_taxa.csv"))

## Draw heatmap of non-omitted taxa --
library(viridis)
head(df)

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

ggplotdata<-data


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
g1 <- ggplot(ggplotdata, aes(x = Site, y=Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  # As our values are continuous, we'll use scale_fill_continuous instead of scale_fill_manual
  #scale_fill_gradient2(low = "aquamarine", mid = "white", high = "violetred") + 
  geom_text(aes(label=asterisk)) +
  #scale_fill_stepsn(low = muted("aquamarine"), mid = "white", high = muted("violetred"), midpoint = 0) 
  scale_fill_stepsn(n.breaks=9, values = NULL, colors = cols) 
g1

### Heatmap2
data <- unique(data)
data <- data %>% select(-c("X","value", "metadata"))
data <- unique(data)
data$coef_d <- as.numeric(data$coef_d)
data_long<-pivot_wider(data, id_cols=Phylum_Genus, names_from = Site, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,Duo,Jej, Ile,Cec,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "Blues"))(5)#, the number (25) represents the number of shades of the gradient
dev.new(width=15, height=10)
heatmap.2(matrix.data, colv= NA, rowv =TRUE, dendrogram ="row", scale="row",density.info="none", trace="none", cexCol = 1, cexRow = 1, margins=c(9,20), col=coul, keysize=1.5, )
#looks like 6 clusters could work

#construct heatmap using heatmap2 with hierarchical clustering
library(dendextend) #Creating a color palette & color breaks

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

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
nrow(matrix.data)

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
obj + ggtitle("UCLA O. SPF Transverse")
