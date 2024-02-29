library(gplots)
library(dendextend) 
library(here)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_4_Genus_Site_Heatmaps_Clustering.R")


### Heatmap ---

## 4A. UCLA O. SPF Luminal --
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


## CS SPF Luminal --
path_to_significant_results <- "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv"
path_to_annotation_file <- "CS-Facility-Analysis/differential_genera_site/Genus_Luminal_taxonomy.csv"
path_to_all_results <- "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"

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

## HUM V Gavage Luminal --
path_to_significant_results <- "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"
path_to_annotation_file <- "Donors-Analysis/differential_genera_site/Genus_Luminal_taxonomy.csv"
path_to_all_results <- "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
matrix.data <- results$matrix
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple") # Set the colors of branches
cols_branches <- c("purple", "forestgreen", "royalblue", "firebrick") # Set the colors of branches


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

