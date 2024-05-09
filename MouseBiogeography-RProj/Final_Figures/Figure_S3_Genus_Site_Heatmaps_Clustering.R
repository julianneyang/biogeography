library(plyr)
library(rlang)
library(here)
library(dplyr)
library(tidyr)
library(gplots)
library(dendextend)
library(ggplot2)
library(cowplot)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_4_Genus_Site_Heatmaps_Clustering.R")

### Upset Plot ---

file_paths <- c("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                "UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite-1-MsID/all_results.tsv",
                "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")

cohort_prefixes <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_V_Gavage",
                     "UCLA_V_SPF",
                     "HUM_Gavage",
                     "SPF_Gavage")

all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes)

all_taxa <- all_taxa %>% select(c("feature", "Cohort")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)



df_wide <- as.data.frame(df_wide)
all_datasets <- names(df_wide)[-1]
taxa_upset <- ComplexUpset::upset(df_wide, all_datasets,
                                  base_annotations=list(
                                    'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
                                      scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')),
                                  themes=list(
                                    default=theme(
                                      axis.ticks.x=element_blank(),
                                      axis.text.x=element_blank(),
                                    ),
                                    intersections_matrix=theme(
                                      axis.ticks.x=element_blank(),
                                      axis.text.x=element_blank(),
                                    )
                                  ))

### Heatmap ---

## S3A. UCLA O. SPF Mucosal --
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


## CS SPF Mucosal --
path_to_significant_results <- "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv"
path_to_annotation_file <- "CS-Facility-Analysis/differential_genera_site/Genus_Mucosal_taxonomy.csv"
path_to_all_results <- "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
matrix.data <- results$matrix
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("royalblue", "firebrick","purple" , "forestgreen") # Set the colors of branches

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

## HUM V Gavage Mucosal --
path_to_significant_results <- "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"
path_to_annotation_file <- "Donors-Analysis/differential_genera_site/Genus_Mucosal_taxonomy.csv"
path_to_all_results <- "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
matrix.data <- results$matrix
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("royalblue", "firebrick", "forestgreen", "purple") # Set the colors of branches

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

## UCLA V SPF Mucosal --
path_to_significant_results <- "UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite-1-MsID/significant_results.tsv"
path_to_annotation_file <- "UCLA_V_SPF_Analysis/differential_genera_site/Genus_taxonomy.csv"
path_to_all_results <- "UCLA_V_SPF_Analysis/differential_genera_site/L6-DCvsAll-CLR-Muc-SeqRunSexSite-1-MsID/all_results.tsv"

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

## HUM Gavage Mucosal --
path_to_significant_results <- "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv"
path_to_annotation_file <- "Humanized-Biogeography-Analysis/differential_genera_site/HUM_Genus_Mucosal_taxonomy.csv"
path_to_all_results <- "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
data <- results$dataframe
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 
y = tapply(data$coef_d, data$Phylum_Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Phylum_Genus= factor(as.character(data$Phylum_Genus), levels = names(y))

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data,aes(x = value, y=Phylum_Genus)) + 
  geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = coul) +
  theme_cowplot(12) +
  theme(legend.position="none") +
  xlab("")+
  ylab("") +
  guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))+
  ggtitle("")

## SPF Gavage Mucosal --
path_to_significant_results <- "Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv"
path_to_annotation_file <- "Humanized-Biogeography-Analysis/differential_genera_site/SPF_Genus_Mucosal_taxonomy.csv"
path_to_all_results <- "Humanized-Biogeography-Analysis/differential_genera_site/SPF_L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv"

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
results <- generate_matrix_for_heatmap_clustering(path_to_significant_results = path_to_significant_results,
                                                  path_to_annotation_file = path_to_annotation_file,
                                                  path_to_all_results = path_to_all_results)
data <- results$dataframe
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 
y = tapply(data$coef_d, data$Phylum_Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Phylum_Genus= factor(as.character(data$Phylum_Genus), levels = names(y))

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(data,aes(x = value, y=Phylum_Genus)) + 
  geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = coul) +
  theme_cowplot(12) +
  theme(legend.position="none") +
  xlab("")+
  ylab("") +
  guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))+
  ggtitle("")
