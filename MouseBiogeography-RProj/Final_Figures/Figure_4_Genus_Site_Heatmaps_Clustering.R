library(gplots)
library(dendextend) 
library(here)
library(ComplexUpset)
library(UpSetR)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_4_Genus_Site_Heatmaps_Clustering.R")

### Functions ---
process_results_for_upset_plot <- function(file_paths, cohort_prefixes) {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site" & qval<0.05)
    
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    data <- data %>% mutate(Cohort = cohort)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  
  return(data_all)
}
### UpSet Plot --

file_paths <- c("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                "CS-Facility-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                "Humanized-Biogeography-Analysis/differential_genera_site/HUM_L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")

cohort_prefixes <- c("UCLA_O_SPF_Luminal",
                     "CS_SPF_Luminal",
                     "HUM_V_Gavage_Luminal",
                     "HUM_Gavage_Luminal")

all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                                      cohort_prefixes = cohort_prefixes)

all_taxa <- all_taxa %>% select(c("feature", "Cohort")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

df_wide <- as.data.frame(df_wide)
df_wide <- df_wide %>% mutate(SPF_Gavage_Luminal = 0)
all_datasets <- names(df_wide)[-1]
taxa_upset <- ComplexUpset::upset(df_wide, all_datasets,
                                  base_annotations=list(
  'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
      scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')))+
    theme_cowplot(12)
UpSetR::upset(df_wide, sets= all_datasets,main.bar.color = "skyblue", order.by = "freq")

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

