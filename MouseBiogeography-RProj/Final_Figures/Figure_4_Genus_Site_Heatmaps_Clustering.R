library(gplots)
library(dendextend) 
library(here)
library(ComplexUpset)
library(tidyverse)
library(UpSetR)
library(cowplot)

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

cohort_prefixes <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_MD_Gavage",
                     "HUM_SD_Gavage")

all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                                      cohort_prefixes = cohort_prefixes)

id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","Cohort","coef_dir")) %>% unique()

id_features <- all_taxa 
id_features <- id_features%>% select(c("feature","Cohort")) %>% unique()

id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)
id_df_wide <- id_df_wide %>% mutate(SPF_Gavage = 0)

all_taxa <- all_taxa %>% select(c("feature", "Cohort")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

df_wide <- as.data.frame(df_wide)
df_wide <- df_wide %>% mutate(SPF_Gavage = 0)
all_datasets <- names(df_wide)[-1]
taxa_upset <- ComplexUpset::upset(df_wide, 
                                  all_datasets,
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

id_df_wide$count_ones <- rowSums(id_df_wide[, c(3:7)])
df_filtered <- id_df_wide[id_df_wide$count_ones >= 3, ]
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
df_filtered$feature<-gsub(".*f__","f__",df_filtered$feature)
df_filtered$feature

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


### Shotgun barplots ---

## UCLA O SPF
result2 <- generate_interregional_taxa_barplot_shotgun_only_named_species(
  path_to_significant_results_tsv = "Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_LineSexSite-1-MsID/significant_results.tsv",
  titlestring="UCLA O. SPF",
  colorvector = cols)

df <- result2$dataframe
phylum_names <- df$Phylum

select_cols <- c("Firmicutes"="#aa0000ff", "Bacteroidetes"="#800080ff","Actinobacteria"="#008000ff",
                 "Bacteria_unclassified"="black", "Candidatus_Saccharibacteria"="#808000ff","Proteobacteria"="#00ffffff")
seecolor::print_color(select_cols)
phylum_colors <- select_cols
names(select_cols)

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

ucla_o_shotgun_species <- result2$plot+
  theme(axis.text.y = element_text(colour = color_mapping))+
  theme(legend.position = "right")
ucla_o_shotgun_species

## CS SPF 

cs_result <- generate_interregional_taxa_barplot_shotgun_only_named_species(path_to_significant_results_tsv = "Shotgun/CS_SPF/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                      titlestring="CS SPF",
                                                      colorvector = cols)

df <- cs_result$dataframe 
phylum_names <- df$Phylum

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

cs_shotgun_species <- cs_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
cs_shotgun_species

## SPF Gavage 
cols=c("#FDE725FF")

spf_result <- generate_interregional_taxa_barplot_shotgun_only_named_species(path_to_significant_results_tsv = "Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                       titlestring="SPF Gavage",
                                                       colorvector = cols)

df <- spf_result$dataframe 
phylum_names <- df$Phylum


# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

spf_shotgun_species <- spf_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
spf_shotgun_species

## HUM Gavage --
cols=c("#440154FF", "#FDE725FF")

hum_result <- generate_interregional_taxa_barplot_shotgun_only_named_species(path_to_significant_results_tsv = "Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                       titlestring="HUM SD Gavage",
                                                       colorvector = cols)

df <- hum_result$dataframe 
phylum_names <- df$Phylum

color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

hum_shotgun_species <- hum_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
hum_shotgun_species 

### Final Figure bottom half
ucla_o_shotgun_species<-ucla_o_shotgun_species+theme(legend.position="none")
plot_grid(taxa_upset, ucla_o_shotgun_species,
          labels=c("F","G"),
          label_size = 20)
plot_grid(cs_shotgun_species, hum_shotgun_species, spf_shotgun_species,
          labels=c("H","I","J"),nrow=3,ncol=1,
          rel_heights = c(1,0.8, 0.4),
          label_size = 20)

### Functions---
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
  data$value = plyr::revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
  data$value = factor(data$value, levels=c("D", "J", "I", "C", "PC", "DC"))
  ggplot_data <- unique(data)
  ggplot_data$Phylum_Genus<-paste(ggplot_data$Phylum,ggplot_data$annotation,sep=" : ")
  
  
  #construct heatmap using heatmap2 with dendrogram
  data_long<-pivot_wider(ggplot_data, id_cols=annotation, names_from = value, values_from =coef_d)
  data_long_final<-data_long[,-1]
  data_long_final<-select(data_long_final,D,J, I,C,PC,DC)
  row.names(data_long_final)= data_long$annotation
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
