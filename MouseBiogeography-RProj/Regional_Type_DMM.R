library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggplot2)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(cluster)
library(factoextra)
library(purrr)
BiocManager::install("DirichletMultinomial")

here::i_am("MouseBiogeography-RProj/Regional_Type_DMM.R")
input_data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/Type Subsets_ KO Abundance Counts - Colon_Ko_Counts.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

#Grab the KOs that are significantly differentially abundant which map to GSEA 
significant_data <- read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/Maaslin2_KO_Type/KO-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID/significant_results.tsv", header=TRUE) 
significant_data <- significant_data %>% filter(qval<0.05 & metadata == "Type")
target <- significant_data$feature

df_input_data <- df_input_data %>% filter(row.names(df_input_data) %in% target)

#run DMM modeling
count <- as.matrix(t(df_input_data))
set.seed(123)

# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(count, k, nstart = 10 )$tot.withinss
}

# Compute and plot wss for k = 1 to k = 15
k.values <- 1:15

# extract wss for 2-15 clusters
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")

#apply DMM Modeling
fit <- lapply(1:3, dmn, count = count, verbose=TRUE)

lplc <- base::sapply(fit, DirichletMultinomial::laplace) # AIC / BIC / Laplace
aic  <- base::sapply(fit, DirichletMultinomial::AIC) # AIC / BIC / Laplace
bic  <- base::sapply(fit, DirichletMultinomial::BIC) # AIC / BIC / Laplace
#plot(lplc, type="b", xlab="Number of Dirichlet Components", ylab="Model Fit")
#lines(aic, type="b", lty = 2)
#lines(bic, type="b", lty = 3)

best <- fit[[which.min(unlist(lplc))]]

#Sample Assignments
assn <- apply(mixture(best), 1, which.max)

#Contribution of each taxonomic group to each component

for (k in seq(ncol(fitted(best)))) {
  d <- melt(fitted(best))
  colnames(d) <- c("OTU", "cluster", "value")
  d <- subset(d, cluster == k) %>%
    # Arrange OTUs by assignment strength
    arrange(value) %>%
    mutate(OTU = factor(OTU, levels = unique(OTU))) %>%
    # Only show the most important drivers
    filter(abs(value) > quantile(abs(value), 0.8))     
  
  p <- ggplot(d, aes(x = OTU, y = value)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(title = paste("Top drivers: community type", k))
  print(p)
}