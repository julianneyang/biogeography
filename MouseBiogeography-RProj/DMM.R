##Dirichlet Multinomial Mixed Models
#Retrieved from: https://microbiome.github.io/tutorials/DMM.html
library(ggplot2)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(cluster)
library(factoextra)
library(purrr)
BiocManager::install("DirichletMultinomial")

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/")
input_data <- read.csv("Type Subsets_ KO Abundance Counts - Colon_Ko_Counts.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

#apply prevalence filter 
t_df_input_data<-as.data.frame(t(df_input_data))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}

df_input_data$prevalence<-prevalence
df_input_data <- df_input_data %>% filter(prevalence>42)

#run DMM modeling
df_input_data<-select(df_input_data, -prevalence)
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
