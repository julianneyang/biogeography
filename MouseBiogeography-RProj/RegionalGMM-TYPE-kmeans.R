library(factoextra)
library(tidyr)
library(cluster)
library(fpc)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/")
duojejilececpcdc<-read.csv("GMM-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

annotation <- read.csv("Revised_Module_Key.csv", header=TRUE)
data<- (merge(gmm_heatmap_final, annotation, by = 'feature'))
data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")


data_wide<-pivot_wider(data, id_cols=hierachy_annotations, names_from = Site, values_from =coef)
names=data_wide$hierachy_annotations
data_wide<-data_wide[,-1]
row.names(data_wide)=names

?fviz_nbclust()
optim_clust<-fviz_nbclust(data_wide, FUNcluster=kmeans, method="wss")
dev.new(width=15, height=10)
optim_clust

# Compute k-means with k = 6 from optim clust
set.seed(123)
km.res <- kmeans(data_wide, 6, nstart = 25)

# Print the results
print(km.res)

# Cluster number for each of the observations
km.res$cluster
head(km.res$cluster, 6)
plotcluster(data_wide, km.res$cluster)
clusplot(data_wide, km.res$cluster, color=TRUE, shade=TRUE, 
         labels = 1, lines=0)

?fvis_cluster()
fvis_cluster(data_wide)
