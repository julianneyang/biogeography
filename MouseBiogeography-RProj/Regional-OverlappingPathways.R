library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/Pathway Heatmaps/LuminalvMucosal/")

#Purpose: create a dataframe with ASVs found only in both batches and call this df_overlap 

df_batch1 <- read.csv ("Pathway Heatmaps/Pathways-Maaslin2-MucosalSI - DuodenumMucosalSI-PWY-heatmap.csv",row.names=1)
df_batch2 <- read.csv ("Pathway Heatmaps/Pathways-Maaslin2-MucosalSI - Jejunum-MucosalSI-PWY-heatmap.csv", row.names = 1)

df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

#df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

row.names(df_overlap3)==row.names(df_overlap4)
finalbound <- cbind(df_overlap3,df_overlap4)

write.csv(finalbound, "Regional-Overlap-ASV.csv")
write.csv(bound_df, "MucosalSI-Overlap.csv")

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/")


#For Luminal vs Mucosal data in each of six intestinal sites: 

#df_batch1 <- read.csv ("Heatmap-LuminalvsMucosal-Colon - Cecum.csv",row.names=1)
#df_batch2 <- read.csv ("Heatmap-LuminalvsMucosal-Colon - Distal_Colon.csv", row.names = 1)
#df_batch3 <- read.csv ("Heatmap-LuminalvsMucosal-Colon - Proximal_Colon.csv", row.names=1)

df_batch1 <- read.csv ("Heatmap-LuminalvsMucosal-SI - Duodenum.csv",row.names=1)
df_batch2 <- read.csv ("Heatmap-LuminalvsMucosal-SI - Ileum.csv", row.names = 1)
df_batch3 <- read.csv ("Heatmap-LuminalvsMucosal-SI - Jejunum.csv", row.names=1)

df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

#df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

row.names(df_overlap3)==row.names(df_overlap4)
finalbound <- cbind(df_overlap3,df_overlap4)

write.csv(finalbound, "Combined-SI-LuminalvsMucosal.csv")


#Purpose: create a dataframe with ASVs found only in both batches and call this df_overlap 

df_batch1 <- read.csv ("Pathway Heatmaps/Pathways-Maaslin2-MucosalSI - DuodenumMucosalSI-PWY-heatmap.csv",row.names=1)
df_batch2 <- read.csv ("Pathway Heatmaps/Pathways-Maaslin2-MucosalSI - Jejunum-MucosalSI-PWY-heatmap.csv", row.names = 1)
df_batch3 <- read.csv ("s57_min10k_Hiseq_2019.csv", row.names=1)

df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

#df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

row.names(df_overlap3)==row.names(df_overlap4)
finalbound <- cbind(df_overlap3,df_overlap4)

write.csv(finalbound, "Regional-Overlap-ASV.csv")
write.csv(bound_df, "MucosalSI-Overlap.csv")

