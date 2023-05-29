library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Type_Modules_Results/Heatmap/")

#Purpose: create a dataframe with ASVs found only in both batches and call this df_overlap 

df_batch1 <- read.csv ("OMIXER TYPE Module Test Results - Cecum_MucosalvsLuminal.csv",row.names=2)
#df_batch1 <- read.csv ("OMIXER TYPE Module Test Results - Duodenum_Mucosal_vs_Luminal.csv",row.names=2)
  df_batch1$Feature <- row.names(df_batch1)
  df_batch1<-filter(df_batch1, FDR < 0.05)

#df_batch2 <- read.csv ("OMIXER TYPE Module Test Results - Jejunum_Mucosal_vs_Luminal.csv", row.names = 2)
df_batch2 <- read.csv ("OMIXER TYPE Module Test Results - DC_MucosalvsLuminal.csv", row.names = 2)
  df_batch2$Feature <- row.names(df_batch2)
  df_batch2<-filter(df_batch2, FDR < 0.05)
  
df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))
df_overlap2 <- as.data.frame(c(Sort(df_overlap2, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
row.names(bound_df) <- bound_df$Feature

df_batch3 <- read.csv ("OMIXER TYPE Module Test Results - PC_MucosalvsLuminal.csv", row.names = 2)
#df_batch3 <- read.csv ("OMIXER TYPE Module Test Results - Ileum_Mucosal_vs_Luminal.csv", row.names = 2)
  df_batch3$Feature <- row.names(df_batch3)
  df_batch3 <- filter(df_batch3, FDR < 0.05)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

df_overlap3 <- as.data.frame(c(Sort(df_overlap3, row.names)))
df_overlap4 <- as.data.frame(c(Sort(df_overlap4, row.names)))

row.names(df_overlap3)==row.names(df_overlap4)

####
finalboundSI <- cbind(df_overlap3,df_overlap4)
finalboundCOL <- cbind(df_overlap3,df_overlap4)
row.names(finalboundCOL) <- finalboundCOL$Feature

df_overlap <- subset(finalboundSI, row.names(finalboundSI) %in% row.names(finalboundCOL)) 
  row.names(df_overlap) <- df_overlap$Feature
df_overlap2 <- subset(finalboundCOL, row.names(finalboundCOL) %in% row.names(finalboundSI))
  row.names(df_overlap2) <-df_overlap2$Feature
df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))
df_overlap2 <- as.data.frame(c(Sort(df_overlap2, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

finalfinalbound <- cbind(df_overlap,df_overlap2)

write.csv(finalfinalbound, "Combined_SIandColon_OMIXER_MucosalvsLuminal.csv")


