library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/")

#Purpose: create a dataframe with ASVs found only in both batches and call this df_overlap 

df_batch1 <- read.csv ("s24_min10k_Novaseq_Jan2020.csv",row.names=1)
df_batch2 <- read.csv ("s3_min10k_Novaseq_Mar2020.tsv.csv", row.names = 1)
df_batch3 <- read.csv ("s57_min10k_Hiseq_2019.csv", row.names=1)

df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

row.names(df_overlap3)==row.names(df_overlap4)
finalbound <- cbind(df_overlap3,df_overlap4)

write.csv(finalbound, "Regional-Overlap-ASV.csv")
