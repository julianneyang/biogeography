#####################Overlapping things for the immune deficiency subset 

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/ImmDef-Mouse-Biogeography-Analysis/")

#Purpose: create a dataframe with ASVs found only in both batches and call this df_overlap 

df_batch1 <- read.csv ("ImmDef-Overlap-Combat-ASV - s23-2014Nov.csv",row.names=1)
df_batch2 <- read.csv ("ImmDef-Overlap-Combat-ASV - s23-2014Sep.csv", row.names = 1)
df_batch3 <- read.csv ("ImmDef-Overlap-Combat-ASV - s7-2015Sep.csv", row.names=1)
df_batch4 <- read.csv ("ImmDef-Overlap-Combat-ASV - s9-2017Apr.csv", row.names=1)

df_overlap <- subset(df_batch1, row.names(df_batch1) %in% row.names(df_batch2)) 
df_overlap2 <- subset(df_batch2, row.names(df_batch2) %in% row.names(df_batch1))

#df_overlap <- as.data.frame(c(Sort(df_overlap, row.names)))

row.names(df_overlap)==row.names(df_overlap2)

bound_df <- cbind(df_overlap,df_overlap2)
df_overlap3 <- subset(df_batch3, row.names(df_batch3) %in% row.names(bound_df))
df_overlap4 <- subset(bound_df, row.names(bound_df) %in% row.names(df_overlap3))

row.names(df_overlap3)==row.names(df_overlap4)
bound_df2 <- cbind(df_overlap3,df_overlap4)

df_overlap6 <- subset(df_batch4, row.names(df_batch4) %in% row.names(bound_df2))
df_overlap5 <- subset(bound_df2, row.names(bound_df2) %in% row.names(df_overlap6))


row.names(df_overlap5)==row.names(df_overlap6)
finalbound <- cbind(df_overlap5,df_overlap6)


write.csv(finalbound, "ImmDef-Overlap-ASV.csv")

==
