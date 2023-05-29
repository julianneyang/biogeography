rm(list = ls())
library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/Maaslin2-All-Intestinal-Sites/")


input_data <- read.csv("PathwayCounts-AllSites - LuminalColon-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("PathwayCounts-AllSites - LuminalSI-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("PathwayCounts-AllSites - MucosalColon-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("PathwayCounts-AllSites - MucosalSI-PWY (1).csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("duodenum-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("ileum-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("jejunum-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("cecum-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("proxcolon-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("distcolon-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file



df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("LumCol-Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-All Sites - LumSI.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-All Sites - MucCol.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-All Sites - MucSI.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-Duodenum.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-Ileum.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-Jejunum.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-Cecum.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-ProximalColon.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-DistalColon.csv",header=TRUE, row.names=1) #mapping file

row.names(input_metadata)==colnames(df_input_data)

df_input_metadata <- as.data.frame(input_metadata)

df_input_metadata$Sex <- factor(df_input_metadata$Sex)


write.csv(Relative_Abundance,"Relative_Abundance-Duodenum-ComBat.csv")

