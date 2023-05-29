rm(list = ls())
library(Maaslin2)
library(funrar)
if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Maaslin2")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2/")

set.seed(123)
input_data <- read.csv("LumCol ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("LumSI ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("MucCol ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("MucSI ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("LumCol Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("LumSI Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("MucCol Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("MucSI Metadata.csv",header=TRUE, row.names=1) #mapping file
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sex <- factor(df_input_metadata$Sex, levels = c("Male", "Female"))
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run, levels = c("Hiseq_April_Nineteen", "NovaSeq_Mar_Twenty", "NovaSeq_Jan_Twenty"))
#df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon","Cecum"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum","Jejunum","Duodenum"))
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
#df_input_metadata$Site_General <- factor(df_input_metadata$Site_General)
sapply(df_input_metadata,levels)


fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-MucosalSI-Maaslin-SeqRunSexSite_General-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG", cores=8)
#?Maaslin2

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)
write.csv(Relative_Abundance, "MucSI-Relative-ASV.csv")


