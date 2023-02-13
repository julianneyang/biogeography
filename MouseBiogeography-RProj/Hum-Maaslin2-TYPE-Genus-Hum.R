library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - DC_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - PC_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - Cec_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - Ile_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - Jej_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - Duo_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - Colon_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Source RPCA/Hum/Maaslin2_TYPE_Genus/Humanized_Maaslin2_L6 SITE and TYPE - SI_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

sapply(df_input_metadata,levels)

#Distal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

#Proximal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE )

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Duodenum-ComBat-SexType-1-MsID", fixed_effects = c("Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE, plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2_TYPE_Genus/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-SI-L6.csv")
