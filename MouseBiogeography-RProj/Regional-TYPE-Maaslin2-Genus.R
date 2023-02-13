library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/")
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_Colon.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_SI.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_DistalColon.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_ProximalColon.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_Cecum.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_Ileum.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_Jejunum.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Type Subsets- Combat-Adjusted Genera.xlsx - genus_Duodenum.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Regional-Metadata-All.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

?Maaslin2

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
#df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Mucosal", "Luminal"))
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

sapply(df_input_metadata,levels)

#Distal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

#Proximal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE )

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE, plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-SI-L6.csv")
