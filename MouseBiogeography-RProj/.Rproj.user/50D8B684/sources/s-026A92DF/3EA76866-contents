library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis")

input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - Cecum_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - Colon_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - SI_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - Duodenum_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - Ileum_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - Jejunum_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - ProximalColon_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - DistalColon_ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_LuminalvsMucosal/Type Subsets_ All Combat-Adjusted ASVs - DistalColon_ASV_names.csv", header=TRUE,row.names=1) # choose filtered non rarefied csv file


df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Maaslin2_LuminalvsMucosal/Regional-Metadata-All.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
#df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Mucosal", "Luminal"))
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

sapply(df_input_metadata,levels)

#Cecum 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Duodenum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Duodenum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Ileum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Ileum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Jejunum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Jejunum-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Proximal
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "MucRef-CLR-Proximal-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Proximal-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Distal
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "2MucRef-CLR-Distal-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LumRef-CLR-Distal-ComBat-SeqRunLineSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex","Type"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)
write.csv(Relative_Abundance,"Relative_Abundance-SI-ComBat.csv")
