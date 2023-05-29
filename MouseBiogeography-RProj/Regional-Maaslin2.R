
library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysi21-8-Microbiome-Batch-Correction-Analysis")

input_data <- read.csv("Maaslin2/LumCol ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2/LumSI ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2/MucCol ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2/MucSI ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2/All Subsets_ Combat Adjusted ASV and Metadata  - Luminal ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2/All Subsets_ Combat Adjusted ASV and Metadata  -  Mucosal ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file


df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Maaslin2/LumCol Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Maaslin2/LumSI Metadata.csv",header=TRUE, row.names=1) #mapping files/20
input_metadata <-read.csv("Maaslin2/MucCol Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Maaslin2/MucSI Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Maaslin2/All Subsets_ Combat Adjusted ASV and Metadata  - Luminal Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Maaslin2/All Subsets_ Combat Adjusted ASV and Metadata  - Mucosal Metadata.csv",header=TRUE, row.names=1) #mapping file

row.names(input_metadata)==colnames(input_data)

df_input_metadata <- as.data.frame(input_metadata)

df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
sapply(df_input_metadata,levels)

#Luminal 
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "SIRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "SIRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Luminal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CecRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "DCRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Luminal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "DuodRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "IleumRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Mucosal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CecRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "DCRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

#Mucosal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "DuodRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "IleumRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none")

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-Mucosal-ComBat.csv")
write.csv(Relative_Abundance,"Relative_Abundance-Luminal-ComBat.csv")
?Maaslin2

#Maaslin version- 1.2.0