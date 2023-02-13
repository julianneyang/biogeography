library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/")
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_Luminal.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_Mucosal.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_LumCol.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_LumSI.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_MucCol.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets- Combat-Adjusted Genera.xlsx - genus_MucSI.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

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

#Luminal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-SIRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum","Proximal_Colon","Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DuodvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Luminal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-CecRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DCRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Luminal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DuodRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-IleumRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_SIRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum","Proximal_Colon","Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6-DuodvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_CecRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_DCRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_DuodRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "L6_IleumRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-Mucosal-L6.csv")
