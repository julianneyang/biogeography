library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction")
input_data <- read.csv("Site Subsets_ KO Abundance Counts - Luminal_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets_ KO Abundance Counts - LuminalColon_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets_ KO Abundance Counts - LuminalSI_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets_ KO Abundance Counts - Mucosal_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets_ KO Abundance Counts - MucosalColon_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Site Subsets_ KO Abundance Counts - MucosalSI_KO.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file


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
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO-SIRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Luminal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO-CecRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO-DCRef-CLR-LumCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Luminal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "DuodRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "IleumRef-CLR-LumSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_SIRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Colon
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Cecum", "Proximal_Colon", "Distal_Colon"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_CecRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_DCRef-CLR-MucCol-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal SI
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_DuodRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "KO_IleumRef-CLR-MucSI-ComBat-SeqRunLineSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-Mucosal-KO.csv")
