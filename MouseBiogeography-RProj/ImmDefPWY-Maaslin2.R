rm(list = ls())
library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/ImmDef-Mouse-Biogeography-Analysis/")


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
input_data <- read.csv("PathwayCounts-AllSubsets - Luminal-PWY.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

input_data <- read.csv("ImmDef-PostCombat-PWY - ImmDef-PWY-ASV.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file


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
input_metadata <-read.csv("Metadata-Luminal.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("Metadata-All Subsets - Mucosal.csv",header=TRUE, row.names=1) #mapping file
input_metadata <-read.csv("ImmDef-PostCombat-PWY - ImmDef-PWY-metadata.csv",header=TRUE, row.names=1) #mapping file


row.names(input_metadata)==colnames(df_input_data)

df_input_metadata <- as.data.frame(input_metadata)

df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Genotype <- factor(df_input_metadata$Genotype, "WT", "TCR_KO", "RAGROR")
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal","Mucosal"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon","Proximal_Colon" , "Cecum", "Ileum", "Jejunum", "Duodenum"))

#df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon","Cecum"))
#$Site <- factor(df_input_metadata$Site, levels=c("Ileum", "Jejunum","Duodenum"))
sapply(df_input_metadata,levels)

warnings()
#fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-LumCol-ComBat-Maaslin-SeqRunSexSite_General-1-Line-MsID", fixed_effects = c("SequencingRun","Sex", "Site"), random_effects = c("Line", "MouseIDLine"),normalization="clr", transform ="LOG")
#fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-LumSI-ComBat-Maaslin-SeqRunSexSite-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
#fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-MucCol-ComBat-Maaslin-SeqRunSexSite-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
#fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-MucSI-ComBat-Maaslin-SeqRunSexSite-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Duodenum-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Ileum-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Jejunum-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Cecum-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-ProximalColon-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-DistalColon-ComBat-Maaslin-SeqRunSexType-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Luminal-ComBat-Maaslin-SeqRunSexSite_General-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Mucosal-ComBat-Maaslin-SeqRunSexSite_General-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Luminal-ComBat-Maaslin-SeqRunSexSite-1-Line-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("Line", "MouseID_Line"),normalization="clr", transform ="LOG")
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "LogCLR-Mucosal-ComBat-Maaslin-SeqRunSexSiteGenotype-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "GenotypeRAGROR.SiteDuodenum", "GenotypeTCR_KO.SiteDuodenum", "GenotypeRAGROR.SiteJejunum", "GenotypeTCR_KO.SiteJejunum", "GenotypeRAGROR.SiteIleum", "GenotypeTCR_KO.SiteIleum","GenotypeRAGROR.SiteCecum", "GenotypeTCR_KO.SiteCecum", "GenotypeRAGROR.SiteProximal_Colon", "GenotypeTCR_KO.SiteProximal_Colon"), random_effects = c("MouseID"),normalization="clr", transform ="LOG")


transposed_input_data <- t(df_input_data)
transposed_input_data <- as.matrix(transposed_input_data) #taxa are now columns, samples are rows. 
df_relative_ASV <- make_relative(transposed_input_data)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)

write.csv(Relative_Abundance,"Relative_Abundance-Duodenum-ComBat.csv")

