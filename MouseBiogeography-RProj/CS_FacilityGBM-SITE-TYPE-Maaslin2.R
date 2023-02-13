library(Maaslin2)
library(funrar)
library(dplyr)


setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/CS_FacilityGMM-SITE-TYPE-Maaslin2.R")
input_data <- read.table("CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/modules.tsv", header=TRUE, row.names=1) 
df_input_data<-as.data.frame(input_data)
input_metadata <-read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv",header=TRUE, row.names=1) #mapping file
row.names(input_metadata)

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata$SampleID <- row.names(input_metadata)
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)

## Get Luminal and Mucosal datasets to test DC vs all ---
lum_samples <- df_input_metadata %>% filter(Type =="Luminal" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
lum_input_data <- df_input_data[, lum_samples]

muc_samples <- df_input_metadata %>% filter(Type =="Mucosal" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
muc_input_data <- df_input_data[, muc_samples]

## Get Site datasets to test Muc vs Lum ---
si_samples <- df_input_metadata %>% filter(Site_General =="SI" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
si_input_data <- df_input_data[, si_samples]

col_samples <- df_input_metadata %>% filter(Site_General =="Colon" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
col_input_data <- df_input_data[, col_samples]

duo_samples <- df_input_metadata %>% filter(Site =="Duodenum" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
duo_input_data <- df_input_data[, duo_samples]

jej_samples <- df_input_metadata %>% filter(Site =="Jejunum" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
jej_input_data <- df_input_data[, jej_samples]

ile_samples <- df_input_metadata %>% filter(Site =="Ileum" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
ile_input_data <- df_input_data[, ile_samples]

cec_samples <- df_input_metadata %>% filter(Site =="Cecum" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
cec_input_data <- df_input_data[, cec_samples]

PC_samples <- df_input_metadata %>% filter(Site =="Proximal_Colon" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
PC_input_data <- df_input_data[, PC_samples]

DC_samples <- df_input_metadata %>% filter(Site =="Distal_Colon" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
DC_input_data <- df_input_data[, DC_samples]

## Run Maaslin2 ---

#Luminal 
target <- colnames(lum_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=lum_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=lum_input_data, input_metadata=df_input_metadata,
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
target <- colnames(muc_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=muc_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=muc_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

## Run Maaslin2 ---

#Distal_Colon -overfitting
target <- colnames(DC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=DC_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

#Proximal_Colon -overfitting
target <- colnames(PC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
target <- colnames(cec_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE )

#Ileum
target <- colnames(ile_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
target <- colnames(jej_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
target <- colnames(duo_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon
target <- colnames(col_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE, plot_scatter = FALSE)

#SI
target <- colnames(si_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)


