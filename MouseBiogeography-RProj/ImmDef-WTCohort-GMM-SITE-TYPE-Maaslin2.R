library(Maaslin2)
library(funrar)
library(dplyr)
library(tidyr)

here::i_am("MouseBiogeography-RProj/ImmDef-WTCohort-GMM-SITE-TYPE-Maaslin2.R")
here::here()

input_data <- read.table("ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/modules.tsv", header=TRUE, row.names=1) 
df_input_data<-as.data.frame(input_data)
input_metadata <-read.csv("ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv",header=TRUE, row.names=1) #mapping file
input_metadata$SampleID <- row.names(input_metadata)

target <- colnames(df_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)


df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID_Original)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)

## Get Mucosal datasets to test DC vs all ---
muc_samples <- df_input_metadata %>% filter(Type =="Mucosal" & Genotype =="WT", SampleID %in% names(df_input_data)) %>% pull(SampleID)
muc_input_data <- df_input_data[, muc_samples]

## Run Maaslin2 ---
#Mucosal
target <- colnames(muc_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=muc_input_data, input_metadata=df_input_metadata, 
                    output = "ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM_ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=muc_input_data, input_metadata=df_input_metadata, 
                    output = "ImmDef-Mouse-Biogeography-Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)


