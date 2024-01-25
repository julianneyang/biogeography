library(Maaslin2)
library(funrar)
library(dplyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/")
first_filepath <- "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/"
input_data <- readr::read_delim(here(paste0(first_filepath,"site_subsets/export_L6_Luminal_UCLA-ComBat-Adjusted-ASV/feature-table.tsv")), delim="\t") # choose filtered non rarefied csv file

input_data <- as.data.frame(input_data)
row.names(input_data)<-input_data$OTU.ID
df_input_data <- select(input_data, -c("taxonomy","OTU.ID"))

input_metadata <-readr::read_delim(here(paste0(first_filepath,"starting_files/Regional-Combat-Metadata.tsv")),delim="\t") #mapping file
input_metadata <- as.data.frame(input_metadata)
row.names(input_metadata) <- input_metadata$SampleID
input_metadata <- input_metadata %>% select(-c("SampleID"))

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
?Maaslin2

## Luminal --
outdir <- paste0(first_filepath,"differential_genera_site/")
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = paste0(outdir,"L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Hiseq_April_Nineteen','Site_General,Colon','Line,JJWT'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = paste0(outdir,"L6-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Hiseq_April_Nineteen','Site_General,Colon','Line,JJWT','Site,Distal_Colon'))

## Mucosal --
input_data <- readr::read_delim(here(paste0(first_filepath,"site_subsets/export_L6_Mucosal_UCLA-ComBat-Adjusted-ASV/feature-table.tsv")), delim="\t") # choose filtered non rarefied csv file
input_data <- as.data.frame(input_data)
row.names(input_data)<-input_data$OTU.ID
df_input_data <- select(input_data, -c("taxonomy","OTU.ID"))

input_metadata <-readr::read_delim(here(paste0(first_filepath,"starting_files/Regional-Combat-Metadata.tsv")),delim="\t") #mapping file
input_metadata <- as.data.frame(input_metadata)
row.names(input_metadata) <- input_metadata$SampleID
input_metadata <- input_metadata %>% select(-c("SampleID"))

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
?Maaslin2

## Mucosal --
outdir <- paste0(first_filepath,"differential_genera_site/")
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = paste0(outdir,"L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Hiseq_April_Nineteen','Site_General,Colon','Line,JJWT'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = paste0(outdir,"L6-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Hiseq_April_Nineteen','Site_General,Colon','Line,JJWT','Site,Distal_Colon'))
