library(Maaslin2)
library(funrar)
library(dplyr)
library(tidyr)
library(here)

here::i_am("MouseBiogeography-RProj/Maaslin2_L2_SITE_all.R")
here::here()

### Site Comparisons at the Phylum Level---
## HUM V Gavage --
input_data <- read.delim(here("Donors-Analysis/site_subsets/export_L2_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), header=TRUE, row.names=1) 
df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))
names(df_input_data) <- gsub("X","", names(df_input_data))

input_metadata <- read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"), header=TRUE) 
input_metadata$SampleID <- gsub("-",".",input_metadata$SampleID)
row.names(input_metadata) <- input_metadata$SampleID

target <- names(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
df_input_metadata$MouseID<- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Site<- factor(df_input_metadata$Site)
df_input_metadata$Site_General<- factor(df_input_metadata$Site_General)

# Luminal 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site_General,Colon'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site,Distal_Colon'))

# Mucosal -
input_data <- read.delim(here("Donors-Analysis/site_subsets/export_L2_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), header=TRUE, row.names=1) 
df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))
names(df_input_data) <- gsub("X","", names(df_input_data))

input_metadata <- read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"), header=TRUE) 
input_metadata$SampleID <- gsub("-",".",input_metadata$SampleID)
row.names(input_metadata) <- input_metadata$SampleID

target <- names(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Donor_ID<- factor(df_input_metadata$Donor_ID)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site_General,Colon'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Maaslin2_L2/HUM_V_Gavage/L2-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site,Distal_Colon'))


## UCLA O. SPF --
input_data <- read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/collapsed_taxa/export_L2_UCLA-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Mucosal_L2.csv", header=TRUE, row.names=1) 

df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))

input_metadata <-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/Regional-Metadata-All.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

sapply(df_input_metadata,levels)

#Luminal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-ColonRef-CLR-Lum-ComBat-SeqRunLineSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID_Line"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

## UCLA V. SPF --
input_data <- read.csv("Maaslin2_L2/UCLA_V_SPF/Mucosal_L2.csv", header=TRUE, row.names=1) 
df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))
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

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_V_SPF/L2_ColonRef-CLR-Muc-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_V_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

## CS SPF -- 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Site_L2/Maaslin2_L2 SITE and TYPE  - Luminal_L2.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_L2/CS_SPF/Site_L2/Maaslin2_L2 SITE and TYPE  - Mucosal_L2.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv",header=TRUE, row.names=1) #mapping file
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)
names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Luminal 
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal 
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/Site_L2/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

## HUM Gavage --
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Luminal_L2.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Mucosal_L2.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)
names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

#Luminal
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Source RPCA/Hum/Maaslin2 Site Genus/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"),
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

## SPF Gavage --
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Luminal.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Mucosal.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

#Luminal Site
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Site
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
R.version
