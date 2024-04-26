library(Maaslin2)
library(funrar)
library(dplyr)
library(tidyr)

here::i_am("MouseBiogeography-RProj/Maaslin2_L2_TYPE_all.R")
here::here()

### Function to run Maaslin2 analysis ---

run_maaslin2 <- function(input_data, output_path, site, site_fixed_effects, random, set_refs) {
  # Get target columns from input data
  target <- colnames(input_data)
  
  # Subset metadata for matching samples
  df_input_metadata_subset <- input_metadata[match(target, row.names(input_metadata)),]
  
  # Run Maaslin2
  fit_data <- Maaslin2(input_data = input_data, 
                       input_metadata = df_input_metadata_subset, 
                       output = output_path, 
                       fixed_effects = {{site_fixed_effects}}, 
                       random_effects = {{random}},
                       min_prevalence = 0.15,
                       reference = {{set_refs}},
                       normalization = "clr", 
                       transform = "none",
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
  
  return(fit_data)
}

### Type Comparisons at the Phylum Level---
## UCLA O. SPF --
input_data <- read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/export_L2_UCLA-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
input_metadata <-read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

samples<- input_metadata %>%
  filter(SampleID %in% names(df_input_data)) %>%
  pull(SampleID)

df_input_data <- df_input_data[,samples]
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata$SampleID <- row.names(input_metadata)
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)


# List of sites and corresponding file paths

site_paths <- list(
  Duodenum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID",
  Jejunum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID",
  Ileum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID",
  Cecum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID",
  Proximal_Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID",
  Distal_Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID",
  Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L2-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID"
)


# Iterate over each site
site_fe <- c("Sequencing_Run", "Line","Sex", "Type")
ranef <- c("MouseID_Line")
refs <- c("Sequencing_Run,Hiseq_April_Nineteen","Line,ItgCre","Site,Distal_Colon")

for (site in names(site_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_fe,ranef,refs)
}

site_general_fe <- c("Sequencing_Run","Line", "Sex", "Site","Type")
for (site in names(site_general_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site_General == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_general_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_general_fe,ranef,refs)
}

## CS SPF --
input_data <- read.delim(here("CS-Facility-Analysis/starting_files/export_L2_CS-Facility-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
input_metadata <-read.delim(here("CS-Facility-Analysis/starting_files/CS_Facility_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

samples<- input_metadata %>%
  filter(SampleID %in% names(df_input_data)) %>%
  pull(SampleID)

df_input_data <- df_input_data[,samples]
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata$SampleID <- row.names(input_metadata)
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Site <- factor(df_input_metadata$Site)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)

# List of sites and corresponding file paths
site_paths <- list(
  Duodenum = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID",
  Jejunum = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID",
  Colon = "CS-Facility-Analysis/differential_genera_type/L2-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"
)

# Iterate over each site
site_fe <- c("Sequencing_Run", "Sex", "Type")
ranef <- c("MouseID")
refs <- c("Sequencing_Run,One","Site,Distal_Colon")

for (site in names(site_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_fe,ranef,refs)
}

site_general_fe <- c("Sequencing_Run","Sex", "Site","Type")
for (site in names(site_general_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site_General == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_general_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_general_fe,ranef,refs)
}

## SPF Gavage --
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Colon.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-SI.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Duo.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Jejunum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Ileum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Cecum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-PC.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-DC.csv", header=TRUE, row.names=1) 

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
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Colon-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-SI-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
## HUM Gavage --
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Colon_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - SI_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Duodenum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Jejunum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Ileum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Cecum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - PC_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - DC_L2.csv", header=TRUE, row.names=1) 

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
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Colon-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-SI-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Duodenum-CLR-ComBat-SexType-1-MsID", 
                    fixed_effects = c("Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
