library(Maaslin2)
library(funrar)
library(dplyr)
library(tidyr)

here::i_am("MouseBiogeography-RProj/Maaslin2_L6_TYPE_all.R")
here::here()

### Function to run Maaslin2 analysis ---

run_maaslin2 <- function(input_data, output_path, site, site_fixed_effects, random, set_refs) {
  # Get target columns from input data
  target <- colnames({{input_data}})
  
  # Subset metadata for matching samples
  df_input_metadata_subset <- input_metadata[match(target, row.names(input_metadata)),]
  
  # Run Maaslin2
  fit_data <- Maaslin2(input_data = input_data, 
                       input_metadata = df_input_metadata_subset, 
                       output = {{output_path}}, 
                       fixed_effects = {{site_fixed_effects}}, 
                       random_effects = {{random}},
                       reference = {{set_refs}},
                       normalization = "clr", 
                       min_prevalence = 0,
                       transform = "none",
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
  
  return(fit_data)
}

### Type Comparisons at the Phylum Level---
## UCLA O. SPF --
input_data <- read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/export_L6_UCLA-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
rows_to_remove <- grep("Mitochondria|Chloroplast", row.names(df_input_data))
df_input_data <- df_input_data[-rows_to_remove, ]

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
  Duodenum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID",
  Jejunum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID",
  Ileum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID",
  Cecum = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID",
  Proximal_Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunLineSexType-1-MsID",
  Distal_Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunLineSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunLineSexSiteType-1-MsID",
  Colon = "Regional-Mouse-Biogeography-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunLineSexSiteType-1-MsID"
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
input_data <- read.delim(here("CS_SPF/starting_files/export_L6_CS-Facility-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
rows_to_remove <- grep("Mitochondria|Chloroplast", row.names(df_input_data))
df_input_data <- df_input_data[-rows_to_remove, ]
input_metadata <-read.delim(here("CS_SPF/starting_files/CS_Facility_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
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
  Duodenum = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID",
  Jejunum = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID",
  Colon = "CS-Facility-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"
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
input_data <- read.delim(here("Humanized-Biogeography-Analysis/starting_files/export_L6_min10000_Cedars_SPF_Colonized-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c(taxonomy))
input_metadata <-read.delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
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
  Duodenum = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID",
  Jejunum = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID",
  Colon = "Humanized-Biogeography-Analysis/differential_genera_type/SPF_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"
)

# Iterate over each site
site_fe <- c("Sequencing_Run", "Sex", "Type")
ranef <- c("MouseID")
refs <- c("Sequencing_Run,2014_Nov","Site,Distal_Colon")

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
input_metadata$Sequencing_Run <- factor(input_metadata$Sequencing_Run, levels =c("2014_Sept","2014_Nov", "2015_Sept"))
input_metadata$Site <- factor(input_metadata$Site, levels =c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
input_metadata$Type <- factor(input_metadata$Type, levels =c("Luminal", "Mucosal"))


for (site in names(site_general_paths)) {
  samples <- df_input_metadata %>% 
    filter(Site_General == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_general_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_general_fe,ranef,refs)
}


## HUM Gavage --
input_data <- read.delim(here("Humanized-Biogeography-Analysis/starting_files/export_L6_min10000_Humanized_Colonized-ComBat-Adjusted-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
rows_to_remove <- grep("Mitochondria|Chloroplast", row.names(df_input_data))
#df_input_data <- df_input_data[-rows_to_remove, ]

input_metadata <-read.delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
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

# HUM SD Gavage -- List of sites and corresponding file paths
site_paths <- list(
  Jejunum = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID",
  Colon = "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"
)


# Iterate over each site
site_fe <- c("Sequencing_Run", "Sex", "Type")
duod_fe <- c("Sex","Type")
ranef <- c("MouseID")
refs <- c("Sequencing_Run,2014_Nov","Site,Distal_Colon")

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

# For duodenum only because sex and seq run are conflated
samples <- df_input_metadata %>% 
  filter(Site == "Duodenum", SampleID %in% names(df_input_data)) %>% 
  pull(SampleID)

input_data <- df_input_data[, samples]

output_path <- "Humanized-Biogeography-Analysis/differential_genera_type/HUM_L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID"
fit_data <- run_maaslin2(input_data, output_path, site, duod_fe,ranef,refs)


site_general_fe <- c("Sequencing_Run","Sex", "Site","Type")
refs<- c("Sequencing_Run,2014_Nov","Site,Distal_Colon")

input_metadata$Sequencing_Run <- factor(input_metadata$Sequencing_Run, levels =c("2014_Sept","2014_Nov", "2015_Sept"))
input_metadata$Site <- factor(input_metadata$Site, levels =c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
input_metadata$Type <- factor(input_metadata$Type, levels =c("Luminal", "Mucosal"))


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

## HUM V Gavage --
input_data <- read.delim(here("Donors-Analysis/starting_files/collapsed_taxa/export_L6_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
rows_to_remove <- grep("Mitochondria|Chloroplast", row.names(df_input_data))
df_input_data <- df_input_data[-rows_to_remove, ]
input_metadata <-read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
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
df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)

site_paths <- list(
  Duodenum = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID",
  Jejunum = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-PC-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-DC-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID",
  Colon = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"
)

# Iterate over each site
site_fe <- c("Sequencing_Run", "Sex", "Type")
ranef <- c("MouseID")
refs <- c("Sequencing_Run,Jan_2017", "Type,Luminal","Site,Distal_Colon")

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
