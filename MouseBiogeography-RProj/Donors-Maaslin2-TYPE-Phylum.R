library(Maaslin2)
library(funrar)
library(dplyr)
library(here)

here::i_am("MouseBiogeography-RProj/Donors-Maaslin2-TYPE-Genus.R")

input_data <- read.delim(here("Donors-Analysis/starting_files/collapsed_taxa/export_L2_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
df_input_data <- df_input_data %>% select(-c(taxonomy))
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



# Function to run Maaslin2 analysis

run_maaslin2 <- function(input_data, output_path, site, site_fixed_effects) {
  # Get target columns from input data
  target <- colnames(input_data)
  
  # Subset metadata for matching samples
  df_input_metadata_subset <- input_metadata[match(target, row.names(input_metadata)),]
  
  # Run Maaslin2
  fit_data <- Maaslin2(input_data = input_data, 
                       input_metadata = df_input_metadata_subset, 
                       output = output_path, 
                       fixed_effects = {{site_fixed_effects}}, 
                       random_effects = c("MouseID"),
                       min_prevalence = 0,
                       reference = c("Sequencing_Run,Jan_2017", "Site,Distal_Colon"),
                       normalization = "clr", 
                       transform = "none",
                       plot_heatmap = FALSE,
                       plot_scatter = FALSE)
  
  return(fit_data)
}

# List of sites and corresponding file paths

site_paths <- list(
  Duodenum = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID",
  Jejunum = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID",
  Ileum = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID",
  Cecum = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID",
  Proximal_Colon = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID",
  Distal_Colon = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID"
)

site_general_paths <- list(
  SI = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-SI-ComBat-SeqRunSexType-1-MsID",
  Colon = "Donors-Analysis/differential_genera_type/L2-LumRef-CLR-Colon-ComBat-SeqRunSexType-1-MsID"
)


# Iterate over each site
site_fe <- c("Sequencing_Run", "Sex", "Type")
for (site in names(site_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_fe)
}

site_general_fe <- c("Sequencing_Run", "Sex", "Site","Type")
for (site in names(site_general_paths)) {
  # Get sample IDs for the current site
  samples <- df_input_metadata %>% 
    filter(Site_General == site, SampleID %in% names(df_input_data)) %>% 
    pull(SampleID)
  
  # Subset input data for the current site
  input_data <- df_input_data[, samples]
  
  # Run Maaslin2 analysis
  output_path <- paste0(site_general_paths[[site]])
  fit_data <- run_maaslin2(input_data, output_path, site, site_general_fe)
}

