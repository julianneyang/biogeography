library(melonnpan)
library(here)
library(stringr)
library(Maaslin2)

here::i_am("MouseBiogeography-RProj/Melonnpan_Predict.R")

### Custom Function ---
predict_metabolites <- function(df_input, weights_path, output_dir, train_metag, threshold = 0.01/100, sample_fraction = 0.10) {
  
  # Check if df_input is a file path or a dataframe
  if (is.character(df_input)) {
    # If it's a file path, read the file
    df <- read.delim(df_input, row.names = 1)
  } else if (is.data.frame(df_input)) {
    # If it's already a dataframe, use it directly
    df <- df_input
  } else {
    stop("df_input must be either a file path or a dataframe")
  }
  # Read in the normalized data and transpose it
  shotgun_dat <- as.data.frame(t(df))
  shotgun_dat <- shotgun_dat[,-1]
  
  # Read in the trained weights
  weights <- read.delim(weights_path)
  
  # Find overlapping IDs
  overlap <- intersect(names(shotgun_dat), weights$ID)
  
  # Filter data based on abundance threshold and prevalence across samples
  binary_abundance <- shotgun_dat > threshold
  proportion_samples <- colMeans(binary_abundance)
  filtered_data <- shotgun_dat[, proportion_samples >= sample_fraction]
  
  # Predict metabolites using MelonnPan
  metabolites <- melonnpan::melonnpan.predict(metag = filtered_data,
                                              output = output_dir,
                                              weight.matrix = weights_path,
                                              train.metag = train_metag)
  
  return(metabolites)
}

### Training Data 
microbiome_train <- readRDS("../melonnpan_data/microbiome_training_data.RDS")
?melonnpan.predict()

### Shotgun Data ---
# Define a reusable function for processing each dataset
process_dataset <- function(dataset_name, ko, input_metadata, output_dir) {
  # Filter metadata for the specified dataset
  metadata <- input_metadata %>%
    filter(Dataset == dataset_name)
  
  # Select samples that exist in both metadata and KO data
  samples <- metadata %>%
    filter(humann_sampleid %in% names(ko)) %>%
    pull(humann_sampleid)
  
  # Subset the KO data for the selected samples
  dataset_ASV <- ko[, samples]
  
  # Match metadata with KO data
  metadata <- metadata[match(colnames(dataset_ASV), metadata$humann_sampleid),]
  
  # Predict metabolites using the MelonnPan model
  predict_metabolites(
    df_input = dataset_ASV,
    weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
    output_dir = here(output_dir),
    train_metag = microbiome_train
  )
}

# Load KO data and metadata
ko <- read.delim(here("Shotgun/relab_normalized/ko_final.tsv"), row.names=1) %>%
  select(-KO)  # Drop the KO column

input_metadata <- readr::read_delim(here("Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"), delim = "\t")

# Convert specific columns in metadata to factors
df_input_metadata <- input_metadata %>%
  mutate(across(c(MouseID, Sex, Type), factor))

# Process the UCLA_O_SPF dataset
process_dataset(
  dataset_name = "UCLA_O_SPF",
  ko = ko,
  input_metadata = df_input_metadata,
  output_dir = "Shotgun/melonnpan/UCLA_O_SPF/"
)

# Process the CS_SPF dataset
process_dataset(
  dataset_name = "CS_SPF",
  ko = ko,
  input_metadata = df_input_metadata,
  output_dir = "Shotgun/melonnpan/CS_SPF/"
)

# Process the SPF_Gavage dataset
process_dataset(
  dataset_name = "SPF_Gavage",
  ko = ko,
  input_metadata = df_input_metadata,
  output_dir = "Shotgun/melonnpan/SPF_Gavage/"
)

# Process the HUM_SD_Gavage dataset
process_dataset(
  dataset_name = "HUM_Gavage",
  ko = ko,
  input_metadata = df_input_metadata,
  output_dir = "Shotgun/melonnpan/HUM_SD_Gavage/"
)

### UCLA O. SPF ---
microbiome <- read.delim(here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_counts.tsv"), row.names=1)
microbiome <- microbiome %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))

input_metadata <-read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv"),header=TRUE, row.names=1) #mapping file

target <- colnames(microbiome)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(microbiome)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- microbiome[,samples]

#df_input_data <- filter_features(df_input_data) #2102

# Predict Metabolite Compostion - 
ucla_o_spf_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Regional-Mouse-Biogeography-Analysis/melonnpan/"),
  train_metag = microbiome_train
)

### CS SPF ---
ko <- read.delim(here("CS_SPF/picrust_output/export_ko_metagenome/feature-table.tsv"), row.names=1)
ko <- ko %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
ko <- ko %>% select(-c("taxonomy"))

input_metadata <-read.delim(here("CS_SPF/starting_files/CS_Facility_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
target <- colnames(ko)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(ko)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- ko[,samples]

#df_input_data <- filter_features(df_input_data) #2326

# Predict Metabolite Compostion after feature filtering - 
cs_spf_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
  output_dir = here("CS_SPF/melonnpan/"),
  train_metag = microbiome_train
)

### HUM Gavage --- 
ko <- read.delim(here("Humanized-Biogeography-Analysis/picrust_output/picrust2_output_Humanized_min10k_Humanized_Combat_Adjusted_ASV.qza/export_ko_metagenome/feature-table.tsv"), row.names=1)
ko <- ko %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
ko <- ko %>% select(-c("taxonomy"))

input_metadata <-read.delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
target <- colnames(ko)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(ko)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- ko[,samples]

#df_input_data <- filter_features(df_input_data) #2138

# Predict Metabolite Compostion after feature filtering - 
hum_sd_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Humanized-Biogeography-Analysis/melonnpan/HUM_SD_Gavage/"),
  train_metag = microbiome_train
)

### SPF Gavage ---
ko <- read.delim(here("Humanized-Biogeography-Analysis/picrust_output/picrust2_output_Cedars_SPF_min10k_Humanized_Combat_Adjusted_ASV.qza/export_ko_metagenome/feature-table.tsv"), row.names=1)
ko <- ko %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
ko <- ko %>% select(-c("taxonomy"))

input_metadata <-read.delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
target <- colnames(ko)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(ko)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- ko[,samples] 

#df_input_data <- filter_features(df_input_data)#2122 features 

# Predict Metabolite Compostion after feature filtering - 
spf_gavage_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Humanized-Biogeography-Analysis/melonnpan/SPF_Gavage/"),
  train_metag = microbiome_train
)

spf_gavage_result$RTSI

### HUM MD Gavage ---
ko <- read.delim(here("Donors-Analysis/picrust_output/export_ko_metagenome/feature-table.tsv"), row.names=1)
ko <- ko %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
ko <- ko %>% select(-c("taxonomy"))
names(ko) <- gsub("X","",names(ko))

input_metadata <-read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata) <- gsub("-",".",row.names(input_metadata))
target <- colnames(ko)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(ko)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- ko[,samples]

#df_input_data <- filter_features(df_input_data) #3620

# Predict Metabolite Compostion after feature filtering - 
hum_md_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("melonnpan_model/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Donors-Analysis/melonnpan/"),
  train_metag = microbiome_train
)

hum_md_result$RTSI

