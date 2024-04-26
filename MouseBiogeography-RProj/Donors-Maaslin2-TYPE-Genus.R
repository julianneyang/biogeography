library(Maaslin2)
library(funrar)
library(dplyr)
library(here)

here::i_am("MouseBiogeography-RProj/Donors-Maaslin2-TYPE-Genus.R")

input_data <- read.delim(here("Donors-Analysis/starting_files/collapsed_taxa/export_L6_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), header=TRUE,row.names=1) 
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


## Get Site datasets to test Lum v Muc ---
si_samples <- df_input_metadata %>% filter(Site_General =="SI", SampleID %in% names(df_input_data)) %>% pull(SampleID)
si_input_data <- df_input_data[, si_samples]

col_samples <- df_input_metadata %>% filter(Site_General =="Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
col_input_data <- df_input_data[, col_samples]

duo_samples <- df_input_metadata %>% filter(Site =="Duodenum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
duo_input_data <- df_input_data[, duo_samples]

jej_samples <- df_input_metadata %>% filter(Site =="Jejunum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
jej_input_data <- df_input_data[, jej_samples]

ile_samples <- df_input_metadata %>% filter(Site =="Ileum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
ile_input_data <- df_input_data[, ile_samples]

cec_samples <- df_input_metadata %>% filter(Site =="Cecum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
cec_input_data <- df_input_data[, cec_samples]

PC_samples <- df_input_metadata %>% filter(Site =="Proximal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
PC_input_data <- df_input_data[, PC_samples]

DC_samples <- df_input_metadata %>% filter(Site =="Distal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
DC_input_data <- df_input_data[, DC_samples]

## Run Maaslin2 ---
filepath <- "Donors-Analysis/differential_genera_type/"

#Distal_Colon 

target <- colnames(DC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=DC_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex","Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017","Type,Luminal"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)



#Proximal_Colon -
target <- colnames(PC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=PC_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex","Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
target <- colnames(cec_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=cec_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
target <- colnames(ile_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=ile_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum 
target <- colnames(jej_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=jej_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
target <- colnames(duo_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=duo_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon - model fails to converge
target <- colnames(col_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=col_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-Colon-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)


#SI
target <- colnames(si_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=si_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"L6-LumRef-CLR-SI-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0,
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

perform_maaslin_analysis <- function(tissue_type, file_path) {
  # Construct file path for the input data
  input_data_path <- here::here(file_path, paste0("export_L6_", tissue_type, "_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"))
  
  # Read input data
  input_data <- read.delim(input_data_path, header=TRUE, row.names=1)
  df_input_data <- as.data.frame(input_data)
  df_input_data <- dplyr::select(df_input_data, -c("taxonomy"))
  names(df_input_data) <- gsub("X", "", names(df_input_data))
  
  # Read input metadata
  input_metadata <- read.delim(here::here("Donors-Analysis/starting_files/Donors_Metadata.tsv"), header=TRUE)
  input_metadata$SampleID <- gsub("-", ".", input_metadata$SampleID)
  row.names(input_metadata) <- input_metadata$SampleID
  
  # Match metadata with target
  target <- colnames(df_input_data)
  input_metadata <- input_metadata[match(target, row.names(input_metadata)), ]
  
  # Prepare metadata
  df_input_metadata <- input_metadata
  df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
  df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
  df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
  df_input_metadata$Sex <- factor(df_input_metadata$Sex)
  df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
  
  # Perform Maaslin2 analysis
  fit_data <- Maaslin2(
    input_data = df_input_data,
    input_metadata = df_input_metadata,
    output = here::here(paste0("Donors-Analysis/differential_genera_type/L6-LumRef-CLR-", tissue_type, "-ComBat-SeqRunSexType-1-MsID-DonorID")),
    fixed_effects = c("Sequencing_Run", "Sex", "Type"),
    random_effects = c("MouseID", "DonorID"),
    normalization = "clr",
    transform = "none",
    plot_heatmap = FALSE,
    plot_scatter = FALSE,
    min_prevalence = 0.15,
    reference = c("Sequencing_Run,Jan_2017", "Type,Luminal")
  )
}

tissue_types <- c("Cecum", "Proximal_Colon", "Distal_Colon", "Duodenum", "Jejunum", "Ileum","Colon","SI")
filepath <- "Donors-Analysis/type_subsets/collapsed_taxa/genus_level/"


for (tissue_type in tissue_types) {
  perform_maaslin_analysis(tissue_type, filepath)
}




