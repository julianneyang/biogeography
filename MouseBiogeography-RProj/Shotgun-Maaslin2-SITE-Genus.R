library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(here)
library(Microbiome.Biogeography)

here::i_am("MouseBiogeography-RProj/Shotgun-Maaslin2-SITE-Genus.R")

### Select Luminal UCLA O SPF ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV - Shotgun_metagenomics_biogeography_species_counts.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata - Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the UCLA_O_SPF Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "UCLA_O_SPF")
samples <- metadata %>%
  filter(BioGeo_SampleID %in% names(df_input_data)) %>%
  pull(BioGeo_SampleID)

ucla_o_spf_ASV <- df_input_data[, samples]
target <- colnames(ucla_o_spf_ASV)
metadata = metadata[match(target, metadata$BioGeo_SampleID),]
target == metadata$BioGeo_SampleID

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata$Line <- factor(metadata$Line)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$BioGeo_SampleID
metadata <- metadata %>% select(-c("BioGeo_SampleID"))
ucla_o_shotgun_filepath <- here("Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_SeqRunLineSexSite-1-MsID")
fit_data = Maaslin2(input_data=ucla_o_spf_ASV, 
                    input_metadata=metadata, 
                    output = ucla_o_shotgun_filepath, 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Line,JJWT','Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Heatmap ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

# UCLA O SPF 
lumtarget <- find_concordant_features_across_sites("../Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_SeqRunLineSexSite-1-MsID/significant_results.tsv")
print(lumtarget)
modify_names <- gsub(".*\\.f__", "", lumtarget)

modified_results <- readr::read_delim(here("Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_SeqRunLineSexSite-1-MsID/all_results.tsv"))
modified_results$feature <- gsub(".*\\.f__", "",modified_results$feature)
readr::write_delim(modified_results, "../Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_SeqRunLineSexSite-1-MsID/modified_all_results.tsv",delim = "\t")
modified_results$feature
lumtarget
ucla_o_shotgun_species_heatmap <- generate_taxa_heat_map_by_site("../Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_SeqRunLineSexSite-1-MsID/modified_all_results.tsv",
                                                        modify_names,
                                                        "UCLA O. SPF Luminal",
                                                        cols,
                                                        bk)
dev.new(width=10,height=10)
ucla_o_shotgun_species_heatmap


### Select CS SPF ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV - Shotgun_metagenomics_biogeography_species_counts.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata - Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the CS SPF Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "CS_SPF")
samples <- metadata %>%
  filter(BioGeo_SampleID %in% names(df_input_data)) %>%
  pull(BioGeo_SampleID)

CS_SPF_ASV <- df_input_data[, samples]
target <- colnames(CS_SPF_ASV)
metadata = metadata[match(target, metadata$BioGeo_SampleID),]
target == metadata$BioGeo_SampleID

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$BioGeo_SampleID
metadata <- metadata %>% select(-c("BioGeo_SampleID"))
CS_shotgun_filepath <- here("Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID")
fit_data = Maaslin2(input_data=CS_SPF_ASV, 
                    input_metadata=metadata, 
                    output = CS_shotgun_filepath, 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Heatmap ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

# CS SPF 
lumtarget <- find_concordant_features_across_sites("../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/significant_results.tsv")
print(lumtarget)
modify_names <- gsub(".*\\.f__", "", lumtarget)

modified_results <- readr::read_delim(here("Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/all_results.tsv"))
modified_results$feature <- gsub(".*\\.f__", "",modified_results$feature)
readr::write_delim(modified_results, "../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",delim = "\t")
modified_results$feature
lumtarget
cs_shotgun_species_heatmap <- generate_taxa_heat_map_by_site("../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",
                                                                 modify_names,
                                                                 "CS SPF Luminal",
                                                                 cols,
                                                                 bk)
dev.new(width=10,height=10)
cs_shotgun_species_heatmap

# CS SPF 
lumtarget <- find_concordant_features_across_sites("../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/significant_results.tsv")
print(lumtarget)
modify_names <- gsub(".*\\.f__", "", lumtarget)

modified_results <- readr::read_delim(here("Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/all_results.tsv"))
modified_results$feature <- gsub(".*\\.f__", "",modified_results$feature)
readr::write_delim(modified_results, "../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",delim = "\t")
modified_results$feature
lumtarget
cs_shotgun_species_heatmap <- generate_taxa_heat_map_by_site("../Shotgun/CS_SPF/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",
                                                             modify_names,
                                                             "CS SPF Luminal",
                                                             cols,
                                                             bk)
dev.new(width=10,height=10)
cs_shotgun_species_heatmap

### Select SPF Gavage ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV - Shotgun_metagenomics_biogeography_species_counts.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata - Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the SPF Gavage Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "SPF_Gavage")
samples <- metadata %>%
  filter(BioGeo_SampleID %in% names(df_input_data)) %>%
  pull(BioGeo_SampleID)

SPF_Gavage_ASV <- df_input_data[, samples]
target <- colnames(SPF_Gavage_ASV)
metadata = metadata[match(target, metadata$BioGeo_SampleID),]
target == metadata$BioGeo_SampleID

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$BioGeo_SampleID
metadata <- metadata %>% select(-c("BioGeo_SampleID"))
SPF_gavage_shotgun_filepath <- here("Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID")
fit_data = Maaslin2(input_data=SPF_Gavage_ASV, 
                    input_metadata=metadata, 
                    output = SPF_gavage_shotgun_filepath, 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Heatmap ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

# SPF Gavage
lumtarget <- find_concordant_features_across_sites("../Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/significant_results.tsv")
print(lumtarget)
modify_names <- gsub(".*\\.f__", "", lumtarget)

modified_results <- readr::read_delim(here("Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/all_results.tsv"))
modified_results$feature <- gsub(".*\\.f__", "",modified_results$feature)
readr::write_delim(modified_results, "../Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",delim = "\t")
modified_results$feature
lumtarget
spf_gavage_shotgun_species_heatmap <- generate_taxa_heat_map_by_site("../Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",
                                                             modify_names,
                                                             "SPF Gavage Luminal",
                                                             cols,
                                                             bk)
dev.new(width=10,height=10)
spf_gavage_shotgun_species_heatmap

### Select Hum Gavage ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV - Shotgun_metagenomics_biogeography_species_counts.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata - Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the SPF Gavage Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "HUM_Gavage")
samples <- metadata %>%
  filter(BioGeo_SampleID %in% names(df_input_data)) %>%
  pull(BioGeo_SampleID)

HUM_Gavage_ASV <- df_input_data[, samples]
target <- colnames(HUM_Gavage_ASV)
metadata = metadata[match(target, metadata$BioGeo_SampleID),]
target == metadata$BioGeo_SampleID

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$BioGeo_SampleID
metadata <- metadata %>% select(-c("BioGeo_SampleID"))
HUM_Gavage_shotgun_filepath <- here("Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID")
fit_data = Maaslin2(input_data=HUM_Gavage_ASV, 
                    input_metadata=metadata, 
                    output = HUM_Gavage_shotgun_filepath, 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Heatmap ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

# HUM Gavage
lumtarget <- find_concordant_features_across_sites("../Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/significant_results.tsv")
print(lumtarget)
modify_names <- gsub(".*\\.f__", "", lumtarget)

modified_results <- readr::read_delim(here("Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/all_results.tsv"))
modified_results$feature <- gsub(".*\\.f__", "",modified_results$feature)
readr::write_delim(modified_results, "../Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",delim = "\t")
modified_results$feature
lumtarget
HUM_Gavage_shotgun_species_heatmap <- generate_taxa_heat_map_by_site("../Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SeqRunSexSite-1-MsID/modified_all_results.tsv",
                                                                     modify_names,
                                                                     "HUM Gavage Luminal",
                                                                     cols,
                                                                     bk)
dev.new(width=10,height=10)
HUM_Gavage_shotgun_species_heatmap
