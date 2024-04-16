library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(here)

setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Shotgun-Maaslin2-SITE-Genus.R")


### Select Luminal UCLA O SPF ---
input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the UCLA_O_SPF Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "UCLA_O_SPF")
samples <- metadata %>%
  filter(sampleid %in% names(df_input_data)) %>%
  pull(sampleid)

ucla_o_spf_ASV <- df_input_data[, samples]
target <- colnames(ucla_o_spf_ASV)
metadata = metadata[match(target, metadata$sampleid),]
target == metadata$sampleid

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata$Line <- factor(metadata$Line)
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sampleid
metadata <- metadata %>% select(-c("sampleid"))
ucla_o_shotgun_filepath <- here("Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_LineSexSite-1-MsID")
fit_data = Maaslin2(input_data=ucla_o_spf_ASV, 
                    input_metadata=metadata, 
                    output = ucla_o_shotgun_filepath, 
                    fixed_effects = c("Line","Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Line,JJWT','Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Barplot ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

# UCLA O SPF 

result <- generate_interregional_taxa_barplot_shotgun(
                                         path_to_significant_results_tsv = "Shotgun/UCLA_O_SPF/Species_DCvsJej_CLR_LineSexSite-1-MsID/significant_results.tsv",
                                        titlestring="UCLA O. SPF",
                                        colorvector = cols)

df <- result$dataframe
phylum_names <- df$Phylum

select_cols <- paletteer::paletteer_d("basetheme::royal",6)
seecolor::print_color(select_cols)
names(select_cols) <- unique(phylum_names)
phylum_colors <- select_cols
names(select_cols)

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

ucla_o_shotgun_species <- result$plot+
  theme(axis.text.y = element_text(colour = color_mapping))+
  theme(legend.position = "right")
dev.new(width=10,height=10)
ucla_o_shotgun_species


### Select CS SPF ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the CS SPF Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "CS_SPF")
samples <- metadata %>%
  filter(sampleid %in% names(df_input_data)) %>%
  pull(sampleid)

CS_SPF_ASV <- df_input_data[, samples]
target <- colnames(CS_SPF_ASV)
metadata = metadata[match(target, metadata$sampleid),]
target == metadata$sampleid

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sampleid
metadata <- metadata %>% select(-c("sampleid"))
CS_shotgun_filepath <- here("Shotgun/CS_SPF/Species_DCvsJej_CLR_SexSite-1-MsID")
fit_data = Maaslin2(input_data=CS_SPF_ASV, 
                    input_metadata=metadata, 
                    output = CS_shotgun_filepath, 
                    fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2
## Barplots --

# CS SPF 

cs_result <- generate_interregional_taxa_barplot_SITE(path_to_significant_results_tsv = "Shotgun/CS_SPF/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                      titlestring="CS SPF",
                                                      colorvector = cols)

df <- cs_result$dataframe 
phylum_names <- df$Phylum

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

cs_shotgun_species <- cs_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
dev.new(width=10,height=10)
cs_shotgun_species

### Select SPF Gavage ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the SPF Gavage Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "SPF_Gavage")
samples <- metadata %>%
  filter(sampleid %in% names(df_input_data)) %>%
  pull(sampleid)

SPF_Gavage_ASV <- df_input_data[, samples]
target <- colnames(SPF_Gavage_ASV)
metadata = metadata[match(target, metadata$sampleid),]
target == metadata$sampleid

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sampleid
metadata <- metadata %>% select(-c("sampleid"))
SPF_gavage_shotgun_filepath <- here("Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID")
fit_data = Maaslin2(input_data=SPF_Gavage_ASV, 
                    input_metadata=metadata, 
                    output = SPF_gavage_shotgun_filepath, 
                    fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
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
spf_result <- generate_interregional_taxa_barplot_SITE(path_to_significant_results_tsv = "Shotgun/SPF_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                               titlestring="SPF Gavage",
                                                               colorvector = cols)
   
df <- spf_result$dataframe 
phylum_names <- df$Phylum

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

spf_shotgun_species <- spf_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
dev.new(width=10,height=10)
spf_shotgun_species

### Select Hum Gavage ---

input_data <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV.tsv"))
df_input_data <- as.data.frame(input_data)
rownames(df_input_data)<-input_data$Species
df_input_data <- df_input_data %>% select(-c("Species"))

# Ensure samples are listed in the same order 
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim = "\t") #mapping file

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type)
sapply(df_input_metadata,levels)

#Pull out the HUM Gavage Samples 
metadata <- df_input_metadata %>%
  filter(Dataset == "HUM_Gavage")
samples <- metadata %>%
  filter(sampleid %in% names(df_input_data)) %>%
  pull(sampleid)

HUM_Gavage_ASV <- df_input_data[, samples]
target <- colnames(HUM_Gavage_ASV)
metadata = metadata[match(target, metadata$sampleid),]
target == metadata$sampleid

metadata$Site <- factor(metadata$Site, levels=c("Distal_Colon","Jejunum"))
metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$sampleid
metadata <- metadata %>% select(-c("sampleid"))
HUM_Gavage_shotgun_filepath <- here("Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID")
fit_data = Maaslin2(input_data=HUM_Gavage_ASV, 
                    input_metadata=metadata, 
                    output = HUM_Gavage_shotgun_filepath, 
                    fixed_effects = c("Sex", "Site"), random_effects = c("MouseID"),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    reference= c('Site,Distal_Colon'),
                    plot_heatmap = FALSE,plot_scatter = FALSE)
?Maaslin2

## HUM Gavage --
hum_result <- generate_interregional_taxa_barplot_SITE(path_to_significant_results_tsv = "Shotgun/HUM_Gavage/Species_DCvsJej_CLR_SexSite-1-MsID/significant_results.tsv",
                                                                       titlestring="HUM Gavage",
                                                                       colorvector = cols)

df <- hum_result$dataframe 
phylum_names <- df$Phylum

# Create a named vector of colors using the phylum color vector
color_mapping <- phylum_colors[phylum_names]
print(color_mapping)

hum_shotgun_species <- hum_result$plot +
  theme(axis.text.y = element_text(colour = color_mapping))
dev.new(width=10,height=10)
hum_shotgun_species
