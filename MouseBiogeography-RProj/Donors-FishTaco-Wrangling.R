library(here)
library(dplyr)
library(tidyverse)
library(readr)

here::i_am("MouseBiogeography-RProj/Donors-FishTaco-Wrangling.R")

# Wrangle PICRUSt KO data into relabun values 
dat <- read.delim(here("Donors-Analysis/picrust_output/export_ko_metagenome/feature-table.tsv"), header=T, sep="\t", row.names=1)
dat <- dat %>% select(-c(taxonomy))
dat <- dat %>% mutate_all(as.numeric)
names(dat)<-gsub("X","",names(dat))

column_sums <- colSums(dat)
relative_abundances <- sweep(dat, 2, column_sums, FUN = "/")
colSums(relative_abundances)

# Wrangle Taxon data into relabun values 
tax <- read.delim(here("Donors-Analysis/starting_files/Donors-Mice-1xPrev0.15-ComBat-ASV.tsv"), header=T, sep="\t", row.names=1)
tax <- tax %>% mutate_all(as.numeric)
names(tax)<-gsub("X","",names(tax))

column_sums <- colSums(tax)
tax_rel_abun <- sweep(tax, 2, column_sums, FUN = "/")
colSums(tax_rel_abun)

# Subset KO and Taxon data to retain only pairwise comparison samples
input_metadata <-read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

# Sanity check that names match 
target <- colnames(dat)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

# Save output files for fishtaco
duo_dc <- input_metadata %>% 
  filter(Site =="Duodenum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% pull(SampleID)
duo_dc_ko <- relative_abundances[, duo_dc]
duo_dc_ko <- duo_dc_ko %>% 
  rownames_to_column(var="Function") 
readr::write_delim(duo_dc_ko, here("Donors-Analysis/fish_taco/duo_dc_ko_relab.tsv"),delim="\t")

duo_dc_tax <- tax_rel_abun[, duo_dc]
duo_dc_tax <- duo_dc_tax %>% 
  rownames_to_column(var="Taxa")
write.table(duo_dc_tax, here("Donors-Analysis/fish_taco/duo_dc_ko_tax.tsv"))

duo_dc_labels <- input_metadata %>% 
  filter(Site =="Duodenum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Duodenum", 1, 0)) %>% 
  select("Label")
duo_dc_labels <- duo_dc_labels %>% 
  rownames_to_column(var = "Sample")
write.table(duo_dc_labels, here("Donors-Analysis/fish_taco/duo_dc_labels.tsv"))
