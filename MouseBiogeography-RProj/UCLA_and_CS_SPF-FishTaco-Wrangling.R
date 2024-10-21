library(here)
library(dplyr)
library(tidyverse)
library(readr)

here::i_am("MouseBiogeography-RProj/UCLA_and_CS_SPF-FishTaco-Wrangling.R")

# Wrangle PICRUSt KO data into relabun values 
ucla_o_spf <- read.delim(here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_relabun.tsv"), header=T, sep="\t", row.names=1)
cs_spf <- read.delim(here("CS_SPF/picrust_output/export_ko_metagenome/feature-table.tsv"), header=T, sep="\t", row.names=1)

cs_spf <- cs_spf %>% select(-c(taxonomy))
cs_spf <- cs_spf %>% mutate_all(as.numeric)
names(cs_spf)<-gsub("X","",names(cs_spf))

column_sums <- colSums(cs_spf)
cs_spf <- sweep(cs_spf, 2, column_sums, FUN = "/")
colSums(relative_abundances)



# Wrangle Taxon data into relabun values 
tax_ucla_o_spf <- read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/UCLA-ComBat-Adjusted-ASV.tsv"), header=T, sep="\t", row.names=1)
tax_ucla_o_spf <- tax_ucla_o_spf %>% mutate_all(as.numeric)
names(tax_ucla_o_spf)<-gsub("X","",names(tax_ucla_o_spf))

column_sums <- colSums(tax_ucla_o_spf)
tax_ucla_o_spf_rel_abun <- sweep(tax_ucla_o_spf, 2, column_sums, FUN = "/")
colSums(tax_ucla_o_spf_rel_abun)

tax_cs_spf <- read.delim(here("CS_SPF/starting_files/CS-Facility-ComBat-Adjusted-ASV.tsv"), header=T, sep="\t", row.names=1)
tax_cs_spf <- tax_cs_spf %>% mutate_all(as.numeric)
names(tax_cs_spf)<-gsub("X","",names(tax_cs_spf))

column_sums <- colSums(tax_cs_spf)
tax_cs_spf_rel_abun <- sweep(tax_cs_spf, 2, column_sums, FUN = "/")
colSums(tax_cs_spf_rel_abun)

tax_ucla_o_spf <- tax_ucla_o_spf %>% rownames_to_column("ASV")
tax_cs_spf <- tax_cs_spf %>% rownames_to_column("ASV")



# Subset KO and Taxon data to retain only pairwise comparison samples
input_metadata <-read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

target <- colnames(ucla_o_spf)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

# Save output files for fishtaco
duo_dc <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Duodenum" | Site == "Distal_Colon", SampleID %in% names(ucla_o_spf)) %>% pull(SampleID)
duo_dc_ko <- ucla_o_spf[, duo_dc]
duo_dc_ko <- duo_dc_ko %>%  
  rownames_to_column(var="Function")  
readr::write_delim(duo_dc_ko, here("Regional-Mouse-Biogeography-Analysis/fish_taco/muc_duo_dc_ko.tsv"),delim="\t")

jej_dc <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Jejunum" | Site == "Distal_Colon", SampleID %in% names(ucla_o_spf)) %>% pull(SampleID)
jej_dc_ko <- ucla_o_spf[, jej_dc]
jej_dc_ko <- jej_dc_ko %>%  
  rownames_to_column(var="Function")  
readr::write_delim(jej_dc_ko, here("Regional-Mouse-Biogeography-Analysis/fish_taco/muc_jej_dc_ko.tsv"),delim="\t")

ile_dc <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Ileum" | Site == "Distal_Colon", SampleID %in% names(ucla_o_spf)) %>% pull(SampleID)
ile_dc_ko <- ucla_o_spf[, ile_dc]
ile_dc_ko <- ile_dc_ko %>%  
  rownames_to_column(var="Function")  
readr::write_delim(ile_dc_ko, here("Regional-Mouse-Biogeography-Analysis/fish_taco/muc_ile_dc_ko.tsv"),delim="\t")

#min_samples <- ceiling(ncol(duo_dc_ko) * 0.10)
#duo_dc_ko <- duo_dc_ko %>%  
#filter(rowSums(. > 0) >= min_samples) %>%
#rownames_to_column(var="Function")  


duo_dc_tax <- tax_ucla_o_spf_rel_abun[, duo_dc]
duo_dc_tax <- duo_dc_tax %>% 
  rownames_to_column(var="Taxa")
readr::write_delim(duo_dc_tax, here("Regional-Mouse-Biogeography-Analysis/fish_taco/duo_dc_tax.tsv"),delim="\t")

jej_dc_tax <- tax_ucla_o_spf_rel_abun[, jej_dc]
jej_dc_tax <- jej_dc_tax %>% 
  rownames_to_column(var="Taxa")
readr::write_delim(jej_dc_tax, here("Regional-Mouse-Biogeography-Analysis/fish_taco/jej_dc_tax.tsv"),delim="\t")

ile_dc_tax <- tax_ucla_o_spf_rel_abun[, ile_dc]
ile_dc_tax <- ile_dc_tax %>% 
  rownames_to_column(var="Taxa")
readr::write_delim(ile_dc_tax, here("Regional-Mouse-Biogeography-Analysis/fish_taco/ile_dc_tax.tsv"),delim="\t")

names(duo_dc_tax)==names(duo_dc_ko)
duo_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Duodenum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Duodenum", 0, 1)) %>% 
  select("Label")
duo_dc_labels <- duo_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(duo_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/duo_dc_0_1_labels.tsv"),delim="\t")

duo_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Duodenum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Duodenum", 1, 0)) %>% 
  select("Label")
duo_dc_labels <- duo_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(duo_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/duo_dc_1_0_labels.tsv"),delim="\t")

jej_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Jejunum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Jejunum", 0, 1)) %>% 
  select("Label")
jej_dc_labels <- jej_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(jej_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/jej_dc_0_1_labels.tsv"),delim="\t")

jej_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Jejunum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Jejunum", 0, 1)) %>% 
  select("Label")
jej_dc_labels <- jej_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(jej_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/jej_dc_1_0_labels.tsv"),delim="\t")

ile_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Ileum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Ileum", 0, 1)) %>% 
  select("Label")
ile_dc_labels <- ile_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(ile_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/ile_dc_0_1_labels.tsv"),delim="\t")

ile_dc_labels <- input_metadata %>% 
  filter(Type =="Mucosal") %>%
  filter(Site =="Ileum" | Site == "Distal_Colon", SampleID %in% names(dat)) %>% 
  mutate(Label = ifelse(Site == "Ileum", 1, 0)) %>% 
  select("Label")
ile_dc_labels <- ile_dc_labels %>% 
  rownames_to_column(var = "Sample")
readr::write_delim(ile_dc_labels, here("Regional-Mouse-Biogeography-Analysis/fish_taco/ile_dc_1_0_labels.tsv"),delim="\t")

