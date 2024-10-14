library(devtools)
devtools::install_github("borenstein-lab/mimosa2")
library(mimosa)
test_m2_analysis(test_vsearch = T)
library(here)
here::i_am("MouseBiogeography-RProj/MIMOSA.R")
library(tidyverse)
library(data.table)

library(KEGGREST)
KEGGREST::
# The following takes apart the pipeline written by Borenstein lab for easier handling of predicting metabolites from microbiome data:
# https://rdrr.io/github/borenstein-lab/mimosa2/src/R/core_mimosa2_funcs.R

### Read in KO data
dat <- read.delim(here("Donors-Analysis/picrust_output/export_ko_metagenome/feature-table.tsv"), header=T, sep="\t", row.names=1)
dat <- dat %>% select(-c(taxonomy))
dat <- dat %>% mutate_all(as.numeric)
names(dat)<-gsub("X","",names(dat))

### Download reference database ---
download_reference_data(seq_db = "Sequence variants (ASVs)", target_db = "AGORA genomes and models")

### Map KOs to Kegg ---
simThreshold <- 0.05
ASV_list <- row.names(tax)
vsearch_loc <- "/home/julianne/vsearch-2.29.0-linux-x86_64/bin/vsearch"
?map_seqvar
seq_results <-  map_seqvar(
  seqs = ASV_list, 
  repSeqDir = here("MouseBiogeography-RProj/ASV_AGORA/data/AGORA/"), 
  repSeqFile = "agora_NCBI_16S.udb", 
  method = "vsearch", 
  file_prefix = "Donors", 
  seqID = simThreshold, 
  add_agora_names = TRUE, 
  vsearch_path = vsearch_loc)

### Merge output of map_seqvar with ASV table ---
species <- tax %>%
  mutate(seqID = paste0("seq", row_number()))
samps <- names(species)[!names(species) %in% c("OTU", "seqID")]
new_species <- merge(species, seq_results, by = "seqID", all.x = TRUE)

# Aggregate species abundances by dbID
new_species <- new_species %>%
  group_by(dbID) %>%
  summarize(across(all_of(samps), ~ sum(.x, na.rm = TRUE))) %>%
  ungroup()

# Replace missing dbID with "Other"
new_species <- new_species %>%
  mutate(dbID = if_else(is.na(dbID), "Other", dbID))
new_species <- new_species %>%
  rename(OTU = dbID)

species <- new_species
mod_list <- unique(species$OTU)

### Build predicted metabolic network ---
network = build_species_networks_w_agora(mod_list, agora_path = here("MouseBiogeography-RProj/ASV_AGORA/data/AGORA/RxnNetworks/"))
species_table <- data.table(species)
setnames(network, "Rev", "Reversible")
cmp_results <- get_species_cmp_scores(
  species_table = species_table,
  network = network,
  normalize=TRUE,
  relAbund = TRUE,
  manual_agora=FALSE,
  humann2=FALSE,
  leave_rxns=FALSE,
  kos_only = FALSE,
  remove_rev = TRUE,
  fill_zeros = TRUE
  
)

cmp_results_agora <- get_species_cmp_scores(
  species_table = species_table,
  network = network,
  normalize=TRUE,
  relAbund = TRUE,
  manual_agora=TRUE,
  humann2=FALSE,
  leave_rxns=FALSE,
  kos_only = FALSE,
  remove_rev = TRUE,
  fill_zeros = TRUE
  
)


head(cmp_results,n = 100)
compound_list <- as.character(cmp_results$compound[1:20])
get_compound_format(compound_list)
