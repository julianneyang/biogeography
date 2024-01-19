
library(devtools)
#install_github("omixer/omixer-rpmR")
library(omixerRpm)
library(here)
library(dplyr)

here::i_am("MouseBiogeography-RProj/omixer-RPM.R")
#load the database (select between gut-brain and gut-metabolic)
listDB()


# Wrangle the samples so it's only what was used for metaphlan
shotgun_dat <- read.table("../Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv", header=T, sep="\t", row.names=1)
shotgun_dat <- shotgun_dat %>% mutate_all(as.numeric)
names(shotgun_dat) <- gsub("_merged__kneaddata_paired_Abundance.RPKs","",names(shotgun_dat))
names(shotgun_dat) <- gsub("merged_","",names(shotgun_dat))
names(shotgun_dat) <- gsub("__kneaddata_paired_Abundance.RPKs","",names(shotgun_dat))
summary(colSums(shotgun_dat))

meta <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t")
meta$sampleid
meta$humann_sampleid <- gsub("_merged_R1_001.fastq_profiled_metagenome","",meta$sampleid)
meta$humann_sampleid <- gsub("_R1_001.fastq_profiled_metagenome","",meta$humann_sampleid)
meta$humann_sampleid <- gsub("-",".",meta$humann_sampleid)
readr::write_csv(meta, here("Shotgun/BioGeo_Shotgun_Metadata.tsv"))

samples <- meta %>%
  filter(humann_sampleid %in% names(shotgun_dat)) %>%
  pull(humann_sampleid)

shotgun_dat <- shotgun_dat[, samples]
shotgun_dat$KO <- row.names(shotgun_dat)
shotgun_dat <-shotgun_dat %>% select(KO,everything())


#Utilize GOMIZER to find optimal module coverage value
db <- loadDB("GBMs.v1.0")
mods <- rpm(shotgun_dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present

db<-loadDefaultDB() # by default is metabolic
mods <- rpm(shotgun_dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present


# get the abundance|coverage as a data.frame with module id and description
coverage <- asDataFrame(mods, "coverage")


