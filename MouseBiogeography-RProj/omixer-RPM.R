
library(devtools)
#install_github("omixer/omixer-rpmR")
library(omixerRpm)
library(here)
library(dplyr)

here::i_am("MouseBiogeography-RProj/omixer-RPM.R")
#load the database (select between gut-brain and gut-metabolic)
listDB()

db <- loadDB("GBMs.v1.0")
db<-loadDefaultDB() # by default is metabolic

# Wrangle the samples so it's only what was used for metaphlan
shotgun_dat <- read.table("../Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv", header=T, sep="\t", row.names=1)
shotgun_dat <- shotgun_dat %>% mutate_all(as.numeric)
names(shotgun_dat) <- gsub("_merged__kneaddata_paired_Abundance.RPKs","",names(shotgun_dat))
names(shotgun_dat) <- gsub("merged_","",names(shotgun_dat))
names(shotgun_dat) <- gsub("__kneaddata_paired_Abundance.RPKs","",names(shotgun_dat))


meta <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t")
meta$sampleid
meta$humann_sampleid <- gsub("_merged_R1_001.fastq_profiled_metagenome","",meta$sampleid)
meta$humann_sampleid <- gsub("_R1_001.fastq_profiled_metagenome","",meta$humann_sampleid)
meta$humann_sampleid <- gsub("-",".",meta$humann_sampleid)

shotgun_dat$KO <- row.names(shotgun_dat)
shotgun_dat <-shotgun_dat %>% select(KO,everything())


dat <- read.table("../Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/OMIXER_Full_KO_Counts - KO_all_counts.tsv", header=T, sep="\t", row.names=1)
db
?rpm()

#Utilize GOMIZER to find optimal module coverage value
mods <- rpm(shotgun_dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present

?is.character()
!is.character(dat)

# get the abundance|coverage as a data.frame with module id and description
coverage <- asDataFrame(mods, "coverage")


