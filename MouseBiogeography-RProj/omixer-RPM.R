
library(devtools)
#install_github("omixer/omixer-rpmR")
library(omixerRpm)
library(here)
library(dplyr)

here::i_am("MouseBiogeography-RProj/omixer-RPM.R")
#load the database (select between gut-brain and gut-metabolic)
listDB()



# Wrangle the samples so it's only what was used for metaphlan
dat <- read.table(here("Donors-Analysis/picrust_output/export_ko_metagenome/feature-table.tsv"), header=T, sep="\t", row.names=1)
dat <- dat %>% select(-c(taxonomy))
dat <- dat %>% mutate_all(as.numeric)
dat$KO <- row.names(dat)
dat <-dat %>% select(KO,everything())


#Utilize GOMIZER to find optimal module coverage value
db <- loadDB("GBMs.v1.0")
mods <- rpm(dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present
db<-loadDefaultDB() # by default is metabolic
mods <- rpm(dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present

?is.character()
!is.character(dat)

# get the abundance|coverage as a data.frame with module id and description
coverage <- asDataFrame(mods, "coverage")


