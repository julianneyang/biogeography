
library(devtools)
#install_github("omixer/omixer-rpmR")
library(omixerRpm)
library(here)
library(dplyr)

here::i_am("MouseBiogeography-RProj/omixer-RPM.R")
#load the database (select between gut-brain and gut-metabolic)
listDB()
db<-loadDefaultDB() # by default is metabolic
db <- loadDB("GBMs.v1.0")
#read table in special format 
dat <- read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/OMIXER_Full_KO_Counts - KO_all_counts.tsv", header=T, sep="\t"), row.names=1
dat <- select(dat, -id)
?rpm()

#Utilize GOMIZER to find optimal module coverage value
mods <- rpm(dat, minimum.coverage=0.6, annotation = 1, module.db=db) #coverage- what percent of KO need to be present in one module for one mod to be present
?is.character()
!is.character(dat)


