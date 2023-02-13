###Purpose: Use only features from first dataset present in this dataset, then remove 2017_April and 2015_September Sequencing Run
###then perform CombatSeq2
library(dplyr)
library(sva)
###for Immune Deficiency Dataset
## Generate a target vector with the features from the main dataset to query against the immune deficiency dataset 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/ImmDef-Mouse-Biogeography-Analysis/")
original_features <-read.csv("Immune Deficiency ASV - target.csv", header=TRUE, row.names=1)
target <- original_features$feature


immdef_ASV<-read.csv("Immune Deficiency ASV - ImmDef-ASV-table.csv", header=TRUE)

## Retrieve only WT Mucosal samples 
metadata<- read.csv("Full-Metadata.csv")
metadata$SampleID <- gsub("-",".",metadata$SampleID)
names(immdef_ASV)
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata %>% filter(Type =="Mucosal" & Genotype =="WT")
samples <- metadata %>% filter(Type =="Mucosal" & Genotype =="WT", SampleID %in% names(immdef_ASV)) %>% pull(SampleID)

final_immdef_ASV <- immdef_ASV[, samples]
row.names(final_immdef_ASV) <- immdef_ASV$feature

write.table(final_immdef_ASV,"WTCohort_ASV_for_alpha_diversity.tsv", sep="\t")

final_immdef_ASV$feature<-immdef_ASV$feature

## Query the target vector against immunodeficiency ASV
final_immdef_ASV<-final_immdef_ASV[final_immdef_ASV$feature %in% target,]
final_immdef_ASV<-select(final_immdef_ASV,-feature)

## Perform CombatSeq2 by Type and Site
newmeta <- metadata %>% filter(Type =="Mucosal" & Genotype =="WT")

batch<- as.character(newmeta$Sequencing_Run)
modcombat=model.matrix(~ Sex + Site,data=newmeta)

input=as.matrix(final_immdef_ASV)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  

write.table(final_immdef_ASV, "WTCohort-ImmDef-ASV.tsv", sep= "\t")
write.table(combat_adjusted_counts,"WTCohort-ImmDef-ComBat-Adjusted-ASV.tsv", sep="\t") 
