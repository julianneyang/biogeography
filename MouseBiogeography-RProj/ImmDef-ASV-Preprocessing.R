###Purpose: Use only features from first dataset present in this dataset, then remove 2017_April and 2015_September Sequencing Run
###then perform CombatSeq2
library(dplyr)
library(sva)
library(dada2)
###for Immune Deficiency Dataset
## Generate a target vector with the features from the main dataset to query against the immune deficiency dataset 

original_features <-readr::read_delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/UCLA-ComBat-Adjusted-ASV.tsv"))
target <- original_features$`#OTU.ID`
target <- gsub('.{1}$', '', target)


immdef_ASV<-readr::read_csv(here("UCLA_V_SPF_Analysis/Immune Deficiency ASV - ImmDef-ASV-table.csv"))
immdef_ASV <- as.data.frame(immdef_ASV)
row.names(immdef_ASV)<-immdef_ASV$feature
immdef_ASV <- immdef_ASV %>% select(-c("#OTU ID","feature","taxonomy"))
summary(colSums(immdef_ASV))
immdef_ASV<- immdef_ASV[colSums(immdef_ASV) >= 10000]

## Retrieve only WT Mucosal samples 
metadata<- readr::read_csv(here("UCLA_V_SPF_Analysis/Full-Metadata.csv"))
metadata$SampleID <- gsub("-",".",metadata$SampleID)
names(immdef_ASV)<-gsub("-",".",names(immdef_ASV))
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata %>% filter(Type =="Mucosal" & Genotype =="WT")
samples <- metadata %>% filter(Type =="Mucosal" & Genotype =="WT", SampleID %in% names(immdef_ASV)) %>% pull(SampleID)

final_immdef_ASV <- immdef_ASV[, samples]
#final_immdef_ASV$`#OTU.ID` <- row.names(final_immdef_ASV)
#final_immdef_ASV <- final_immdef_ASV %>% select("#OTU.ID", everything())

#readr::write_delim(final_immdef_ASV,here("UCLA_V_SPF_Analysis/WTCohort_ASV_for_alpha_diversity.tsv"), delim="\t")


## Query the target vector against immunodeficiency ASV
final_immdef_ASV$feature <- row.names(final_immdef_ASV)

final_immdef_ASV<-final_immdef_ASV[final_immdef_ASV$feature %in% target,]
final_immdef_ASV<-select(final_immdef_ASV,-feature)

final_immdef_ASV$`#OTU.ID` <- row.names(final_immdef_ASV)
final_immdef_ASV <- final_immdef_ASV %>% select("#OTU.ID", everything())

readr::write_delim(final_immdef_ASV,here("UCLA_V_SPF_Analysis/starting_files/UCLA_V_SPF_min10k_ASV.tsv"), delim="\t")

## Generate a taxonomy key ---
sequences<- readr::read_csv(here("UCLA_V_SPF_Analysis/starting_files/UCLA_V_SPF_Original_Taxonomy_Key.csv"))
seqs <- sequences$OTU
set.seed(100)
taxa <- assignTaxonomy(seqs, 
                       "../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
output$ASV <- row.names(output)
taxonomy_for_qza <- output %>% select(c("ASV","taxonomy"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Feature ID" = "ASV"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Taxon" = "taxonomy"))
readr::write_delim(taxonomy_for_qza,here("UCLA_V_SPF_Analysis/starting_files/ucla_v_taxonomy.tsv"),delim="\t")

#write.table(final_immdef_ASV, "WTCohort-ImmDef-ASV.tsv", sep= "\t")
