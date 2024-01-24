##################Assign Taxonomy after running DeBlur for denoising reads########################
library(dada2)

set.seed(100) # Initialize random number generator for reproducibility
#export your rep-seqs.fna file output from deblur denoise-16S

here::i_am("../MouseBiogeography-RProj/Deblur_taxonomy_assignment_UCLA.R")
taxa <- assignTaxonomy("../Regional-Mouse-Biogeography-Analysis/repseqs-Regional-Combat-ASV.tsv", 
                       "../../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
output$ASV <- row.names(output)
readr::write_csv(output, here("Regional-Mouse-Biogeography-Analysis/UCLA_taxonomy_assignments.csv"))

