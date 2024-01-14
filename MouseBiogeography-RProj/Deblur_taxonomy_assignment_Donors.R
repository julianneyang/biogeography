### Assign Taxonomy after running DeBlur for denoising reads: Donors ---
library(dada2)
library(here)

BiocManager::install("dada2")

set.seed(100) # Initialize random number generator for reproducibility

## export your rep-seqs.fna file output from deblur denoise-16S and use it to assign taxonomy --
# for paired end reads

taxa <- assignTaxonomy("../Donors-Analysis/export_Donors-paired-end-rep-seqs/dna-sequences.fasta", "../../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
write.csv(output, "../Donors-Analysis/Donors_PE_taxonomy_assignments.csv")

## Split rep-seqs FASTA file into two columns (so you have the QIIME gibberish in column A and the ASV sequence in column B) ---
# for paired end reads 
library(dplyr)
library(tidyr)
fasta<- read.table("../Donors-Analysis/export_Donors-paired-end-rep-seqs/dna-sequences.fasta")
fasta2<-fasta %>% 
  separate(V1, c("ASV","QIIME_sequence"), sep = "([>])")
fasta2[fasta2==""]<-NA

ASV <- fasta2$ASV
fasta3 <- data.frame(ASV)
fasta3 <-as.data.frame(fasta3[-1,])
QIIME_seqs <- data.frame(fasta2$QIIME_sequence)
QIIME_seqs <- data.frame(QIIME_seqs[-10656,])  #modify the number to be the number of the last row
fasta3$QIIME_seqs<- QIIME_seqs$QIIME_seqs..10656... #modify the number to be the number of the last row
finalfasta<-na.omit(fasta3)
write.csv(finalfasta, "../Donors-Analysis/Donors_PE_aligned_DNA_sequences.csv")


## export your rep-seqs.fna file output from deblur denoise-16S and use it to assign taxonomy --
# for single end reads

taxa <- assignTaxonomy("../Donors-Analysis/export_Donors-single-filtered-rep-seqs/dna-sequences.fasta", "../../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
write.csv(output, "../Donors-Analysis/Donors_single_taxonomy_assignments.csv")

## Split rep-seqs FASTA file into two columns (so you have the QIIME gibberish in column A and the ASV sequence in column B) ---
# for single end reads 
library(dplyr)
library(tidyr)
fasta<- read.table("../Donors-Analysis/export_Donors-single-filtered-rep-seqs/dna-sequences.fasta")
fasta2<-fasta %>% 
  separate(V1, c("ASV","QIIME_sequence"), sep = "([>])")
fasta2[fasta2==""]<-NA

ASV <- fasta2$ASV
fasta3 <- data.frame(ASV)
fasta3 <-as.data.frame(fasta3[-1,])
QIIME_seqs <- data.frame(fasta2$QIIME_sequence)
QIIME_seqs <- data.frame(QIIME_seqs[-1528,])  #modify the number to be the number of the last row
fasta3$QIIME_seqs<- QIIME_seqs$QIIME_seqs..1528... #modify the number to be the number of the last row
finalfasta<-na.omit(fasta3)
write.csv(finalfasta, "../Donors-Analysis/Donors_single_aligned_DNA_sequences.csv")

## Combine taxonomy assignments --
PE_taxonomy <- read.csv("")
