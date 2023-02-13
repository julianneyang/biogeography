##################Assign Taxonomy after running DeBlur for denoising reads########################
library(dada2)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/All SLC Cre Biogeography FastQ Files/532Itgdel-allSLCCre-table-export")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/All SLC Cre Biogeography FastQ Files/extracted-rep-seqs-2/")
set.seed(100) # Initialize random number generator for reproducibility
#export your rep-seqs.fna file output from deblur denoise-16S
taxa <- assignTaxonomy("dna-sequences.fasta", "C:/Users/Jacobs Laboratory/Desktop/silva_nr_v132_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "C:/Users/Jacobs Laboratory/Desktop/silva_species_assignment_v132.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
write.csv(output, "SLC_Deblur_all_Silva_taxonomy_assignments.csv")

#################Split rep-seqs FASTA file into two columns (so you have the QIIME gibberish in column A and the ASV sequence in column B)########################
library(dplyr)
library(tidyr)
fasta<- read.table("dna-sequences-fasta.txt")
fasta2<-fasta %>% separate(V1, c("ASV","QIIME_sequence"), sep = "([>])")
fasta2[fasta2==""]<-NA

ASV <- fasta2$ASV
fasta3 <- data.frame(ASV)
fasta3 <-as.data.frame(fasta3[-1,])
QIIME_seqs <- data.frame(fasta2$QIIME_sequence)
QIIME_seqs <- data.frame(QIIME_seqs[-5792,])  #modify the number to be the number of the last row
fasta3$QIIME_seqs<- QIIME_seqs$QIIME_seqs..5792... #modify the number to be the number of the last row
finalfasta<-na.omit(fasta3)
write.csv(finalfasta, "SLC_CRE_aligned_dna_sequences.csv")
