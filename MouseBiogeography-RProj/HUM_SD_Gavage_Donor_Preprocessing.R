BiocManager::install("dada2")
library(dada2)

## Fix mismatch read IDS in hum_SD_Gavage_Donor --
path<- "../HUM_SD_Gavage_Donor"
fnFs <- sort(list.files(path, pattern="R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="R2.fastq", full.names = TRUE))
file_in <- c(fnFs, fnRs)
sample.names <- sapply(strsplit(basename(fnFs), "_R1"), `[`, 1)

outF <- paste0(path,"/",sample.names,"_matched_R1.fastq")
outR <- paste0(path,"/",sample.names,"_matched_R2.fastq")
file_out <- c(outF,outR)
fastqPairedFilter(fn = file_in, fout = file_out, matchIDs=TRUE, compress=FALSE)

## Run Deblur from command line --

## Then, assign taxonomy here --
set.seed(100) # Initialize random number generator for reproducibility
#export your rep-seqs.fna file output from deblur denoise-16S
taxa <- assignTaxonomy("../HUM_SD_Gavage_Donor/export_SD_Donor_rep_seqs/dna-sequences.fasta", "../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
write.csv(output, "../HUM_SD_Gavage_Donor/HUM_SD_Gavage_taxonomy_assignments.csv")

## Align DNA sequences --
library(dplyr)
library(tidyr)
fasta<- read.table("../HUM_SD_Gavage_Donor/export_SD_Donor_rep_seqs/dna-sequences.fasta")
fasta2<-fasta %>% separate(V1, c("ASV","QIIME_sequence"), sep = "([>])")
fasta2[fasta2==""]<-NA

ASV <- fasta2$ASV
fasta3 <- data.frame(ASV)
fasta3 <-as.data.frame(fasta3[-1,])
QIIME_seqs <- data.frame(fasta2$QIIME_sequence)
QIIME_seqs <- data.frame(QIIME_seqs[-200,])
fasta3$QIIME_seqs<- QIIME_seqs$QIIME_seqs..200...
finalfasta<-na.omit(fasta3)
write.csv(finalfasta, "../HUM_SD_Gavage_Donor/HUM_SD_Gavage_aligned_dna_sequences.csv")

## Create a taxonomy file that can be imported into QIIME --
finalfasta$ASV <-finalfasta$`fasta3[-1, ]`
output$ASV <- row.names(output)
taxonomy<- merge(finalfasta, output, by="ASV")
taxonomy <- taxonomy %>% select(c("ASV","taxonomy"))

taxonomy_for_qza <- rename(taxonomy, c("Feature ID" = "ASV"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Taxon" = "taxonomy"))

readr::write_delim(taxonomy_for_qza, "../HUM_SD_Gavage_Donor/HUM_SD_taxonomy.tsv", delim ="\t")

## Replace QIIMEID with ASV sequence in the table --

SD_table <- read.delim("../HUM_SD_Gavage_Donor/export_SD_Donor_table/feature-table.tsv")
row.names(SD_table) <- SD_table$OTU.ID
SD_table <- rename(SD_table, c("QIIME_seqs" = "OTU.ID"))

SD_table <- SD_table %>% select(-c(taxonomy))
SD_table <- merge(SD_table,finalfasta, by="QIIME_seqs")
SD_table <- SD_table %>% select(c("ASV", "Hum.SD.Donor.feces"))

## Read in HUM SD Gavage tsv file and restrict the Donor ASV to that 
hum_sd_mouse_table <- read.delim(here("Humanized-Biogeography-Analysis/starting_files/Colonized-ComBat-Adjusted-ASV.tsv"),row.names=1) 
ASV_mouse <- row.names(hum_sd_mouse_table)

SD_table <-SD_table[SD_table$ASV %in% ASV_mouse,]
SD_table <- rename(SD_table, c("#OTU.ID" = "ASV"))
SD_table <- rename(SD_table, c("A001" = "Hum.SD.Donor.feces"))
readr::write_delim(SD_table, "../HUM_SD_Gavage_Donor/HUM_SD_Hoomin.tsv", delim="\t")

