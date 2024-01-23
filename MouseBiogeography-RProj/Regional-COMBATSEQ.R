BiocManager::install(c("sva"))
library("sva")
library("dplyr")
library("plyr")
##NOTE: Metadata file needs to exactly match the data file in the order of samples
here::i_am("MouseBiogeography-RProj/Regional-COMBATSEQ.R")
newmeta=read.csv("../Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Sequencing Run Overlap ASV and Metadata - Regional-Overlap-Metadata.csv", header=T,row.names=1)

newmeta=as.data.frame(newmeta)
filtereddata=read.csv("../Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Sequencing Run Overlap ASV and Metadata - Regional-Overlap-ASV.csv", row.names=1)
#filtereddata= filtereddata[,-c(90)] #omit taxonomy column 
warnings()
batch=newmeta$Sequencing_Run
modcombat=model.matrix(~ Type + Site,data=newmeta)
#group=newmeta$Group


input=as.matrix(filtereddata)
combat_adjusted_counts <- ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  
combat_adjusted_counts <- as.data.frame(combat_adjusted_counts)

### Make a taxonomy file for QIIME 2 taxa-barplots and for collapsing ASV ---
taxonomy <- readr::read_csv(here("Regional-Mouse-Biogeography-Analysis/UCLA_taxonomy_assignments.csv"))
aligned<- readr::read_csv(here("Regional-Mouse-Biogeography-Analysis/UCLA_aligned_dna_sequences.csv"))
aligned$ASV <- aligned$`fasta3[-1, ]`
taxonomy_for_qza <- merge(taxonomy, aligned, by="ASV")
taxonomy_for_qza <- taxonomy_for_qza %>% select(c("ASV","taxonomy"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Feature ID" = "ASV"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Taxon" = "taxonomy"))
readr::write_delim(taxonomy_for_qza,here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/UCLA_taxonomy.tsv"),delim="\t")

### Switch QIIME SEqs for actual ASV sequences --
combat_adjusted_counts$QIIME_seqs <- rownames(combat_adjusted_counts)
taxonomy <- aligned %>% select(c("QIIME_seqs","ASV"))
combat_adjusted_counts_ASV <-merge(combat_adjusted_counts, taxonomy, by="QIIME_seqs")
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>% select(-c("QIIME_seqs"))
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>%
  select(ASV, everything())
combat_adjusted_counts_ASV <- rename(combat_adjusted_counts_ASV, c("#OTU.ID" = "ASV"))

readr::write_delim(combat_adjusted_counts_ASV,here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/UCLA-ComBat-Adjusted-ASV.tsv"), delim="\t") 

