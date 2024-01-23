###Purpose: Use only features from first dataset present in this dataset, then remove 2017_April and 2015_September Sequencing Run
###then perform CombatSeq2
library(dplyr)
library(sva)
library(here)
library(dplyr)
library(dada2)

here::i_am("MouseBiogeography-RProj/Hum-ASV-Preprocessing.R")
ASV<-readr::read_csv(here("Humanized-Biogeography-Analysis/Humanized ASV and Taxonomy Key - feature-table.csv"))
ASV <- as.data.frame(ASV)
row.names(ASV) <- ASV$OTU.ID
ASV <- ASV %>% select(-c("OTU.ID"))
summary(colSums(ASV))
ASV<- ASV[colSums(ASV) >= 10000]


## Remove DSS Samples
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim = "\t")
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$DSS_Treatment <- factor(metadata$DSS_Treatment)

samples <- metadata %>% filter(DSS_Treatment =="No", SampleID %in% names(ASV)) %>% pull(SampleID)
metadata <- metadata %>% filter(DSS_Treatment =="No", SampleID %in% names(ASV))
metadata<- readr::write_delim(metadata, here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim = "\t")

final_ASV <- ASV[, samples]
final_ASV <- as.data.frame(final_ASV)

write.table(final_ASV,"Humanized_ASV_for_alpha_diversity.tsv", sep="\t")
## Filter by Sequencing Run
Nov2014 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2014_Nov", SampleID %in% names(ASV)) %>%pull(SampleID)
Nov2014 <- final_ASV[,Nov2014]
Sep2014 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2014_Sept", SampleID %in% names(ASV)) %>%pull(SampleID)
Sep2014 <- final_ASV[,Sep2014]

Sep2015 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2015_Sept", SampleID %in% names(ASV)) %>%pull(SampleID)
Sep2015 <- final_ASV[,Sep2015]

## Apply feature prevalence filter of 15% to each Sequencing Run dataset
#Nov2014
t_df_input_data<-as.data.frame(t(Nov2014))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}

Nov2014$prevalence<-prevalence #features present in at least 15 samples our of 100
Nov2014<- Nov2014 %>% filter(prevalence>=15) #[Result: 509 features, 100 samples]

#Sep2014
t_df_input_data<-as.data.frame(t(Sep2014))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}
70*0.15
Sep2014$prevalence<-prevalence#features present in at least 11 samples our of 70
Sep2014<- Sep2014 %>% filter(prevalence>=11) #[Result: 488 features, 70 samples]

#Sep2015
t_df_input_data<-as.data.frame(t(Sep2015))

ctr= 0
prevalence <- vector(mode="numeric")

for(i in 1:ncol(t_df_input_data)){
  v<-t_df_input_data %>% pull(i)
  for(j in 1:length(v)){
    if (t_df_input_data[j,i]>0){
      ctr=1+ctr
    }
    else {
      ctr=ctr
    }
  }
  prevalence<-append(prevalence,ctr)
  ctr=0
}

Sep2015$prevalence<-prevalence#features present in at least 6 samples out of 43
Sep2015<- Sep2015 %>% filter(prevalence>=6) #[Result: 363 features, 43 samples]

##Find the intersection of all features then query target vector against each dataset, then merge datasets by feature
target1 <- intersect(row.names(Nov2014),row.names(Sep2014)) #[452 features]
target <- intersect(target1, row.names(Sep2015)) #[263 features]

Nov2014<-Nov2014[row.names(Nov2014) %in% target,]
Sep2014<-Sep2014[row.names(Sep2014) %in% target,]
Sep2015<-Sep2015[row.names(Sep2015) %in% target,]

row.names(Nov2014) ==row.names(Sep2014)
row.names(Nov2014) ==row.names(Sep2015)

final_humanized_ASV <- cbind(Nov2014,Sep2014,Sep2015)
final_humanized_ASV <-select(final_humanized_ASV, -prevalence) #[263 features, 213 samples]


## Perform CombatSeq2 by Type and Site
newmeta <- metadata %>% filter(DSS_Treatment =="No")

batch<- as.character(newmeta$Sequencing_Run)
modcombat=model.matrix(~ Microbiota + Type + Sex + Site,data=newmeta)

input=as.matrix(final_humanized_ASV)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  

combat_adjusted_counts <- as.data.frame(combat_adjusted_counts)

### Redo taxonomy assignments with Silva 138.1

seqs <-readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized_aligned_dna_sequences.tsv"),delim="\t")
seqs <- seqs$ASV
set.seed(100)
taxa <- assignTaxonomy(seqs, 
                       "../16S_Taxonomy_Classifiers/silva_nr99_v138.1_wSpecies_train_set.fa.gz", multithread=TRUE)  # update directory with the deblur all.seqs.fa output file and with the most recent Silva 99% OTU database  
taxa <- addSpecies(taxa, "../16S_Taxonomy_Classifiers/silva_species_assignment_v138.1.fa.gz")
taxa[is.na(taxa)] <- ""
taxonomy<-paste("k__",taxa[,1],"; ","p__",taxa[,2],"; ","c__",taxa[,3],"; ","o__",taxa[,4],"; ","f__",taxa[,5],"; ","g__",taxa[,6],"; ","s__",taxa[,7],sep="")
taxonomy <- as.data.frame(taxonomy)
output<-cbind(taxa, taxonomy)
output$ASV <- row.names(output)
readr::write_csv(output, here("Humanized-Biogeography-Analysis/starting_files/Humanized_taxonomy_assignments.csv"))

### Make a taxonomy file for QIIME 2 taxa-barplots and for collapsing ASV ---
taxonomy <- readr::read_csv(here("Humanized-Biogeography-Analysis/starting_files/Humanized_taxonomy_assignments.csv"))
aligned<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized_aligned_dna_sequences.tsv"),delim="\t")
taxonomy_for_qza <- merge(taxonomy, aligned, by="ASV")
taxonomy_for_qza <- taxonomy_for_qza %>% select(c("ASV","taxonomy"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Feature ID" = "ASV"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Taxon" = "taxonomy"))
readr::write_delim(taxonomy_for_qza,here("Humanized-Biogeography-Analysis/starting_files/Humanized_taxonomy.tsv"),delim="\t")

### Switch QIIME SEqs for actual ASV sequences --
combat_adjusted_counts$QIIME_seqs <- rownames(combat_adjusted_counts)
taxonomy <- aligned %>% select(c("QIIME_seqs","ASV"))
combat_adjusted_counts_ASV <-merge(combat_adjusted_counts, taxonomy, by="QIIME_seqs")
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>% select(-c("QIIME_seqs"))
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>%
  select(ASV, everything())
combat_adjusted_counts_ASV <- rename(combat_adjusted_counts_ASV, c("#OTU.ID" = "ASV"))

readr::write_delim(combat_adjusted_counts_ASV,here("Humanized-Biogeography-Analysis/starting_files/Colonized-ComBat-Adjusted-ASV.tsv"), delim="\t") 
