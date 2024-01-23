###Purpose: Perform CombatSeq2 after prevalence filtering to 15% for each of two sequencing runs. 
library(dplyr)
library(sva)
library(here)

### For CS Facility Microbiota ---

here::i_am("MouseBiogeography-RProj/CS_Facility_ASV_Preprocessing.R")
ASV<-readr::read_delim(here("CS-Facility-Analysis/CS_Facility_min10k_ASV.tsv"), delim="\t")
names(ASV)
names(ASV)<-gsub("-",".",names(ASV))
ASV <- as.data.frame(ASV)
row.names(ASV)<-ASV$OTU.ID
ASV <- select(ASV, -c("taxonomy","OTU.ID"))
#readr::write_delim(ASV,here("CS-Facility-Analysis/CS_Facility_min10k_ASV.tsv"),delim = "\t")

## Edit Metadata: Remove Samples from Tg mice in Metadata and Drop Samples not in ASV, match ASV sample order 

metadata<- readr::read_csv(here("CS-Facility-Analysis/CS_Facility_Metadata.csv"))

metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Genotype <- factor(metadata$Genotype)
metadata$SampleID
metadata$Site <- factor(metadata$Site)
metadata$Type <- factor(metadata$Type)
metadata$Sex <- factor(metadata$Sex)

metadata <- metadata %>%
  filter(Genotype=="WT") %>%
  filter(SampleID %in% names(ASV))

metadata <- as.data.frame(metadata)
row.names(metadata)<-metadata$SampleID

target <- names(ASV)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

## Filter by Sequencing Run ---
srun1 <- metadata%>% 
  filter(Sequencing_Run=="One", SampleID %in% names(ASV)) %>%pull(SampleID)
srun1 <- ASV[,srun1]


srun2 <- metadata%>% 
  filter(Sequencing_Run=="Two", SampleID %in% names(ASV)) %>%pull(SampleID)
srun2 <- ASV[,srun2]

## Apply feature prevalence filter of 15% to each Sequencing Run dataset
#One
t_df_input_data<-as.data.frame(t(srun1))

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

72*0.15
srun1$prevalence<-prevalence #features present in at least 11 samples our of 72
srun1<- srun1 %>% filter(prevalence>=11) #[Result: 391 features, 72 samples]

#Two
t_df_input_data<-as.data.frame(t(srun2))

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
36*0.15
srun2$prevalence<-prevalence#features present in at least 5 samples out of 36
srun2<- srun2 %>% filter(prevalence>=5) #[Result: 422 features, 36 samples]

##  Find the intersection of all features then query target vector against each dataset, then merge datasets by feature
target <- intersect(row.names(srun1),row.names(srun2)) #[310 features]

srun1<-srun1[row.names(srun1) %in% target,]
srun2<-srun2[row.names(srun2) %in% target,]

row.names(srun1) ==row.names(srun2)

final_ASV <- cbind(srun1, srun2)
final_ASV <-select(final_ASV, -prevalence) #[310 features, 108 samples]

target <- names(final_ASV)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

## Perform CombatSeq2 by Type Sex and Site
batch<- as.character(metadata$Sequencing_Run)
modcombat=model.matrix(~ Type + Sex+ Site,data=metadata)

input=as.matrix(final_ASV)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  
combat_adjusted_counts <- as.data.frame(combat_adjusted_counts)
readr::write_delim(combat_adjusted_counts,here("CS-Facility-Analysis/CS-Facility-ComBat-Adjusted-ASV.tsv"), delim="\t") 

### Make a taxonomy file for QIIME 2 taxa-barplots and for collapsing ASV ---
taxonomy <- readr::read_csv(here("CS-Facility-Analysis/CS_Facility_Deblur_outputs/CS_Facility_taxonomy_assignments.csv"))
taxonomy$ASV <- taxonomy$...1 
aligned<- readr::read_csv(here("CS-Facility-Analysis/CS_Facility_Deblur_outputs/CS_Facility_aligned_dna_sequences.csv"))
aligned$ASV <- aligned$`fasta3[-1, ]`
taxonomy_for_qza <- merge(taxonomy, aligned, by="ASV")
taxonomy_for_qza <- taxonomy_for_qza %>% select(c("ASV","taxonomy"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Feature ID" = "ASV"))
taxonomy_for_qza <- rename(taxonomy_for_qza, c("Taxon" = "taxonomy"))
readr::write_delim(taxonomy_for_qza,here("CS-Facility-Analysis/CS_taxonomy.tsv"),delim="\t")
