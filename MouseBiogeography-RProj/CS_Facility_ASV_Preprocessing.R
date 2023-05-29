###Purpose: Perform CombatSeq2 after prevalence filtering to 15% for each of two sequencing runs. 
library(dplyr)
library(sva)
library(here)

### For CS Facility Microbiota ---

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/CS_Facility_ASV_Preprocessing.R")

here()

ASV<-read.csv(here("CS-Facility-Analysis","CS_Facility_ASV - min10k-feature-table.csv"),row.names=1)
names(ASV)

## Edit Metadata: Remove Samples from Tg mice in Metadata and Drop Samples not in ASV, match ASV sample order 

metadata<- read.csv(here("CS-Facility-Analysis/CS_Facility_Metadata.csv"))

metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$Genotype <- factor(metadata$Genotype)
metadata$X.SampleID
metadata$Site <- factor(metadata$Site)
metadata$Type <- factor(metadata$Type)
metadata$Sex <- factor(metadata$Sex)

metadata <- metadata %>%
  filter(Genotype=="WT") %>%
  filter(SampleID %in% names(ASV))

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

73*0.15
srun1$prevalence<-prevalence #features present in at least 11 samples our of 73
srun1<- srun1 %>% filter(prevalence>=11) #[Result: 391 features, 73 samples]

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
37*0.15
srun2$prevalence<-prevalence#features present in at least 6 samples out of 37
srun2<- srun2 %>% filter(prevalence>=6) #[Result: 396 features, 37 samples]

##  Find the intersection of all features then query target vector against each dataset, then merge datasets by feature
target <- intersect(row.names(srun1),row.names(srun2)) #[306 features]

srun1<-srun1[row.names(srun1) %in% target,]
srun2<-srun2[row.names(srun2) %in% target,]

row.names(srun1) ==row.names(srun2)

final_ASV <- cbind(srun1, srun2)
final_ASV <-select(final_ASV, -prevalence) #[306 features, 108 samples]

target <- names(final_ASV)
metadata <- metadata[match(target, row.names(metadata)),]
target == row.names(metadata)

## Perform CombatSeq2 by Type Sex and Site
batch<- as.character(metadata$Sequencing_Run)
modcombat=model.matrix(~ Type + Sex+ Site,data=metadata)

input=as.matrix(final_ASV)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  

write.table(combat_adjusted_counts,"CS-Facility-Analysis/CS-Facility-ComBat-Adjusted-ASV.tsv", sep="\t") 
