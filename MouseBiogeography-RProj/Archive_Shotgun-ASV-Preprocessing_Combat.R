###Purpose: Use only features from first dataset present in this dataset, then remove 2017_April and 2015_September Sequencing Run
###then perform CombatSeq2
library(dplyr)
library(sva)
library(here)
library(dplyr)
library(dada2)

here::i_am("MouseBiogeography-RProj/Hum-ASV-Preprocessing.R")


counts <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_ASV.tsv"))
df_input_data <- as.data.frame(counts)
row.names(df_input_data)<-df_input_data$Species
final_ASV <- df_input_data %>% select(-c("Species"))

#summary(colSums(ASV))
#ASV<- ASV[colSums(ASV) >= 10000]


## Remove DSS Samples
metadata<- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim = "\t")
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)


## Filter by Sequencing Run
Nov2014 <- metadata%>% filter(Sequencing_Run=="One", sampleid %in% names(final_ASV)) %>%pull(sampleid)
Nov2014 <- final_ASV[,Nov2014]
Sep2014 <- metadata%>% filter(Sequencing_Run=="Two", sampleid %in% names(ASV)) %>%pull(sampleid)
Sep2014 <- final_ASV[,Sep2014]

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

70*0.15
Nov2014$prevalence<-prevalence #features present in at least 11 samples our of 70
Nov2014<- Nov2014 %>% filter(prevalence>=11) #[Result: 509 features, 100 samples]

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
34*0.15
Sep2014$prevalence<-prevalence#features present in at least 5 samples our of 34
Sep2014<- Sep2014 %>% filter(prevalence>=5) #[Result: 303 features,34 samples]


##Find the intersection of all features then query target vector against each dataset, then merge datasets by feature
target <- intersect(row.names(Nov2014),row.names(Sep2014)) #[283 features]

Nov2014<-Nov2014[row.names(Nov2014) %in% target,]
Sep2014<-Sep2014[row.names(Sep2014) %in% target,]

row.names(Nov2014) ==row.names(Sep2014)

final_humanized_ASV <- cbind(Nov2014,Sep2014)
row.names(final_humanized_ASV) <-row.names(Sep2014)
final_humanized_ASV <-select(final_humanized_ASV, -prevalence) #[263 features, 213 samples]


## Perform CombatSeq2 by Type and Site
batch<- as.character(metadata$Sequencing_Run)
modcombat=model.matrix(~ Sex + Site,data=metadata)

input=as.matrix(final_humanized_ASV)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  

combat_adjusted_counts <- as.data.frame(combat_adjusted_counts)

combat_adjusted_counts$OTU.ID <- rownames(combat_adjusted_counts)
combat_adjusted_counts <- combat_adjusted_counts %>%
  select(OTU.ID, everything())

readr::write_delim(combat_adjusted_counts,here("Shotgun/Shotgun-ComBat-Adjusted-ASV.tsv"), delim="\t") 
