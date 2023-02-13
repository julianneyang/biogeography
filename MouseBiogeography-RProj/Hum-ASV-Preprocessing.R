###Purpose: Use only features from first dataset present in this dataset, then remove 2017_April and 2015_September Sequencing Run
###then perform CombatSeq2
library(dplyr)
library(sva)
###for Humanized

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

ASV<-read.csv("Humanized ASV and Taxonomy Key - feature-table.csv", header=TRUE)

## Remove DSS Samples
metadata<- read.csv("Humanized Metadata - All-Humanized-Metadata.csv")
metadata$Sequencing_Run <- factor(metadata$Sequencing_Run)
metadata$DSS_Treatment <- factor(metadata$DSS_Treatment)

samples <- metadata %>% filter(DSS_Treatment =="No", X.SampleID %in% names(ASV)) %>% pull(X.SampleID)

final_ASV <- ASV[, samples]

write.table(final_ASV,"Humanized_ASV_for_alpha_diversity.tsv", sep="\t")
## Filter by Sequencing Run
row.names(final_ASV) <- ASV$X.OTU.ID
Nov2014 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2014_Nov", X.SampleID %in% names(ASV)) %>%pull(X.SampleID)
Nov2014 <- final_ASV[,Nov2014]
  
Sep2014 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2014_Sept", X.SampleID %in% names(ASV)) %>%pull(X.SampleID)
Sep2014 <- final_ASV[,Sep2014]

Sep2015 <- metadata%>% filter(DSS_Treatment=="No"&Sequencing_Run=="2015_Sept", X.SampleID %in% names(ASV)) %>%pull(X.SampleID)
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
?intersection
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

write.table(combat_adjusted_counts,"Humanized-ComBat-Adjusted-ASV.tsv", sep="\t") 
