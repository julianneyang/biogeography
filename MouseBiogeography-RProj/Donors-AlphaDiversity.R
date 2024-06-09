library(ggplot2)
library(plyr)
library(rlang)
library(rstatix)  
library(nlme)
library(phyloseq)
library(dplyr)
library(cowplot)
library(ggpubr)

here()
here::i_am("MouseBiogeography-RProj/Donors-AlphaDiversity.R")                 
chao<- readr::read_delim(here("Donors-Analysis/alpha_diversity/alpha_Donors_ASV_min10000/chao1_dir/alpha-diversity.tsv"))
otu <- readr::read_delim(here("Donors-Analysis/alpha_diversity/alpha_Donors_ASV_min10000/otus_dir/alpha-diversity.tsv"))
pe <- readr::read_delim(here("Donors-Analysis/alpha_diversity/alpha_Donors_ASV_min10000/pielou_e_dir/alpha-diversity.tsv"))
shannon <- readr::read_delim(here("Donors-Analysis/alpha_diversity/alpha_Donors_ASV_min10000/shannon_dir/alpha-diversity.tsv"))

temp1<- merge(chao,otu, by="...1")
temp2<- merge(pe,shannon, by="...1")
data<-merge(temp1,temp2,by="...1")
data$SampleID <- data$...1
data$SampleID <- gsub("-",".",data$SampleID)

readr::write_rds(data, "Donors-Analysis/alpha_diversity/alpha_diversity.RDS")

# Merge metadata
data<-readr::read_rds(here("Donors-Analysis/alpha_diversity/alpha_diversity.RDS"))
metadata<- readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"))
metadata$SampleID <- gsub("-",".",metadata$SampleID)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

#Linear mixed effects models, full dataset
names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$MouseID <- factor(data$MouseID)
sapply(data,levels)

luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

# luminal - Site_General
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site_General, random = ~1|(Donor_ID/MouseID), data=luminaldata)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)

# luminal - Site
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|(Donor_ID/MouseID), data=luminaldata)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)


# luminal - Site_General
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)

# mucosal - Site_General
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_evenness~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)

# luminal - Site
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)

# mucosal - Site 
output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)

