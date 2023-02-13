library(ggplot2)
library(plyr)
library(rlang)
library(rstatix)  
library(nlme)
library(phyloseq)
library(dplyr)
library(cowplot)

getwd()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/ImmDef-Mouse-Biogeography-Analysis/")

chao1<- read.table("alpha_min10k_ValCohort_ImmDef_ASV.qza/chao1.qza_dir/alpha-diversity.tsv",header=T, sep="\t")
otus<-read.table("alpha_min10k_ValCohort_ImmDef_ASV.qza/otus.qza_dir/alpha-diversity.tsv", header=T,sep="\t")
pielou<-read.table("alpha_min10k_ValCohort_ImmDef_ASV.qza/pielou_e.qza_dir/alpha-diversity.tsv", header=T,sep="\t")
shannon<-read.table("alpha_min10k_ValCohort_ImmDef_ASV.qza/shannon.qza_dir/alpha-diversity.tsv", header=T,sep="\t")

temp1<- merge(chao1,otus, by="X")
temp2<- merge(pielou,shannon, by="X")
data<-merge(temp1,temp2,by="X")
data$SampleID <- data$X

write.csv(data, "alpha_diversity_ValCohort.csv")

generate_adiv_plots <- function(input_data, input_metadata, X, Y, facetvariable){
  #read in files
  metadata<- read.csv(input_metadata)
  data<-read.csv(input_data)
  metadata$SampleID <-gsub("-",".",metadata$SampleID)
   #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("Colon", "SI"))
  data$Genotype <- factor(data$Genotype, levels = c("WT", "TCR_KO", "RAGROR"))
  data$Site <- factor(data$Site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
  
 shannon <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{facetvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    theme(legend.position = "top")
 
  shannon + facet_grid(data$Genotype~.)
}

#Aggregate Plot by Sites 
shannon<- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site, shannon, Genotype)
otus <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site, observed_otus, Genotype)
chao1 <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site, chao1, Genotype)
pielou_e <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site, pielou_e, Genotype)

dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General
shannon<- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site_General, shannon, Genotype)
otus <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site_General, observed_otus, Genotype)
chao1 <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site_General, chao1, Genotype)
pielou_e <- generate_adiv_plots("alpha_diversity_ValCohort.csv", "Full-Metadata.csv", Site_General, pielou_e, Genotype)

dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)
#Violin Plot for Mucosal ImmDef Dataset, using subsets
cols <- c("SI" = "red", "Colon" = "blue")
ggplot(data=WTdata,aes(x=Site_General,y=faith_pd, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Mucosal Regional Differences between SI and Colon in WT mice") + ylim(0,125)
ggplot(data=RAGRORdata,aes(x=Site_General,y=faith_pd, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Mucosal Regional Differences between SI and Colon in RAGROR mice") +ylim(0,125)
ggplot(data=TCRKOdata,aes(x=Site_General,y=faith_pd, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Mucosal Regional Differences between SI and Colon in TCRKO mice") +ylim(0,125)

ggplot(data=WTdata,aes(x=Site,y=faith_pd, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences Across Six Sites in WT mice")  + ylim(0,125)
ggplot(data=RAGRORdata,aes(x=Site,y=faith_pd, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences Across Six Sites in RAGROR mice") +ylim(0,125)
ggplot(data=TCRKOdata,aes(x=Site,y=faith_pd, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences Across Six Sites in TCRKO mice") +ylim(0,125)


#Using full dataset 
faith_pd<-ggplot(data=data_box_plot,aes(x=Site_General,y=faith_pd, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences between SI and Colon")
faith_pd<-ggplot(data=data_box_plot,aes(x=Site,y=faith_pd, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences Across Six Sites") 
faith_pd<-ggplot(data=data_box_plot,aes(x=Type_Site,y=faith_pd, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Faith's Index: Luminal Regional Differences Across Six Sites") 
faith_pd + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
                 
evenness<-ggplot(data=data_box_plot,aes(x=Site,y=pielou_evenness, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Pielou's Evenness: Luminal Regional Differences Across Six Sites") 
evenness<-ggplot(data=data_box_plot,aes(x=Site_General,y=pielou_evenness, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Pielou's Evenness: Luminal Regional Differences between SI and Colon") 
evenness<-ggplot(data=data_box_plot,aes(x=Type_Site,y=pielou_evenness, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Pielou's Evenness: Luminal Regional Differences Across Six Sites") 
evenness + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

shannon<-ggplot(data=data_box_plot,aes(x=Site,y=shannon_entropy, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Shannon_Entropy: Luminal Regional Differences Across Six Sites") 
shannon<-ggplot(data=data_box_plot,aes(x=Site_General,y=shannon_entropy, shape=Site_General, color=Site_General)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Shannon_Entropy: Luminal Regional Differences between SI and Colon") 
shannon<-ggplot(data=data_box_plot,aes(x=Type_Site,y=shannon_entropy, shape=Site, color=Site)) + geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+mytheme +ggtitle("Shannon_Entropy: Luminal Regional Differences Across Six Sites") 
shannon + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


#Linear mixed effects models, Sequencing Run as fixed effect
names(data)
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Genotype*Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Genotype*Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Genotype*Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Genotype*Site, random = ~1|MouseID_Original, data=data)
summary(output)

names(data)
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Genotype*Site_General, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Genotype*Site_General, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Genotype*Site_General, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Genotype*Site_General, random = ~1|MouseID_Original, data=data)
summary(output)
