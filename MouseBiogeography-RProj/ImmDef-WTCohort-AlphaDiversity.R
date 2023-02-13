library(ggplot2)
library(plyr)
library(rlang)
library(rstatix)  
library(nlme)
library(phyloseq)
library(dplyr)
library(cowplot)

getwd()
here::i_am("MouseBiogeography-RProj/ImmDef-WTCohort-AlphaDiversity.R")

chao1<- read.table("ImmDef-Mouse-Biogeography-Analysis/d10625_alpha_min10000_WTCohort_ASV_for_alpha_diversity/chao1_dir/alpha-diversity.tsv",header=T, sep="\t")
otus<-read.table("ImmDef-Mouse-Biogeography-Analysis/d10625_alpha_min10000_WTCohort_ASV_for_alpha_diversity/otus_dir/alpha-diversity.tsv", header=T,sep="\t")
pielou<-read.table("ImmDef-Mouse-Biogeography-Analysis/d10625_alpha_min10000_WTCohort_ASV_for_alpha_diversity/pielou_e_dir/alpha-diversity.tsv", header=T,sep="\t")
shannon<-read.table("ImmDef-Mouse-Biogeography-Analysis/d10625_alpha_min10000_WTCohort_ASV_for_alpha_diversity/shannon_dir/alpha-diversity.tsv", header=T,sep="\t")

temp1<- merge(chao1,otus, by="X")
temp2<- merge(pielou,shannon, by="X")
data<-merge(temp1,temp2,by="X")
data$SampleID <- data$X

write.csv(data, "ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.RDS")

generate_adiv_plots <- function(input_data, input_metadata, X, Y, facetvariable){
  #read in files
  metadata<- read.csv(input_metadata)
  data<-read.csv(input_data)
  metadata$SampleID <-gsub("-",".",metadata$SampleID)
   #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Genotype <- factor(data$Genotype, levels = c("WT", "TCR_KO", "RAGROR"))
  data$Site <- factor(data$Site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
  
 shannon <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{facetvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    theme(legend.position = "top")

}

#Aggregate Plot by Sites 
shannon<- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, shannon, Site)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.65,label="p.signif")
otus<- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, observed_otus, Site)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.65,label="p.signif")
chao1 <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, chao1, Site)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.65,label="p.signif")
pielou_e <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, pielou_e, Site)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.65,label="p.signif")
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General
shannon<- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site_General, shannon, Site_General)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site_General, observed_otus, Site_General)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site_General, chao1,Site_General)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site_General, pielou_e, Site_General)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 

dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)


#Linear mixed effects models, Sequencing Run as fixed effect
#read in files
metadata<- read.csv("ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv")
data<-read.csv("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv")
metadata$SampleID <-gsub("-",".",metadata$SampleID)
#append metadata
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

#Ensure correct ordering of levels 
data$Site_General <- factor(data$Site_General, levels = c("Colon", "SI"))
data$Site <- factor(data$Site, levels = c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))

names(data)
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site, random = ~1|MouseID_Original, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex*Site, random = ~1|MouseID_Original, data=data)
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
