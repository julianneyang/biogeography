library(ggplot2)
library(plyr)
library(rlang)
library(rstatix)  
library(nlme)
library(phyloseq)
library(dplyr)
library(cowplot)

getwd()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/")

chao1<- read.table("alpha_min10k_Humanized_ASV_for_alpha_diversity.qza/chao1.qza_dir/alpha-diversity.tsv",header=T, sep="\t")
otus<-read.table("alpha_min10k_Humanized_ASV_for_alpha_diversity.qza/otus.qza_dir/alpha-diversity.tsv", header=T,sep="\t")
pielou<-read.table("alpha_min10k_Humanized_ASV_for_alpha_diversity.qza/pielou_e.qza_dir/alpha-diversity.tsv", header=T,sep="\t")
shannon<-read.table("alpha_min10k_Humanized_ASV_for_alpha_diversity.qza/shannon.qza_dir/alpha-diversity.tsv", header=T,sep="\t")

temp1<- merge(chao1,otus, by="X")
temp2<- merge(pielou,shannon, by="X")
data<-merge(temp1,temp2,by="X")
data$SampleID <- data$X

write.csv(data, "alpha_diversity_Humanized.csv")

generate_adiv_plots <- function(input_data, X, Y, facetvariable){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("Colon", "SI"))
  data$Microbiota <- factor(data$Microbiota, levels = c("Cedars_SPF", "Humanized"))
  data$Site <- factor(data$Site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))

  shannon <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    theme(legend.position = "top")
  shannon+facet_grid(~Microbiota)
}

#read in files
metadata<- read.csv("Humanized Metadata - All-Humanized-Metadata (1).csv")
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate
data$Microbiota
luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

#Aggregate Plot by Sites Luminal Data
shannon<- generate_adiv_plots(luminaldata, Site, shannon, Microbiota) +
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
otus<- generate_adiv_plots(luminaldata, Site, observed_otus, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
chao1<- generate_adiv_plots(luminaldata, Site, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
pielou_e<- generate_adiv_plots(luminaldata, Site, pielou_e, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)


#Aggregate Plot by Sites Mucosal Data
shannon<- generate_adiv_plots(mucosaldata, Site, shannon, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
otus<- generate_adiv_plots(mucosaldata, Site, observed_otus, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
chao1<- generate_adiv_plots(mucosaldata, Site, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
pielou_e<- generate_adiv_plots(mucosaldata, Site, pielou_e, Microbiota)+
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General Luminal Data
shannon<- generate_adiv_plots(luminaldata, Site_General, shannon, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(luminaldata, Site_General, observed_otus, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(luminaldata, Site_General, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(luminaldata, Site_General, pielou_e, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 

dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General Mucosal Data

shannon<- generate_adiv_plots(mucosaldata, Site_General, shannon, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(mucosaldata, Site_General, observed_otus, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(mucosaldata, Site_General, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(mucosaldata, Site_General, pielou_e, Microbiota)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Type Colon Data
shannon<- generate_adiv_plots(colondata, Type, shannon, Microbiota)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test")
otus <- generate_adiv_plots(colondata, Type, observed_otus, Microbiota)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test")
chao1 <- generate_adiv_plots(colondata, Type, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test")
pielou_e <- generate_adiv_plots(colondata, Type, pielou_e, Microbiota)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test")
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Type SI Data
shannon<- generate_adiv_plots(SIdata, Type, shannon, Microbiota) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(SIdata, Type, observed_otus, Microbiota) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(SIdata, Type, chao1, Microbiota)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(SIdata, Type, pielou_e, Microbiota)+
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

##Aggregate Plot by Type SI Data by site
SIdata<-filter(data, Site_General =="SI" & Microbiota == "Cedars_SPF")
shannon<- generate_adiv_plots(SIdata, Type, shannon) + facet_grid(~Site) + 
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(SIdata, Type, observed_otus) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(SIdata, Type, chao1) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(SIdata, Type, pielou_e) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, ncol=4)

SIdata<-filter(data, Site_General =="SI" & Microbiota == "Humanized")
shannon<- generate_adiv_plots(SIdata, Type, shannon) + facet_grid(~Site) + 
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(SIdata, Type, observed_otus) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(SIdata, Type, chao1) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(SIdata, Type, pielou_e) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, ncol=4)

##Aggregate Plot by Type Colon Data by site
colondata<-filter(data, Site_General =="Colon" & Microbiota=="Cedars_SPF")
shannon<- generate_adiv_plots(colondata, Type, shannon) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(colondata, Type, observed_otus)+ facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(colondata, Type, chao1)+ facet_grid(~Site)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(colondata, Type, pielou_e)+ facet_grid(~Site)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, ncol=4)

colondata<-filter(data, Site_General =="Colon" & Microbiota=="Humanized")
shannon<- generate_adiv_plots(colondata, Type, shannon) + facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(colondata, Type, observed_otus)+ facet_grid(~Site) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(colondata, Type, chao1)+ facet_grid(~Site)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(colondata, Type, pielou_e)+ facet_grid(~Site)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, ncol=4)

#Linear mixed effects models, full dataset
names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
sapply(data,levels)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Type + Microbiota*Site_General, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type + Microbiota*Site_General, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Type + Microbiota*Site_General, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type + Microbiota*Site_General, random = ~1|MouseID, data=data)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Type + Microbiota*Site, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type + Microbiota*Site, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Type + Microbiota*Site, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type + Microbiota*Site, random = ~1|MouseID, data=data)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site + Microbiota*Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site + Microbiota*Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site + Microbiota*Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site + Microbiota*Type, random = ~1|MouseID, data=data)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Microbiota+ Sex+ Site + Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Microbiota + Sex + Site + Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site + Microbiota + Type, random = ~1|MouseID, data=data)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site + Microbiota + Type, random = ~1|MouseID, data=data)
summary(output)

#Linear mixed effects models, SPF Dataset
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
sapply(data,levels)
spf<-filter(data, Microbiota=="Cedars_SPF")
luminalspf <- filter(luminaldata, Microbiota=="Cedars_SPF")
luminalspf$Site = factor(luminalspf$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))

output=lme(fixed= shannon ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=luminalspf)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=luminalspf)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=luminalspf)
summary(output)

mucosalspf <- filter(mucosaldata, Microbiota=="Cedars_SPF")
mucosalspf$Site = factor(mucosalspf$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))

output=lme(fixed= shannon ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=mucosalspf)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=mucosalspf)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=mucosalspf)
summary(output)

#Linear mixed effects models, Humanized Dataset
humanized<-filter(data, Microbiota =="Humanized")
humlum <- filter(luminaldata, Microbiota=="Humanized")
humlum$Site = factor(humlum$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))

output=lme(fixed= shannon ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=humlum)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=humlum)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=humlum)
summary(output)

hummuc <- filter(mucosaldata, Microbiota=="Humanized")
hummuc$Site = factor(hummuc$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))

output=lme(fixed= shannon ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site_General, random = ~1|MouseID, data=hummuc)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=hummuc)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex  + Site, random = ~1|MouseID, data=hummuc)
summary(output)
