library(ggplot2)
library(plyr)
library(agricolae)
library(rstatix)  
library(nlme)
library(phyloseq)
library(dplyr)
library(nlme)
library(ggsignif)
library(here)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d1123")
filepath <- "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/"
chao1<- read.table(here(paste0(filepath,"chao1_dir/alpha-diversity.tsv")),header=T, sep="\t")
otus<- read.table(here(paste0(filepath,"otus_dir/alpha-diversity.tsv")),header=T, sep="\t")
pielou<-read.table(here(paste0(filepath,"pielou_e_dir/alpha-diversity.tsv")),header=T, sep="\t")
shannon<-read.table(here(paste0(filepath,"shannon_dir/alpha-diversity.tsv")),header=T, sep="\t")

temp1<- merge(chao1,otus, by="X")
temp2<- merge(pielou,shannon, by="X")
data<-merge(temp1,temp2,by="X")
data$SampleID <- data$X

write.csv(data, "C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d1123/alpha_diversity_Regional.csv")

generate_adiv_plots <- function(input_data, X, Y){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
  
 ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    theme(legend.position = "top")

}

#merge metadata with alpha diversity
data$SampleID <- gsub("-",".",data$X)
metadata<- read.csv("Regional-Combat-Metadata.csv", header=TRUE)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Line = factor(data$Line)
sapply(data,levels)

luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

names(data)
#Aggregate Plot by Sites Luminal Data
shannon<- generate_adiv_plots(luminaldata, Site, shannon) + 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6) 
otus<- generate_adiv_plots(luminaldata, Site, observed_otus)+ 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6) 
chao1<- generate_adiv_plots(luminaldata, Site, chao1)+ 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6) 
pielou_e<- generate_adiv_plots(luminaldata, Site, pielou_e)+ 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6) 

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2, align = 'hv')


#Aggregate Plot by Sites Mucosal Data
shannon<- generate_adiv_plots(mucosaldata, Site, shannon) + 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
otus<- generate_adiv_plots(mucosaldata, Site, observed_otus) + 
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
chao1<- generate_adiv_plots(mucosaldata, Site, chao1) +
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)
pielou_e<- generate_adiv_plots(mucosaldata, Site, pielou_e) +
  geom_signif(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                 c("Distal_Colon", "Cecum"),
                                 c("Distal_Colon", "Ileum"),
                                 c("Distal_Colon", "Jejunum"),
                                 c("Distal_Colon", "Duodenum")),map_signif_level = TRUE, textsize = 6)

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2, align = 'hv')

#Aggregate Plot by Site_General Luminal Data
shannon<- generate_adiv_plots(luminaldata, Site_General, shannon)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(luminaldata, Site_General, observed_otus)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(luminaldata, Site_General, chao1)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(luminaldata, Site_General, pielou_e)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2, align = 'hv')

#Aggregate Plot by Site_General Mucosal Data
shannon<- generate_adiv_plots(mucosaldata, Site_General, shannon) +
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(mucosaldata, Site_General, observed_otus)+
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(mucosaldata, Site_General, chao1) +
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(mucosaldata, Site_General, pielou_e) +
  geom_signif(comparisons = list(c("Colon", "SI")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)
?geom_signif
#Aggregate Plot by Type Colon Data
shannon<- generate_adiv_plots(colondata, Type, shannon) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(colondata, Type, observed_otus) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test")
chao1 <- generate_adiv_plots(colondata, Type, chao1) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(colondata, Type, pielou_e)+
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 

dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Type Colon Data by Site
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

#Aggregate Plot by Type SI Data 
shannon<- generate_adiv_plots(SIdata, Type, shannon) + 
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
otus <- generate_adiv_plots(SIdata, Type, observed_otus) + 
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
chao1 <- generate_adiv_plots(SIdata, Type, chao1) + 
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
pielou_e <- generate_adiv_plots(SIdata, Type, pielou_e) +
  geom_signif(comparisons = list(c("Luminal", "Mucosal")),map_signif_level = TRUE, textsize = 6, test="wilcox.test") 
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2,cols=2)


#Aggregate Plot by Type SI Data by site
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

##Linear mixed effects models, Site differences 
data <- read.csv(here(paste0(filepath,"alpha_diversity_Regional.csv")),row.names=1)
metadata<- read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv"),row.names=1)
metadata$SampleID <- row.names(metadata)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
sapply(data,levels)

duodata<-filter(data, Site=="Duodenum")
jejdata<-filter(data, Site =="Jejunum")
iledata<-filter(data, Site =="Ileum")
cecdata<-filter(data, Site =="Cecum")
pcdata<-filter(data, Site =="Proximal_Colon")
dcdata<-filter(data, Site =="Distal_Colon")

# Luminal Data
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Line + Site_General, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Site_General, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Line + Site_General, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Site_General, random = ~1|MouseID_Line, data=luminaldata)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex +  Line + Site, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Site, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Line + Site, random = ~1|MouseID_Line, data=luminaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex +Line + Site, random = ~1|MouseID_Line, data=luminaldata)
summary(output)

# Mucosal Data
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Line + Site_General, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line +Site_General, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex +Line + Site_General, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line +Site_General, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)

output=lme(fixed= shannon ~ Sequencing_Run + Sex + Line + Site, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line +Site, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Line +Site, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex +Line +Site, random = ~1|MouseID_Line, data=mucosaldata)
summary(output)

# Sample Type Differences
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=duodata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=duodata)
summary(output)

output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=jejdata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=jejdata)
summary(output)

output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=iledata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line+ Type, random = ~1|(MouseID_Line), data=iledata)
summary(output)

output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line+ Type, random = ~1|(MouseID_Line), data=cecdata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=cecdata)
summary(output)

output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=pcdata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=pcdata)
summary(output)

output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=dcdata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=dcdata)
summary(output)
