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
here::i_am("MouseBiogeography-RProj/CS_Facility-AlphaDiversity.R")                 
otus<- read.table("CS-Facility-Analysis/alpha_min_10000_table/otus_dir/alpha-diversity.tsv",header=T, sep="\t")
chao1<-read.table("CS-Facility-Analysis/alpha_min_10000_table/chao1_dir/alpha-diversity.tsv", header=T,sep="\t")
pielou<-read.table("CS-Facility-Analysis/alpha_min_10000_table/pielou_e_dir/alpha-diversity.tsv", header=T,sep="\t")
shannon<-read.table("CS-Facility-Analysis/alpha_min_10000_table/shannon_dir/alpha-diversity.tsv", header=T,sep="\t")

temp1<- merge(chao1,otus, by="X")
temp2<- merge(pielou,shannon, by="X")
data<-merge(temp1,temp2,by="X")
data$SampleID <- gsub("-",".",data$X)

write_rds(data, "CS-Facility-Analysis/alpha_diversity_CS_Facility.RDS")

generate_adiv_plots <- function(input_data, X, Y){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))

  shannon <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    theme(legend.position = "top")
}

#read in files
metadata<- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv", row.names=1)
metadata$SampleID<-row.names(metadata)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

#data <- data %>% pivot_longer(cols=c("chao1", "shannon","pielou_e","observed_otus"),names_to="adiv")

luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

#Aggregate Plot by Sites Luminal Data
shannon <- generate_adiv_plots(luminaldata, Site, shannon) +
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(luminaldata, Site, observed_otus)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(luminaldata, Site, chao1)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(luminaldata, Site, pielou_e)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)


#Aggregate Plot by Sites Mucosal Data
shannon <- generate_adiv_plots(mucosaldata, Site, shannon) +
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(mucosaldata, Site, observed_otus)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(mucosaldata, Site, chao1)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(mucosaldata, Site, pielou_e)+
  stat_compare_means(comparisons = list(c("Distal_Colon", "Proximal_Colon"),
                                        c("Distal_Colon", "Cecum"),
                                        c("Distal_Colon", "Ileum"),
                                        c("Distal_Colon", "Jejunum"),
                                        c("Distal_Colon", "Duodenum")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)

dev.new(width=15, height=20)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General Mucosal Data
shannon <- generate_adiv_plots(mucosaldata, Site_General, shannon) +
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(mucosaldata, Site_General, observed_otus)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(mucosaldata, Site_General, chao1)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(mucosaldata, Site_General, pielou_e)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Site_General Luminal Data

shannon <- generate_adiv_plots(luminaldata, Site_General, shannon) +
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(luminaldata, Site_General, observed_otus)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(luminaldata, Site_General, chao1)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(luminaldata, Site_General, pielou_e)+
  stat_compare_means(comparisons = list(c("Colon","SI")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Type Colon Data
shannon <- generate_adiv_plots(colondata, Type, shannon) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(colondata, Type, observed_otus)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(colondata, Type, chao1)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(colondata, Type, pielou_e)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

#Aggregate Plot by Type SI Data
shannon <- generate_adiv_plots(SIdata, Type, shannon) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(SIdata, Type, observed_otus)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(SIdata, Type, chao1)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(SIdata, Type, pielou_e)+
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)


##Aggregate Plot by Type SI Data by site
shannon <- generate_adiv_plots(SIdata, Type, shannon) + facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(SIdata, Type, observed_otus)+ facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(SIdata, Type, chao1)+facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(SIdata, Type, pielou_e)+facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

##Aggregate Plot by Type Colon Data by site
shannon <- generate_adiv_plots(colondata, Type, shannon) + facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
otus<- generate_adiv_plots(colondata, Type, observed_otus)+ facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
chao1<- generate_adiv_plots(colondata, Type, chao1)+facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
pielou_e<- generate_adiv_plots(colondata, Type, pielou_e)+facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Luminal", "Mucosal")),method="wilcox", vjust=0.3,label="p.signif",step.increase=0.05)
dev.new(width=15, height=10)
plot_grid(shannon, otus, chao1, pielou_e, rows=2, cols=2)

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
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=luminaldata)
summary(output)

# mucosal - Site_General
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site_General, random = ~1|MouseID, data=mucosaldata)
summary(output)

# luminal - Site
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=luminaldata)
summary(output)

# mucosal - Site 
output=lme(fixed= shannon ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= chao1 ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)
output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=mucosaldata)
summary(output)

