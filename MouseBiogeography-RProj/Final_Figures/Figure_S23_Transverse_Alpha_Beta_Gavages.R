###### The Biggest and Most Tragic Rearrangement of Mouse Biogeography ---
### Date: 06.15.2024
### Figure Number: Supp 19
### Figure Contents: SPF Gavage and HUM Gavage Transverse alpha and beta diversity 
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(Microbiome.Biogeography)
library(devtools)
library(roxygen2)
find_rtools()

setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")


### Alpha Diversity ---
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_Transverse_Alpha_Beta_Gavages.R")

# HUM V Gavage  
data<-readRDS("Donors-Analysis/alpha_diversity/alpha_diversity.RDS")
metadata<- read.delim("Donors-Analysis/starting_files/Donors_Metadata.tsv")

data <- data %>% filter(observed_features <= 600)

hum_v_gavage_otus_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, observed_features, Type, 0, 600)+ 
  facet_grid(~Site) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM MD Gavage")+
  labs(y="# ASVs", x="")+
  stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)
hum_v_gavage_pe_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, pielou_evenness, Type, 0, 1)+ 
  facet_grid(~Site) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")+
  stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)


# SPF Gavage 
data<-read.csv("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv", header=TRUE, row.names=1)
metadata<- read.delim("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv")
intermediate<- metadata %>% select(c("SampleID", "Microbiota"))
intermediate<- (merge(data, intermediate, by = 'SampleID'))
data<- intermediate
data$Microbiota
data <- data %>% filter(Microbiota=="Cedars_SPF")

spf_gavage_otus_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, observed_otus, Type, 0, 600)+ 
  facet_grid(~Site) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("SPF Gavage")+
  labs(y="# ASVs", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)
spf_gavage_pe_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, pielou_e, Type, 0, 1)+ 
  facet_grid(~Site) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)


# HUM Gavage 
data<-read.csv("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv", header=TRUE, row.names=1)
metadata<- read.delim("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv")
intermediate<- metadata %>% select(c("SampleID", "Microbiota"))
intermediate<- (merge(data, intermediate, by = 'SampleID'))
data<- intermediate
data$Microbiota
data <- data %>% filter(Microbiota=="Humanized")

hum_otus_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, observed_otus, Type, 0, 600)+ 
  facet_grid(~Site) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM SD Gavage")+
  labs(y="# ASVs", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)
hum_pe_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, pielou_e, Type, 0, 1)+ 
  facet_grid(~Site) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF")+
  labs(y="Pielou's evenness", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("L", "M")),
                             method="wilcox",vjust=0.5,
                             label="p.signif", hide.ns = TRUE)

### Beta Diversity ---
Type_cols<-c("L"="#481567FF", "M" = "#3CBB75FF")

# HUM V Gavage
data<-read.csv(here("Donors-Analysis/type_rpca/Donors RPCA - Colon.csv"),header=FALSE)
metadata<- read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"), header=TRUE)

hum_v_gavage_pcoa_colon <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"SPF Gavage", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5)) 
hum_v_gavage_pcoa_colon

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - SI.csv",
               header=FALSE)

hum_v_gavage_pcoa_si <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"SPF Gavage", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5)) 
hum_v_gavage_pcoa_si

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - Duodenum.csv",
               header=FALSE)

hum_v_gavage_pcoa_duo <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
hum_v_gavage_pcoa_duo

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - Jejunum.csv",
               header=FALSE)

hum_v_gavage_pcoa_jej <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"HUM V. Gavage", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
hum_v_gavage_pcoa_jej

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - Ileum.csv",
               header=FALSE)

hum_v_gavage_pcoa_ile <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
hum_v_gavage_pcoa_ile

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - Cecum.csv",
               header=FALSE)

hum_v_gavage_pcoa_cec <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
hum_v_gavage_pcoa_cec 

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - PC.csv",
               header=FALSE)

hum_v_gavage_pcoa_pc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
hum_v_gavage_pcoa_pc

data<-read.csv("Donors-Analysis/type_rpca/Donors RPCA - DC .csv",
               header=FALSE)

hum_v_gavage_pcoa_dc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
  hum_v_gavage_pcoa_dc


# SPF Gavage 
data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Colon.csv",header=FALSE)
metadata<- read.delim("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv", header=TRUE)

spf_gavage_pcoa_colon <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"SPF Gavage", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5)) 
spf_gavage_pcoa_colon

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - SI.csv",
               header=FALSE)

spf_gavage_pcoa_si <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"SPF Gavage", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5)) 
spf_gavage_pcoa_si

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Duodenum.csv",
               header=FALSE)

spf_gavage_pcoa_duo <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_duo

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Jejunum.csv",
               header=FALSE)

spf_gavage_pcoa_jej <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"SPF Gavage", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_jej

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Ileum.csv",
               header=FALSE)

spf_gavage_pcoa_ile <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_ile

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Cecum.csv",
               header=FALSE)

spf_gavage_pcoa_cec <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_cec

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Proximal_Colon.csv",
               header=FALSE)

spf_gavage_pcoa_pc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_pc

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Distal_Colon.csv",
               header=FALSE)

spf_gavage_pcoa_dc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
spf_gavage_pcoa_dc

# HUM Gavage

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Colon.csv",header=FALSE)
metadata<- read.delim("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv", header=TRUE)

hum_gavage_pcoa_colon <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"HUM Gavage", "Type", Type_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - SI.csv",
               header=FALSE)
hum_gavage_pcoa_si <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Type", Type_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Duo.csv",
               header=FALSE)
hum_gavage_pcoa_duo <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Jej.csv",
               header=FALSE)
hum_gavage_pcoa_jej <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Ile.csv",
               header=FALSE)
hum_gavage_pcoa_ile <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Cec.csv",
               header=FALSE)
hum_gavage_pcoa_cec <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - PC.csv",
               header=FALSE)
hum_gavage_pcoa_pc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - DC.csv",
               header=FALSE)
hum_gavage_pcoa_dc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

### Figure Assembly ---
dev.new()
plot_grid(spf_gavage_otus_transverse, hum_otus_transverse, hum_v_gavage_otus_transverse,
          spf_gavage_pe_transverse, hum_pe_transverse, hum_v_gavage_pe_transverse,
          ncol=3, labels=c("A", "B", "C", "", "",""))
dev.new()
plot_grid(spf_gavage_pcoa_si, spf_gavage_pcoa_colon, hum_gavage_pcoa_si, hum_gavage_pcoa_colon, hum_v_gavage_pcoa_si, hum_v_gavage_pcoa_colon,
          spf_gavage_pcoa_duo, spf_gavage_pcoa_cec, hum_gavage_pcoa_duo, hum_gavage_pcoa_cec,hum_v_gavage_pcoa_duo, hum_v_gavage_pcoa_cec,
          nrow=2, ncol=6, labels=c("D", "E", "F", "G", 
                                   "H", "I", "J", "K","L","M","N","O"))

dev.new()
plot_grid(spf_gavage_pcoa_jej, spf_gavage_pcoa_pc, hum_gavage_pcoa_jej, hum_gavage_pcoa_pc,hum_v_gavage_pcoa_jej, hum_v_gavage_pcoa_pc,
                    spf_gavage_pcoa_ile, spf_gavage_pcoa_dc, hum_gavage_pcoa_ile, hum_gavage_pcoa_dc,hum_v_gavage_pcoa_ile, hum_v_gavage_pcoa_dc,
                    nrow=2, ncol=6)
