###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 10.20.2022
### Figure Number: 1, allegedly
### Figure Contents: Luminal alpha, beta, and DAT aggregated across datasets
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
library(here)

remove.packages("Microbiome.Biogeography")
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")


### Alpha Diversity ---

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Luminal_Aggregated.R")

compare_vector <- list(c("DC", "PC"),
                       c("DC", "Cec"),
                       c("DC", "Ile"),
                       c("DC", "Jej"),
                       c("DC", "Duo"))

# HUM V Gavage
data<-readr::read_rds(here("Donors-Analysis/alpha_diversity/alpha_diversity.RDS"))
metadata<- readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"))
type_metadata <- metadata %>% select("Type", "SampleID")
type_metadata$SampleID <- gsub("-",".",type_metadata$SampleID)
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
luminaldata<-data %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-Type)

hum_v_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_features, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hum V. Gavage")+
  labs(y="# ASVs", x="")

hum_v_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_evenness, Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("HUM V. Gavage")+
  labs(y="# ASVs", x="")


# UCLA Original
data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv"), header=TRUE, row.names=1)
metadata<- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv"), header=TRUE,row.names=1)
type_metadata <- metadata %>% select("Type", "SampleID")
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
luminaldata<-data %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-Type)

ucla_o_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA O. SPF")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

ucla_o_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


# CS SPF 
data <- readRDS(here("CS-Facility-Analysis/alpha_diversity_CS_Facility.RDS"))
metadata<- read.csv(here("CS-Facility-Analysis/CS_Facility_Metadata.csv"))
type_metadata <- metadata %>% select("Type", "SampleID")
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
luminaldata<-data %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-Type)

cs_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

cs_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF")+
  labs(y="Pielou's Evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


#Hum Gavage
data<-read.csv(here("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
luminaldata <- data %>% dplyr::filter(Microbiota=="Humanized") %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-c(Microbiota,Type))

hum_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("Hum Gavage")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

hum_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Hum Gavage")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)



#SPF Gavage
data<-read.csv(here("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
luminaldata <- data %>% dplyr::filter(Microbiota=="Cedars_SPF") %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-c(Microbiota,Type))

spf_gavage_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF Gavage")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

spf_gavage_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF Gavage")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
dev.new(width=15, height=10)
alpha_diversity_lum <- plot_grid(ucla_o_otus_lum, cs_otus_lum, spf_gavage_otus_lum, hum_otus_lum,hum_v_otus_lum, nrow=1)
alpha_diversity_pe_lum <- plot_grid(ucla_o_pe_lum, cs_pe_lum, spf_gavage_pe_lum, hum_pe_lum, hum_v_pe_lum,nrow=1)

### Beta Diversity ---
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
SI_cols <- c("Duo" = "firebrick", "Jej"="gold", "Ile" = "forestgreen")
Colon_cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")


# UCLA Original SPF
data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - Luminal ordination.csv"),header=FALSE)
metadata<- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv"), header=TRUE,row.names=1)

ucla_o_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site_General", cols_general)+
  #labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_lum

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - LumCol ordination.csv"),
               header=FALSE)

ucla_o_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", Colon_cols)+
  #labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_lc

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - LumSI ordination.csv"),
               header=FALSE)

ucla_o_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", SI_cols)+
  #labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_lsi

# CS SPF

metadata <- read.csv(here("CS-Facility-Analysis/CS_Facility_Metadata.csv"),
                     header=TRUE)
data<-read.csv(here("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal.csv"),
               header=FALSE)
cs_spf_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))

data<-read.csv(here("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_Colon.csv"),
               header=FALSE)
cs_spf_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
cs_spf_lc 

data<-read.csv(here("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_SI.csv"),
               header=FALSE)
cs_spf_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))

cs_spf_lsi

# SPF Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - Lum.csv"),
               header=FALSE)
spf_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site_General", cols_general) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LC.csv"),
               header=FALSE)
spf_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", Colon_cols) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
spf_pcoa_lc

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LSI.csv"),
               header=FALSE)
spf_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", SI_cols) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
spf_pcoa_lsi

# Hum Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),
                     delim="\t")
data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - Luminal.csv"),
               header=FALSE)
hum_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site_General", cols_general) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_lum

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LC.csv"),
               header=FALSE)
hum_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_lc

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LSI.csv"),
               header=FALSE)
hum_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_lsi

# HUM V Gavage

metadata <- readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),delim="\t")
data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - Luminal.csv"),
               header=FALSE)
hum_v_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"HUM V SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - LC.csv"),
               header=FALSE)
hum_v_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
hum_v_pcoa_lc

data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - LSI.csv"),
               header=FALSE)
hum_v_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))

hum_v_pcoa_lsi

interregional_lum <- plot_grid(ucla_o_pcoa_lum, cs_spf_lum, spf_pcoa_lum, hum_pcoa_lum,hum_v_pcoa_lum,nrow=1)
#interregional_lum <- interregional_lum +labs(title="Luminal") + theme(plot.title = element_text(hjust = 0.5)) 
interregional_lc <- plot_grid(ucla_o_pcoa_lc, cs_spf_lc, spf_pcoa_lc, hum_pcoa_lc,hum_v_pcoa_lc,nrow=1)
interregional_lsi <- plot_grid(ucla_o_pcoa_lsi, cs_spf_lsi, spf_pcoa_lsi, hum_pcoa_lsi,hum_v_pcoa_lsi,nrow=1)

### Compile alpha and beta diversity ---
fig_luminal_top <- plot_grid(alpha_diversity_lum, alpha_diversity_pe_lum, 
                             ncol=1,
                             labels=c("A","B"))
dev.new(width=15, height=10)
fig_luminal_top

fig_luminal_bottom <- plot_grid(interregional_lum,
                         interregional_lc, interregional_lsi, ncol=1, nrow=3,
                         labels=c("C","D","E"))
dev.new(width=15, height=10)
fig_luminal_bottom