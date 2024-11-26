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
library(here)

#devtools::install_github("julianneyang/Microbiome.Biogeography")
library("Microbiome.Biogeography")



### Alpha Diversity ---

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_3_Luminal_Alpha_Beta_Diversity.R")

compare_vector <- list(c("DC", "PC"),
                       c("DC", "Cec"),
                       c("DC", "Ile"),
                       c("DC", "Jej"),
                       c("DC", "Duo"))

# HUM MD Gavage
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
  ggtitle("HUM MD Gavage")+
  labs(y="ASVs", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

hum_v_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_evenness, Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("HUM MD Gavage")+
  labs(y="Pielou's evenness", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))


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
  labs(y="ASVs", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

ucla_o_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


# CS SPF 
data <- readRDS(here("CS_SPF/alpha_min_10000_table/alpha_diversity_CS_Facility.RDS"))
metadata<- read.csv(here("CS_SPF/CS_Facility_Metadata.csv"))
type_metadata <- metadata %>% select("Type", "SampleID")
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
luminaldata<-data %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-Type)

cs_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="ASVs", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

cs_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF")+
  labs(y="Pielou's evenness", x="")+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


#HUM SD Gavage
data<-read.csv(here("Humanized-Biogeography-Analysis/alpha_diversity/alpha_diversity_Humanized.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
luminaldata <- data %>% dplyr::filter(Microbiota=="Humanized") %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-c(Microbiota,Type))

hum_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM SD Gavage")+
  labs(y="ASVs", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

hum_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Hum Gavage")+
  labs(y="Pielou's evenness", x="")+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)



#SPF Gavage
data<-read.csv(here("Humanized-Biogeography-Analysis/alpha_diversity/alpha_diversity_Humanized.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
luminaldata <- data %>% dplyr::filter(Microbiota=="Cedars_SPF") %>% dplyr::filter(Type =="Luminal")
luminaldata <- luminaldata %>% select(-c(Microbiota,Type))

spf_gavage_pe_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF Gavage")+
  labs(y="Pielou's evenness", x="")  +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

spf_gavage_otus_lum <- Microbiome.Biogeography::generate_adiv_plots(luminaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("SPF Gavage")+
  labs(y="ASVs", x="") +
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
dev.new(width=15, height=10)
alpha_diversity_lum <- plot_grid(ucla_o_otus_lum, cs_otus_lum, spf_gavage_otus_lum, hum_otus_lum,hum_v_otus_lum, nrow=1)
alpha_diversity_pe_lum <- plot_grid(ucla_o_pe_lum, cs_pe_lum, spf_gavage_pe_lum, hum_pe_lum, hum_v_pe_lum,nrow=1)

### Beta Diversity ---
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
SI_cols <- c("D" = "firebrick", "J"="gold", "I" = "forestgreen")
Colon_cols <- c("C" = "cyan", "PC" = "blue", "DC" = "magenta")


# UCLA Original SPF
data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - Luminal ordination.csv"),header=FALSE)
metadata<- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv"), header=TRUE,row.names=1)

ucla_o_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site_General", cols_general)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  labs(y="PC 2 (49.0%)")+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

ucla_o_pcoa_lum

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - LumCol ordination.csv"),
               header=FALSE)

ucla_o_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", Colon_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

ucla_o_pcoa_lc

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - LumSI ordination.csv"),
               header=FALSE)

ucla_o_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", SI_cols)+
  labs(y="PC 2 (49.0%)")+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

ucla_o_pcoa_lsi

# CS SPF

metadata <- read.csv(here("CS_SPF/CS_Facility_Metadata.csv"),
                     header=TRUE)
data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal.csv"),
               header=FALSE)
cs_spf_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))


data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_Colon.csv"),
               header=FALSE)
cs_spf_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

cs_spf_lc 

data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_SI.csv"),
               header=FALSE)
cs_spf_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  labs(y="PC 2 (49.0%)")+
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

cs_spf_lsi

# SPF Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - Lum.csv"),
               header=FALSE)
spf_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site_General", cols_general) +
  #labs(title="CS SPF Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))


data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LC.csv"),
               header=FALSE)
spf_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", Colon_cols) +
  #labs(title="CS SPF Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

spf_pcoa_lc

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LSI.csv"),
               header=FALSE)
spf_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", SI_cols) +
  labs(y="PC 2 (46.0%)")+
  #labs(title="CS SPF Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

spf_pcoa_lsi

# Hum SD Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),
                     delim="\t")
data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - Luminal.csv"),
               header=FALSE)
hum_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site_General", cols_general) +
  #labs(title="HUM Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

hum_pcoa_lum

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LC.csv"),
               header=FALSE)
hum_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="HUM Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

hum_pcoa_lc

data<-read.csv(here("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LSI.csv"),
               header=FALSE)
hum_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="HUM Gavage") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

hum_pcoa_lsi

# HUM MD Gavage

metadata <- readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),delim="\t")
data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - Luminal.csv"),
               header=FALSE)
data[,1] <- gsub("-",".",data[,1])
hum_v_pcoa_lum <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"HUM V SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))


data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - LC.csv"),
               header=FALSE)
data[,1] <- gsub("-",".",data[,1])

hum_v_pcoa_lc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))

hum_v_pcoa_lc

data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - LSI.csv"),
               header=FALSE)
data[,1] <- gsub("-",".",data[,1])

hum_v_pcoa_lsi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="CS SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
  theme(axis.text.x=element_text(color="#000000"),
        axis.text.y=element_text(color="#000000"),
        axis.ticks=element_line(color="#000000"),
        plot.title=element_text(face="plain", color="#000000"))


hum_v_pcoa_lsi

### Add R^2 and p-values to beta diversity figures ---
f_ucla_o_pcoa_lum <- ggdraw(add_sub(ucla_o_pcoa_lum, 
                                  bquote(paste(~italic(R)^2,"=0.261  ", 
                                               italic("P"),"<9.99 E-5"))))
f_cs_spf_lum <- ggdraw(add_sub(cs_spf_lum, 
                                  bquote(paste(~italic(R)^2,"=0.366  ", 
                                               italic("P"),"<9.99 E-5"))))

f_spf_pcoa_lum<- ggdraw(add_sub(spf_pcoa_lum, 
                             bquote(paste(~italic(R)^2,"=0.230  ", 
                                          italic("P"),"<9.99 E-5"))))

f_hum_pcoa_lum<- ggdraw(add_sub(hum_pcoa_lum, 
                              bquote(paste(~italic(R)^2,"=0.206  ", 
                                           italic("P"),"<9.99 E-5"))))
f_hum_v_pcoa_lum<- ggdraw(add_sub(hum_v_pcoa_lum, 
                              bquote(paste(~italic(R)^2,"=0.200  ", 
                                           italic("P"),"<9.99 E-5"))))
f_ucla_o_pcoa_lc <- ggdraw(add_sub(ucla_o_pcoa_lc, 
                            bquote(paste(~italic(R)^2,"=0.028  ", 
                                         italic("P"),"<9.99 E-5"))))

f_cs_spf_lc <- ggdraw(add_sub(cs_spf_lc, 
                            bquote(paste(~italic(R)^2,"=0.034  ", 
                                         italic("P"),"=0.078"))))

f_spf_pcoa_lc <- ggdraw(add_sub(spf_pcoa_lc, 
                            bquote(paste(~italic(R)^2,"=0.031  ", 
                                         italic("P"),"=0.082"))))

f_hum_pcoa_lc <- ggdraw(add_sub(hum_pcoa_lc, 
                              bquote(paste(~italic(R)^2,"=0.004  ", 
                                           italic("P"),"=0.939"))))

f_hum_v_pcoa_lc <- ggdraw(add_sub(hum_v_pcoa_lc, 
                              bquote(paste(~italic(R)^2,"=4.50 E-4  ", 
                                           italic("P"),"=0.548"))))

f_ucla_o_pcoa_lsi <- ggdraw(add_sub(ucla_o_pcoa_lsi, 
                                 bquote(paste(~italic(R)^2,"=0.017  ", 
                                              italic("P"),"<5.79 E-3"))))

f_cs_spf_lsi <- ggdraw(add_sub(cs_spf_lsi, 
                            bquote(paste(~italic(R)^2,"=0.218  ", 
                                         italic("P"),"=0.444"))))

f_spf_pcoa_lsi <- ggdraw(add_sub(spf_pcoa_lsi, 
                              bquote(paste(~italic(R)^2,"=0.097  ", 
                                           italic("P"),"=0.230"))))

f_hum_pcoa_lsi <- ggdraw(add_sub(hum_pcoa_lsi, 
                              bquote(paste(~italic(R)^2,"=0.026  ", 
                                           italic("P"),"=0.565"))))

f_hum_v_pcoa_lsi <- ggdraw(add_sub(hum_v_pcoa_lsi, 
                                bquote(paste(~italic(R)^2,"=0.008  ", 
                                             italic("P"),"=0.268"))))

interregional_lum <- plot_grid(f_ucla_o_pcoa_lum, f_cs_spf_lum, f_spf_pcoa_lum, f_hum_pcoa_lum,f_hum_v_pcoa_lum,nrow=1)
#interregional_lum <- interregional_lum +labs(title="Luminal") + theme(plot.title = element_text(hjust = 0.5)) 
interregional_lc <- plot_grid(f_ucla_o_pcoa_lc, f_cs_spf_lc, f_spf_pcoa_lc, f_hum_pcoa_lc,f_hum_v_pcoa_lc,nrow=1)
interregional_lsi <- plot_grid(f_ucla_o_pcoa_lsi, f_cs_spf_lsi, f_spf_pcoa_lsi, f_hum_pcoa_lsi,f_hum_v_pcoa_lsi,nrow=1)

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
