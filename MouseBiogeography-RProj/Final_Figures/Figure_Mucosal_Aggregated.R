###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 10.20.2022
### Figure Number: 1, allegedly
### Figure Contents: Mucosal alpha, beta, and DAT aggregated across datasets
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

### Alpha Diversity ---
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Mucosal_Aggregated.R")

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
mucosaldata<-data %>% dplyr::filter(Type =="Mucosal")
mucosaldata <- mucosaldata %>% select(-Type)

hum_v_otus_muc <- Microbiome.Biogeography::generate_adiv_plots(mucosaldata, metadata,Site, observed_features, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM MD Gavage")+
  labs(y="# ASVs", x="")

hum_v_pe_muc <- Microbiome.Biogeography::generate_adiv_plots(mucosaldata, metadata,Site, pielou_evenness, Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("HUM V. Gavage")+
  labs(y="# ASVs", x="")



# UCLA Original
data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv"), header=TRUE, row.names=1)
metadata<- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv"), header=TRUE,row.names=1)
type_metadata <- metadata %>% select("Type", "SampleID")
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
mucosaldata<-data %>% dplyr::filter(Type =="Mucosal")
mucosaldata <- mucosaldata %>% select(-Type)

ucla_o_otus_muc <- Microbiome.Biogeography::generate_adiv_plots(mucosaldata, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA O. SPF")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
ucla_o_otus_muc

ucla_o_pe_muc <- Microbiome.Biogeography::generate_adiv_plots(mucosaldata, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA V. SPF")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
ucla_o_pe_muc

# CS SPF 
data <- readr::read_rds(here("CS_SPF/alpha_diversity_CS_Facility.RDS"))
metadata<- readr::read_csv(here("CS_SPF/CS_Facility_Metadata.csv"))
type_metadata <- metadata %>% select("Type", "SampleID")
intermediate<- (merge(data, type_metadata, by = 'SampleID'))
data<- intermediate
data<-data %>% dplyr::filter(Type =="Mucosal")
data <- data %>% select(-Type)

cs_otus_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
cs_otus_muc

cs_pe_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
cs_pe_muc


#Hum Gavage
data<-read.csv(here("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
data <- data %>% dplyr::filter(Microbiota=="Humanized") %>% dplyr::filter(Type =="Mucosal")
data <- data %>% select(-c(Microbiota,Type))

hum_otus_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM SD Gavage")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
hum_otus_muc 

hum_pe_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, pielou_e, Site, 0,1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("Hum Gavage")+
  labs(y="Pielou's evenness", x="")
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
hum_pe_muc 

#UCLA V. SPF 
data<- read.csv(here("UCLA_V_SPF_Analysis/alpha_diversity_WTCohort.csv"), header=TRUE, row.names=1)
metadata<- readr::read_delim(here("UCLA_V_SPF_Analysis/starting_files/UCLA_V_SPF_Metadata.tsv"),delim="\t")

ucla_v_otus_muc<- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA V. SPF")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
ucla_v_pe_muc<- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA V. SPF")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#SPF Gavage
data<-read.csv("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv", header=TRUE, row.names=1)
metadata<- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),delim="\t")
source_metadata <- metadata %>% select(c("SampleID","Type", "Microbiota"))
intermediate<- (merge(data, source_metadata, by = 'SampleID'))
data<- intermediate
data <- data %>% dplyr::filter(Microbiota=="Cedars_SPF") %>% dplyr::filter(Type =="Mucosal")
data <- data %>% select(-c(Microbiota,Type))

spf_gavage_otus_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, observed_otus, Site, 0, 600) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("SPF Gavage")+
  labs(y="# ASVs", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
spf_gavage_otus_muc

spf_gavage_pe_muc <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Site, pielou_e, Site, 0, 1) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF Gavage")+
  labs(y="Pielou's evenness", x="")
  #stat_compare_means(comparisons = compare_vector,
                     #method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
spf_gavage_pe_muc

# Compile aggregated figure
alpha_diversity_muc <- plot_grid(ucla_o_otus_muc, ucla_v_otus_muc, cs_otus_muc, spf_gavage_otus_muc, hum_otus_muc, hum_v_otus_muc, nrow=1)
#alpha_diversity_muc <- alpha_diversity_muc +labs(title="Mucosal") + theme(plot.title = element_text(hjust = 0.5)) 

alpha_diversity_pe_muc <- plot_grid(ucla_o_pe_muc, ucla_v_pe_muc, cs_pe_muc, spf_gavage_pe_muc, hum_pe_muc, hum_v_pe_muc,nrow=1)
#alpha_diversity_pe_muc <- alpha_diversity_pe_muc +labs(title="Mucosal") + theme(plot.title = element_text(hjust = 0.5)) 

### Beta Diversity ---
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
SI_cols <- c("Duo" = "firebrick", "Jej"="gold", "Ile" = "forestgreen")
Colon_cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")


# UCLA Original SPF
metadata<- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv"), header=TRUE,row.names=1)

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - Mucosal ordination.csv"),
               header=FALSE)

ucla_o_pcoa_muc<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site_General", cols_general)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #labs(title="UCLA O. SPF") + 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_muc

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - MucCol ordination.csv"),
               header=FALSE)
ucla_o_pcoa_mc<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", Colon_cols)+
  #labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_mc

data<-read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/RPCA for all Sites - MucSI ordination.csv"),
               header=FALSE)
ucla_o_pcoa_msi<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", SI_cols)+
  #labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_msi

# UCLA Validation SPF
metadata<- read.delim(here("UCLA_V_SPF_Analysis/starting_files/UCLA_V_SPF_Metadata.tsv"), header=TRUE)

data<-read.csv(here("UCLA_V_SPF_Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal.csv"),
               header=FALSE)

ucla_v_pcoa_muc<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA V. SPF", "Site_General", cols_general)+
  #labs(title="UCLA V. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_v_pcoa_muc

data<-read.csv(here("UCLA_V_SPF_Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal_Colon.csv"),
               header=FALSE)
ucla_v_pcoa_mc<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA V. SPF", "Site", Colon_cols)+
  #labs(title="UCLA V. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_v_pcoa_mc

data<-read.csv(here("UCLA_V_SPF_Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal_SI.csv"),
               header=FALSE)
ucla_v_pcoa_msi<- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Site", SI_cols)+
  #labs(title="UCLA V. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
ucla_v_pcoa_msi

# CS SPF
metadata <- read.csv(here("CS_SPF/CS_Facility_Metadata.csv"),
                     header=TRUE)
data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal.csv"),
               header=FALSE)
cs_spf_muc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
cs_spf_muc

data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_Colon.csv"),
               header=FALSE)
cs_spf_mc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
cs_spf_mc

data<-read.csv(here("CS_SPF/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_SI.csv"),
               header=FALSE)
cs_spf_msi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
cs_spf_msi

# SPF Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),
                              delim="\t")
data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - Muc.csv",
               header=FALSE)
spf_pcoa_muc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site_General", cols_general) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
spf_pcoa_muc

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MC.csv",
               header=FALSE)
spf_pcoa_mc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", Colon_cols) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
spf_pcoa_mc

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MSI.csv",
               header=FALSE)
spf_pcoa_msi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF Gavage", "Site", SI_cols) +
  #labs(title="CS SPF Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
spf_pcoa_msi

# Hum Gavage
metadata <- readr::read_delim(here("Humanized-Biogeography-Analysis/starting_files/Humanized-Metadata.tsv"),
                              delim="\t")
data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - Mucosal.csv",
               header=FALSE)
hum_pcoa_muc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"", "Site_General", cols_general) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_muc

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - MC.csv",
               header=FALSE)
hum_pcoa_mc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"", "Site", Colon_cols) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_mc

data<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - MSI.csv",
               header=FALSE)
hum_pcoa_msi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"", "Site", SI_cols) +
  #labs(title="HUM Gavage") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
  #theme(plot.title = element_text(hjust = 0.5))
hum_pcoa_msi

# HUM V Gavage
metadata <- readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),delim="\t")
data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - Mucosal.csv"),
               header=FALSE)
hum_v_pcoa_muc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"HUM V SPF", "Site_General", cols_general) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - MC.csv"),
               header=FALSE)
hum_v_pcoa_mc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", Colon_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))
hum_v_pcoa_mc

data<-read.csv(here("Donors-Analysis/site_rpca/Donors RPCA - MSI.csv"),
               header=FALSE)
hum_v_pcoa_msi <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Site", SI_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 
#theme(plot.title = element_text(hjust = 0.5))

hum_v_pcoa_msi

# Compile aggregated figure
interregional_muc <- plot_grid(ucla_o_pcoa_muc, ucla_v_pcoa_muc, cs_spf_muc, spf_pcoa_muc, hum_pcoa_muc, hum_v_pcoa_muc, nrow=1)
#interregional_muc <- interregional_muc +labs(title="Mucosal") + theme(plot.title = element_text(hjust = 0.5)) 
interregional_mc <- plot_grid(ucla_o_pcoa_mc, ucla_v_pcoa_mc, cs_spf_mc, spf_pcoa_mc, hum_pcoa_mc, hum_v_pcoa_mc,nrow=1)
#interregional_mc <- interregional_mc +labs(title="Mucosal") + theme(plot.title = element_text(hjust = 0.5)) 
interregional_msi <- plot_grid(ucla_o_pcoa_msi, ucla_v_pcoa_msi, cs_spf_msi, spf_pcoa_msi, hum_pcoa_msi, hum_v_pcoa_msi,nrow=1)
#interregional_msi <- interregional_msi +labs(title="Mucosal") + theme(plot.title = element_text(hjust = 0.5)) 

#Compile alpha and beta diversity 
fig_mucosal_top <- plot_grid(alpha_diversity_muc, alpha_diversity_pe_muc,
                             ncol=1,
                         labels=c("A","B"))
dev.new(width=15, height=10)
fig_mucosal_top

fig_mucosal_bottom <- plot_grid(interregional_muc,
                             interregional_mc, interregional_msi, 
                             ncol=1, nrow=3,
                             labels=c("C","D","E"))
dev.new(width=15, height=10)
fig_mucosal_bottom
