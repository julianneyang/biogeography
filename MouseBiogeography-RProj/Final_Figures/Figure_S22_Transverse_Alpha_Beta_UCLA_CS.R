###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 10.20.2022
### Figure Number: 
### Figure Contents: UCLA O SPF and CS SPF Transverse alpha and beta diversity 
###### whining ends here ---
library(renv)
renv::snapshot()
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(Microbiome.Biogeography)


remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")


### Alpha Diversity ---
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Transverse_Alpha_Beta_UCLA_CS.R")

compare_vector <- list(c("DC", "PC"),
                       c("DC", "Cec"),
                       c("DC", "Ile"),
                       c("DC", "Jej"),
                       c("DC", "Duo"))

# UCLA Original
data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv", header=TRUE, row.names=1)
metadata<- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv", header=TRUE,row.names=1)

ucla_o_otus_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, observed_otus, Type, 0, 600)+ 
  facet_grid(~Site) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA O. SPF")+
  labs(y="# ASVs", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("Lum", "Muc")),
                     method="wilcox",vjust=0.5,
                     label="p.signif",step.increase=0.08, hide.ns = TRUE)
ucla_o_pe_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, pielou_e, Type, 0, 1)+ 
  facet_grid(~Site) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("Lum", "Muc")),
                             method="wilcox",vjust=0.5,
                             label="p.signif",step.increase=0.08, hide.ns = TRUE)


# CS SPF
data <- readRDS("CS-Facility-Analysis/alpha_diversity_CS_Facility.RDS")
metadata<- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv")

cs_otus_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, observed_otus, Type, 0, 600)+ 
  facet_grid(~Site) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="# ASVs", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("Lum", "Muc")),
                             method="wilcox",vjust=0.5,
                             label="p.signif",step.increase=0.08, hide.ns = TRUE)
cs_pe_transverse <- Microbiome.Biogeography::generate_adiv_plots(data, metadata,Type, pielou_e, Type, 0, 1)+ 
  facet_grid(~Site) +
  #theme(plot.title = element_text(hjust = 0.5)) +
  #ggtitle("CS SPF")+
  labs(y="Pielou's evenness", x="")+
  ggpubr::stat_compare_means(comparisons = list(c("Lum", "Muc")),
                             method="wilcox",vjust=0.5,
                             label="p.signif",step.increase=0.08, hide.ns = TRUE)

### Beta Diversity ---
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")


# UCLA Original SPF
data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Colon ordination.csv",header=FALSE)
metadata<- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv", header=TRUE,row.names=1)

ucla_o_pcoa_colon <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
  #theme(plot.title = element_text(hjust = 0.5)) 
ucla_o_pcoa_colon

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - SI_ordination.csv",
               header=FALSE)

ucla_o_pcoa_si <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="Colon") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5)) 
ucla_o_pcoa_si

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Duodenum ordination.csv",
               header=FALSE)

ucla_o_pcoa_duo <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_duo

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Jejunum ordination.csv",
               header=FALSE)

ucla_o_pcoa_jej <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_jej

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Ileum ordination.csv",
               header=FALSE)

ucla_o_pcoa_ile <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_ile

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Cecum ordination.csv",
               header=FALSE)

ucla_o_pcoa_cec <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_cec

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Proximal_Colon ordination.csv",
               header=FALSE)

ucla_o_pcoa_pc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_pc

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/beta_diversity/Luminal vs Mucosal_ RPCA for all Sites - Distal_Colon ordination.csv",
               header=FALSE)

ucla_o_pcoa_dc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_dc

# CS SPF

metadata <- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv",
                     header=TRUE)
data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Colon.csv",
               header=FALSE)
cs_spf_colon <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Type", Type_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
  #theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - SI.csv",
               header=FALSE)
cs_spf_si <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"CS SPF", "Type", Type_cols) +
  #labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site_General)
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Duodenum.csv",
               header=FALSE)
cs_spf_duo <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Jejunum.csv",
               header=FALSE)
cs_spf_jej <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Ileum.csv",
               header=FALSE)
cs_spf_ile <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Cecum.csv",
               header=FALSE)
cs_spf_cec <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Proximal_Colon.csv",
               header=FALSE)
cs_spf_pc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

data<-read.csv("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Distal_Colon.csv",
               header=FALSE)
cs_spf_dc <- Microbiome.Biogeography::generate_pcoA_plots(data,metadata,"UCLA O. SPF", "Type", Type_cols)+
  #labs(title="UCLA O. SPF") + 
  #theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  facet_grid(~Site)+
  theme(legend.position = "none")
#theme(plot.title = element_text(hjust = 0.5))

### Figure Assembly ---
plot_grid(ucla_o_otus_transverse, cs_otus_transverse, 
          ucla_o_pe_transverse, cs_pe_transverse, 
          ncol=2, labels=c("A", "", "B", ""))
plot_grid(ucla_o_pcoa_si, ucla_o_pcoa_colon, cs_spf_si, cs_spf_colon,
          ucla_o_pcoa_duo, ucla_o_pcoa_cec, cs_spf_duo, cs_spf_cec,
          ucla_o_pcoa_jej, ucla_o_pcoa_pc, cs_spf_jej, cs_spf_pc,
          ucla_o_pcoa_ile, ucla_o_pcoa_dc, cs_spf_ile, cs_spf_dc,
          nrow=4, ncol=4, labels=c("C", "D", "E", "F", 
                                   "G", "H", "I", "J"))
