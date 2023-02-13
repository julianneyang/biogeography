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

### Alpha Diversity ---
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Figure_Luminal_aggregated.R")

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv", header=TRUE, row.names=1)
metadata<- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv", header=TRUE,row.names=1)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate
luminaldata<-data %>% dplyr::filter(Type =="Luminal")

ucla_o_otus <- Microbiome.Biogeography::generate_adiv_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv",
                                   "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv",
                                   Site, observed_otus, Site, 0, 725) +
  facet_wrap(~Type) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
ucla_o_otus
