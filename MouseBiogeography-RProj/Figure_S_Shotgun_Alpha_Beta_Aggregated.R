###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 10.20.2022
### Figure Number: 1, allegedly
### Figure Contents: Luminal alpha, beta, and DAT aggregated across datasets
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)
library(nlme)
library(here)
library(Microbiome.Biogeography)


setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")


### Alpha Diversity ---
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Figure_S_Shotgun_Alpha_Beta_Aggregated.R")

compare_vector <- c("DC", "Jej")

# UCLA Original
chao<- "Shotgun/alpha_diversity/alpha_min_1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv"
otu <- "Shotgun/alpha_diversity/alpha_min_1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv"
pe <- "Shotgun/alpha_diversity/alpha_min_1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv"
shannon <- "Shotgun/alpha_diversity/alpha_min_1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv"
meta <- "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"

adiv_data <- merge_adiv_shotgun(chao1_filepath = chao,
                                otus_filepath = otu,
                                pielou_filepath = pe,
                                shannon_filepath = shannon,
                                metadata_filepath = meta)

output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)

ucla_o_otus_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                               otus_filepath = otu,
                                               pielou_filepath = pe,
                                               shannon_filepath = shannon,
                                               metadata_filepath = meta, 
                                               X=Site, Y=observed_features, 
                                               fillvariable = Site, 0, 300) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA O. SPF")+
  labs(y="# Species", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)

ucla_o_otus_lum


ucla_o_pe_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                             otus_filepath = otu,
                                             pielou_filepath = pe,
                                             shannon_filepath = shannon,
                                             metadata_filepath = meta,
                                             X=Site, Y=pielou_evenness, 
                                             fillvariable = Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("UCLA O. SPF")+
  labs(y="Pielou's evenness", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
ucla_o_pe_lum


# CS SPF 
chao<- "Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv"
otu <- "Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv"
pe <- "Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv"
shannon <- "Shotgun/alpha_diversity/alpha_min_500000_CS_SPF_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv"
meta <- "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"

cs_otus_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                         otus_filepath = otu,
                                         pielou_filepath = pe,
                                         shannon_filepath = shannon,
                                         metadata_filepath = meta,
                                         X=Site, Y=observed_features, 
                                         fillvariable = Site, 0, 300) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="# Species", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
cs_otus_lum

cs_pe_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                           otus_filepath = otu,
                                           pielou_filepath = pe,
                                           shannon_filepath = shannon,
                                           metadata_filepath = meta,
                                           X=Site, Y=pielou_evenness, 
                                           fillvariable = Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("CS SPF")+
  labs(y="Pielou's evenness", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
cs_pe_lum

adiv_data <- merge_adiv_shotgun(chao1_filepath = chao,
                                otus_filepath = otu,
                                pielou_filepath = pe,
                                shannon_filepath = shannon,
                                metadata_filepath = meta)

output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)

#Hum Gavage
chao<- "Shotgun/alpha_diversity/alpha_min_50000_HUM_Gavage_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv"
otu <- "Shotgun/alpha_diversity/alpha_min_50000_HUM_Gavage_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv"
pe <- "Shotgun/alpha_diversity/alpha_min_50000_HUM_Gavage_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv"
shannon <- "Shotgun/alpha_diversity/alpha_min_50000_HUM_Gavage_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv"
meta <- "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"

hum_gavage_otus_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                           otus_filepath = otu,
                                           pielou_filepath = pe,
                                           shannon_filepath = shannon,
                                           metadata_filepath = meta,
                                           X=Site, Y=observed_features, 
                                           fillvariable = Site, 0, 300) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM Gavage")+
  labs(y="# Species", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
hum_gavage_otus_lum

hum_gavage_pe_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                         otus_filepath = otu,
                                         pielou_filepath = pe,
                                         shannon_filepath = shannon,
                                         metadata_filepath = meta,
                                         X=Site, Y=pielou_evenness, 
                                         fillvariable = Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("HUM Gavage")+
  labs(y="Pielou's evenness", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
hum_gavage_pe_lum

adiv_data <- merge_adiv_shotgun(chao1_filepath = chao,
                                otus_filepath = otu,
                                pielou_filepath = pe,
                                shannon_filepath = shannon,
                                metadata_filepath = meta)

output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)

#SPF Gavage
chao<- "Shotgun/alpha_diversity/alpha_min_400000_SPF_Gavage_BioGeo_Shotgun_ASV/chao1_dir/alpha-diversity.tsv"
otu <- "Shotgun/alpha_diversity/alpha_min_400000_SPF_Gavage_BioGeo_Shotgun_ASV/otus_dir/alpha-diversity.tsv"
pe <- "Shotgun/alpha_diversity/alpha_min_400000_SPF_Gavage_BioGeo_Shotgun_ASV/pielou_e_dir/alpha-diversity.tsv"
shannon <- "Shotgun/alpha_diversity/alpha_min_400000_SPF_Gavage_BioGeo_Shotgun_ASV/shannon_dir/alpha-diversity.tsv"
meta <- "Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"

spf_gavage_otus_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                                   otus_filepath = otu,
                                                   pielou_filepath = pe,
                                                   shannon_filepath = shannon,
                                                   metadata_filepath = meta,
                                                   X=Site, Y=observed_features, 
                                                   fillvariable = Site, 0, 300) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("SPF Gavage")+
  labs(y="# Species", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
spf_gavage_otus_lum

spf_gavage_pe_lum <- generate_adiv_plots_shotgun(chao1_filepath = chao,
                                                 otus_filepath = otu,
                                                 pielou_filepath = pe,
                                                 shannon_filepath = shannon,
                                                 metadata_filepath = meta,
                                                 X=Site, Y=pielou_evenness, 
                                                 fillvariable = Site, 0, 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  ggtitle("SPF Gavage")+
  labs(y="Pielou's evenness", x="")+
  stat_compare_means(comparisons = compare_vector,
                     method="wilcox", vjust=0.5,label="p.signif", hide.ns = TRUE)
spf_gavage_pe_lum

adiv_data <- merge_adiv_shotgun(chao1_filepath = chao,
                                otus_filepath = otu,
                                pielou_filepath = pe,
                                shannon_filepath = shannon,
                                metadata_filepath = meta)

output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)
output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Site, random = ~1|MouseID, data=adiv_data)
summary(output)

## Final Figures --

alpha_diversity_lum <- plot_grid(ucla_o_otus_lum, cs_otus_lum, spf_gavage_otus_lum, hum_gavage_otus_lum, nrow=1)
alpha_diversity_pe_lum <- plot_grid(ucla_o_pe_lum, cs_pe_lum, spf_gavage_pe_lum, hum_gavage_pe_lum, nrow=1)
### Beta Diversity ---
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Site_cols <- c("Jej"="gold", "DC" = "magenta")


# UCLA Original SPF
data<-readr::read_csv(here("Shotgun/Site_RPCA/rpca/BioGeo_RPCA - UCLA_O_SPF_ordination.csv"))
data <- rbind(names(data), data)
metadata <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t")

ucla_o_pcoa_lum <- generate_pcoA_plots_shotgun(data,metadata,"UCLA O. SPF", "Site", Site_cols)+
  labs(title="UCLA O. SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
ucla_o_pcoa_lum

# PERMANOVA
data <- read.delim("Shotgun/Site_RPCA/rpca/dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv",row.names = 1)
row.names(data) <- gsub("-",".",row.names(data))
write.table(data,"Shotgun/Site_RPCA/rpca/dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv")


permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

run_repeated_PERMANOVA(path_to_distance_matrix_tsv ="Shotgun/Site_RPCA/rpca/dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv",
                       path_to_metadata_csv = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

metadata <- read.csv("Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",header=TRUE, row.names=1)
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~  Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis$aov.tab

# CS SPF
data<-readr::read_csv(here("Shotgun/Site_RPCA/rpca/BioGeo_RPCA - CS_SPF_ordination.csv"))
data <- rbind(names(data), data)
metadata <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t")

cs_spf_lum <- generate_pcoA_plots_shotgun(data,metadata,"CS SPF", "Site", Site_cols) +
  labs(title="CS SPF") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

# PERMANOVA
data <- read.delim("Shotgun/Site_RPCA/rpca/dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv",row.names = 1)
row.names(data) <- gsub("-",".",row.names(data))
write.table(data,"Shotgun/Site_RPCA/rpca/dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv")


permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

run_repeated_PERMANOVA(path_to_distance_matrix_tsv ="Shotgun/Site_RPCA/rpca/dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv",
                       path_to_metadata_csv = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

metadata <- read.csv("Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",header=TRUE, row.names=1)
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~  Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis$aov.tab

# SPF Gavage
data<-readr::read_csv(here("Shotgun/Site_RPCA/rpca/BioGeo_RPCA - SPF_Gavage_ordination.csv"))
data <- rbind(names(data), data)
metadata <- readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t")

spf_pcoa_lum <- generate_pcoA_plots_shotgun(data,metadata,"CS SPF Gavage", "Site", Site_cols) +
  labs(title="SPF Gavage") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 

# PERMANOVA
data <- read.delim("Shotgun/Site_RPCA/rpca/dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv",row.names = 1)
row.names(data) <- gsub("-",".",row.names(data))
write.table(data,"Shotgun/Site_RPCA/rpca/dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv")


permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Shotgun/Site_RPCA/rpca/dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv",
                       path_to_metadata_csv = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

 target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data.dist, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Site , data=metadata, permutations=10000)
data.adonis$aov.tab

# HUM Gavage
data<-readr::read_csv(here("Shotgun/Site_RPCA/rpca/BioGeo_RPCA - HUM_Gavage_ordination.csv"))
data <- rbind(names(data), data)
metadata <- readr::read_delim(here("Shotgun/starting_files/BioGeo_Shotgun_Metadata.tsv"),delim="\t")

hum_pcoa_lum <- generate_pcoA_plots_shotgun(data,metadata,"Hum Gavage", "Site", Site_cols) +
  labs(title="HUM Gavage") + 
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) 

# PERMANOVA
data <- read.delim("Shotgun/Site_RPCA/rpca/dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza.txt/distance-matrix.tsv",row.names = 1)
row.names(data) <- gsub("-",".",row.names(data))
write.table(data,"Shotgun/Site_RPCA/rpca/dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv")


permute_within <- c("Site")
subject_data <- c("Sequencing_Run", "Sex", "MouseID")

run_repeated_PERMANOVA(path_to_distance_matrix_tsv = "Shotgun/Site_RPCA/rpca/dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza.txt/2_distance-matrix.tsv",
                       path_to_metadata_csv = "Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",
                       permute_columns_vector = permute_within,
                       subject_metadata_vector=subject_data)

metadata <- read.csv("Shotgun/starting_files/BioGeo_Shotgun_Metadata.csv",header=TRUE, row.names=1)
target <- row.names(data)
metadata = metadata[match(target, row.names(metadata)),]
target == row.names(metadata)
data.dist <- as.dist(as(data, "matrix"))

set.seed(11)
data.adonis=adonis(data.dist ~  Sequencing_Run + Sex + Site, data=metadata, permutations=10000)
data.adonis$aov.tab
### Compile alpha and beta diversity ---
fig_luminal_top <- plot_grid(alpha_diversity_lum, alpha_diversity_pe_lum, 
                             ncol=1,
                             labels=c("A","B"))

fig_luminal_bottom <- plot_grid(ucla_o_pcoa_lum,cs_spf_lum,spf_pcoa_lum,hum_pcoa_lum, 
                                ncol=4, 
                                labels=c("C","","",""))
dev.new(width=15, height=10)
fig_luminal

