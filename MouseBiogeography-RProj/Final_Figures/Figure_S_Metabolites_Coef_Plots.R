###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 3.10.2023
### Figure Number: 6
### Figure Contents: GMM Site Heatmaps for all Cohorts 
###### whining ends here ---

library(ggplot2)
library(tidyverse)
#library(rlang)
library(cowplot)
library(viridis)
#library(plyr)
library(gridExtra)
library(paletteer)
library(ComplexUpset)
library(here)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Coef_Plots.R")

### Upset Plot ---

file_paths <- c("Regional-Mouse-Biogeography-Analysis/melonnpan/Mets-SeqRunLineSexSite_General_1-MouseID/all_results.tsv",
                "CS_SPF/melonnpan/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv",
                "Donors-Analysis/melonnpan/Mets-SeqRunSexSite_General-1-MsID-DonorID/all_results.tsv",
                "Humanized-Biogeography-Analysis/melonnpan/HUM_SD_Gavage/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv",
                "Humanized-Biogeography-Analysis/melonnpan/SPF_Gavage/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv")


cohort_prefixes <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_MD_Gavage",
                     "HUM_SD_Gavage",
                     "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site_General")


id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)
#id_df_wide <- id_df_wide %>% mutate(SPF_Gavage = 0)

all_taxa <- all_taxa %>% select(c("feature", "Cohort")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)
df_wide <- as.data.frame(df_wide)
#df_wide <- df_wide %>% mutate(SPF_Gavage = 0)

df_wide <- as.data.frame(df_wide)
all_datasets <- names(df_wide)[-(1)]
mets_upset <- upset(df_wide, all_datasets,width_ratio=0.2,
                                  base_annotations=list(
                                    'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
                                      scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')),
                   #min_degree=3,
                   themes=upset_default_themes(text=element_text(size=16), colour="black"))
                   #themes=upset_modify_themes(
                     #list(
                       #'intersections_matrix'=theme(text=element_text(size=16,colour = "black"))
                     #)
                   #))


id_df_wide$count_ones <- rowSums(id_df_wide %>% select(all_of(cohort_prefixes)))
df_filtered <- id_df_wide[id_df_wide$count_ones >= 5, ]
length(df_filtered$feature)
length(unique(df_filtered$feature))
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
metabs_of_interest <- df_filtered$feature
write_rds(metabs_of_interest,here("Highlighted_metabolites.RDS"))

test <- read.delim(here("Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/all_results.tsv"))

### Make a big barplot for shared_metabs ---
gmm <- readRDS(here("Highlighted_GMM_Fig_6.RDS"))
names(gmm)
gbm <- readRDS(here("Highlighted_Luminal_GBM.RDS"))
names(gbm)
mets <- readRDS(here("Highlighted_metabolites.RDS"))
names(mets) <-  c(
  "Lipids",                 # "LPC.O.16.0"
  "Lipids",                 # "D.erythro.N.stearoylsphingosine"
  "Small Organic Compounds",# "phosphate"
  "Amino Acids",            # "Thr.Ala"
  "Lipids",                 # "Cer.38.3.O4.Cer.24.2.O3.14.1.2OH."
  "Lipids",                 # "phosphoinositol"
  "Lipids",                 # "Cer.36.1.O3.Cer.18.0.O2.18.1"
  "Amino Acids",            # "Ala.Gln"
  "Lipids",                 # "Cer.36.3.O2.Cer.18.1.O2.18.2"
  "Amino Acids",            # "Phe.Thr"
  "Lipids",                 # "Cer.17.1.2O.24.1"
  "Amino Acids",            # "Ala.Ala.Lys"
  "Amino Acids",            # "Ala.Thr"
  "Small Organic Compounds",# "Cholesterol"
  "Amino Acids",            # "Leu.Trp"
  "Amino Acids",            # "Val.Ser"
  "Amino Acids",            # "Ile.Ala"
  "Lipids",                 # "X1..1Z.Hexadecenyl..sn.glycero.3.phosphocholine"
  "Amino Acids",            # "Phe.Val"
  "Amino Acids",            # "Ile.Gln"
  "Small Organic Compounds",# "X3..1.Pyrazolyl..alanine"
  "Small Organic Compounds",# "X3.Methyl_gamma.butyrolactone"
  "Small Organic Compounds",# "piperidone"
  "Amino Acids",            # "Ile.Val.Arg"
  "Small Organic Compounds",# "X2.Hydroxyphenethyl.alcohol"
  "Lipids",                 # "Cer.d44.2"
  "Lipids",                 # "Cer.42.2.O2.Cer.18.1.O2.24.1"
  "Lipids",                 # "Cer.38.2.O2.Cer.18.1.O2.20.1"
  "Lipids",                 # "Cer.38.5.O3.Cer.18.1.O2.20.4.O"
  "Lipids",                 # "SPB.18.1.O2"
  "Small Organic Compounds",# "X5.aminovaleric.acid"
  "Small Organic Compounds",# "glyceric.acid"
  "Lipids",                 # "Cer.d40.2"
  "Small Organic Compounds",# "lactulose"
  "Small Organic Compounds",# "glycerol.3.galactoside"
  "Lipids",                 # "HexCer.42.2.O4.HexCer.18.1.O3.24.1.2OH."
  "Amino Acids",            # "Isoleucine"
  "Amino Acids",            # "Tryptophan"
  "Small Organic Compounds",# "Deoxycholate"
  "Small Organic Compounds",# "guanidinosuccinate"
  "Lipids",                 # "Cer.NDS.d36.0"
  "Lipids",                 # "Cer.36.1.O2.Cer.18.0.O2.18.1"
  "Lipids"                  # "Cer.38.0.O3.Cer.20.0.O2.18.0.O"
)


# Load all the tsv files, add Cohort column, and merge them
all_data <- lapply(1:length(file_paths), function(i) {
  df <- read_tsv(file_paths[i])
  df <- df %>% mutate(Cohort = cohort_prefixes[i])
  return(df)
}) %>% bind_rows()


# Filter the dataframe based on values in "feature" column matching `mets`
filtered_data <- all_data %>% filter(feature %in% mets)%>% 
  filter(metadata == "Site_General") %>%
  filter(qval<0.05)

# Create new column "Site" based on "coef"
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "Colon", "SI"))
filtered_data <- filtered_data %>% 
  mutate(data_type = "16S")
# Create horizontal barplot with faceting by "Cohort"
cols = viridis(2)
metab_plot <- ggplot(filtered_data, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Region-Specific Metabolites") +
  theme_cowplot(16)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))


dev.new(width=10, height=20)
metab_plot

### add shotgun data to big plot ---

shotgun_fp <- c("Shotgun/melonnpan/UCLA_O_SPF/Mets-LineSexSite_General-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/CS_SPF/Mets-SexSite-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/HUM_SD_Gavage/Mets-SexSite-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/all_results.tsv")

shotgun_prefix <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_SD_Gavage",
                     "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site")

# Given 0 significant metabolites pick out the features mapping to mets
all_shotgun_data <- lapply(1:length(shotgun_fp), function(i) {
  df <- read_tsv(shotgun_fp[i])
  df <- df %>% mutate(Cohort = shotgun_prefix[i])
  return(df)
}) %>% bind_rows()

# Filter the dataframe based on values in "feature" column matching `mets`
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% mets)%>% 
  filter(metadata == "Site") 
filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Distal_Colon", "Jejunum"))
filtered_shotgun_data <-filtered_shotgun_data %>% 
  mutate(data_type = "Shotgun")
# Bind the dataframe together
full_df <- rbind(filtered_shotgun_data, filtered_data)
full_df$annotation <- names(mets)[match(full_df$feature, mets)]

aa <- full_df %>% filter(annotation=="Amino Acids")
lipids <- full_df %>% filter(annotation=="Lipids")
org <- full_df %>% filter(annotation=="Small Organic Compounds")

aa_plot <- ggplot(aa, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Amino Acid Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

lipids_plot <- ggplot(lipids, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Lipid Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

soc <- ggplot(org, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Other Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

dev.new(width=10,height=10)

top <- plot_grid(lipids_plot,
                 soc, labels=c("A","B"))
bottom <- plot_grid(aa_plot,
                 mets_upset, labels=c("C","D"), 
                 rel_widths = c(0.9,1))

dev.new(width=10,height=10)
plot_grid(top, bottom, nrow=2)

cols <- viridis::viridis(4)
metab_plot <- ggplot(full_df, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Region-Specific Metabolites") +
  theme_cowplot(16)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))


dev.new(width=10,height=10)
metab_plot

dev.new(width=10,height=10)
metab_plot + facet_wrap(~Cohort +data_type,nrow=1)
