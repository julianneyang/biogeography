###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 5.29.23
### Figure Number: Main Figure GBM  
### Figure Contents: GBM Luminal and Mucosal Intraregional Coef plots
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(RColorBrewer)
#library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(Microbiome.Biogeography)
library(paletteer)
library(here)
library(ComplexUpset)

setwd("/home/julianne/Documents/biogeography/")
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_7_GBM_Coef_Plots.R")

### Wrangle GBM data for target features ---

process_results_for_upset_plot <- function(file_paths, cohort_prefixes) {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site" & qval<0.05)
  
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    data <- data %>% mutate(Cohort = cohort)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  
  return(data_all)
}

process_results_files <- function(file_paths, feature_value, new_value, new_coef, cohort_prefixes) {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site" & feature == feature_value)
    
    # Add a new row to the data frame
    new_row <- list(metadata = "Site",
                    feature = feature_value,
                    value = new_value,
                    coef = new_coef,
                    stderr = 0,
                    N = NA,
                    N.not.0 = NA,
                    pval = NA,
                    qval = 100)
    data <- rbind(data, new_row)
    
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    data <- data %>% mutate(Cohort = cohort)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  
  return(data_all)
}

process_gbm_files_shotgun <- function(file_paths, feature_value, cohort_prefixes,feature_annotation) {
  data_all <- NULL
  
  for (i in seq_along(file_paths)) {
    file_path <- file_paths[i]
    cohort_prefix <- cohort_prefixes[i]
    
    # Read the results file
    results <- read.table(here(file_path), header = TRUE)
    
    # Filter the results for the specified feature
    data <- filter(results, metadata == "Site" & feature == feature_value)
    
    # Add a cohort variable
    cohort <- paste0(cohort_prefix)
    annotation <- feature_annotation
    data <- data %>% mutate(Cohort = cohort)
    data <- data %>% mutate(Annotation=annotation)
    
    # Append to the combined data frame
    if (is.null(data_all)) {
      data_all <- data
    } else {
      data_all <- rbind(data_all, data)
    }
  }
  
  return(data_all)
}

### Upset Plot ---


lum_file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS_SPF/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")
  
lum_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "HUM_SD_Gavage",
                         "SPF_Gavage")

all_taxa <- process_results_for_upset_plot(file_paths = lum_file_paths,
                                           cohort_prefixes = lum_cohort_prefixes,
                                           filter_by = "Site")

module_key <- read.csv(here("GBM_Module_Key.csv"))
anno <- module_key %>% select(c("feature", "annotation"))
all_taxa <- merge(all_taxa, anno, by="feature")

### List GBM of interest ---
id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","annotation","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)
id_df_wide <- id_df_wide %>% mutate(SPF_Gavage = 0)

id_df_wide$count_ones <- rowSums(id_df_wide[, c(4:8)])
df_filtered <- id_df_wide[id_df_wide$count_ones >= 3, ]
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
gbm_of_interest <- df_filtered$feature

names(gbm_of_interest) <-df_filtered$annotation
readr::write_rds(gbm_of_interest, "Highlighted_Luminal_GBM.RDS")

### Make Upset Plot
all_taxa<- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
all_taxa <- all_taxa %>% select(c("feature","annotation","Cohort")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)
df_wide <- as.data.frame(df_wide)
df_wide <- df_wide %>% mutate(SPF_Gavage = 0)

df_wide <- as.data.frame(df_wide)
all_datasets <- names(df_wide)[-(1:2)]
gbm_upset <- ComplexUpset::upset(df_wide, all_datasets,themes=upset_default_themes(axis.title = element_text(color = 'black'),
                                                                                   axis.text=element_text(color='black'),
                                                                                   axis.ticks = element_line(color = 'black')),
                                 base_annotations=list(
                                   'Intersection size'=intersection_size(counts=TRUE,
                                                                         mapping=aes(fill='bars_color')) + 
                                     scale_fill_manual(values=c('bars_color'='skyblue'), guide='none'))) 



### Make GBM Coef plots ---


new_value <- "Distal_Colon"
new_coef <- 0
data_all <- process_results_files(lum_file_paths, feature_value, new_value, new_coef, lum_cohort_prefixes)
final_df <- data_all[FALSE,]
for (i in seq_along(gbm_of_interest)) {
  feature_value <- gbm_of_interest[i]
  data_all <- process_results_files(lum_file_paths, feature_value, new_value, new_coef, lum_cohort_prefixes)
  final_df <- rbind(final_df,data_all)
}

# color legend for coef plots 
lum_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "HUM_SD_Gavage",
                         "SPF_Gavage")
my_palette <- c(paletteer_d("basetheme::brutal",6))
names(my_palette) <-c(lum_cohort_prefixes, "UCLA_V_SPF")
cols <- my_palette[names(my_palette) %in% lum_cohort_prefixes]


map <- module_key %>% select(c("feature", "annotation"))
map_data_all <- merge(final_df,map,by="feature")

data <- map_data_all
data$value <- plyr::revalue(data$value,c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
data$value <- factor(data$value, levels = c("D", "J", "I", "C", "PC", "DC"))

# Plotting
create_plot <- function(data, anno) {
  ggplot(data %>% filter(annotation == anno),
         aes(x = value, y = coef, group = Cohort, color = Cohort)) +
    geom_line(linewidth = 2) +
    geom_hline(yintercept=0, linetype='dotted', col = 'black',linewidth = 1)+
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.1) +
    labs(x = "", y = "") +
    scale_color_manual(values = cols, name = "") +
    theme_cowplot(16) +
    ggtitle("") +
    geom_point(size = 3, shape = 21, fill = "black") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
}

gbm <- unique(data$annotation)

gaba <- create_plot(data, gbm[1]) + facet_wrap(~annotation)+
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))
estradiol <- create_plot(data,gbm[2])+ facet_wrap(~annotation)+
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))+
  labs(y="Effect size")
acetate <- create_plot(data,gbm[3])+ facet_wrap(~annotation)+
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))
trp <- create_plot(data,gbm[4])+ facet_wrap(~annotation)+
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))
butsyn <- create_plot(data,gbm[5])+ facet_wrap(~annotation) +
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))+
  labs(y="Effect size")
  

# Process the results files --

# Define the file paths, cohort prefixes, and other parameters
shotgun_fp <- c("Shotgun/CS_SPF/GBM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/HUM_Gavage/GBM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/SPF_Gavage/GBM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/UCLA_O_SPF/GBM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv")
shotgun_prefix <- c("CS SPF",
                    "UCLA O. SPF",
                    "SPF Gavage",
                    "HUM SD Gavage")

feature_value <- gbm_of_interest[1]
feature_annotation <- names(gbm_of_interest[1])
data <- process_gbm_files_shotgun(shotgun_fp, feature_value, shotgun_prefix,feature_annotation)

GBM_shotgun_df <- data[FALSE,]

for (i in seq_along(gbm_of_interest)){
  feature_value <- gbm_of_interest[i]
  feature_annotation <- names(gbm_of_interest[i])
  data <- process_gbm_files_shotgun(shotgun_fp, feature_value, shotgun_prefix,feature_annotation)
  GBM_shotgun_df <- rbind(data,GBM_shotgun_df)
}

v <-names(gbm_of_interest)
res_plot <- GBM_shotgun_df%>% select(c("coef", "qval","Cohort","Annotation"))
res_plot$Annotation <- factor(res_plot$Annotation, levels=c(v[1], v[2],v[3],v[4],v[5]))
res_plot$Annotation_Cohort <- paste0(res_plot$Annotation,"_",res_plot$Cohort)
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "Distal_Colon", "Jejunum"))

y = tapply(res_plot$coef, res_plot$Annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$Annotation= factor(as.character(res_plot$Annotation), levels = names(y))

names(my_palette) <-levels(res_plot$Cohort)
cols=c("#440154FF", "#FDE725FF")

gbm_shotgun_plot <- res_plot %>%
  # filter(qval < 0.05, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, Cohort, fill = site)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+  
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(12) +
  theme(plot.title = element_text(face = "plain")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (Jejunum/Distal_Colon)",
       y = "",
       fill = "") +
  ggtitle("Shotgun Data") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top",
        legend.justification = "center")+
  facet_wrap(~Annotation)+
  theme(legend.position="none",strip.background=element_rect(colour="black",fill="white"))
  

plot_grid(butsyn,trp, acetate, labels=c("A","B","C"),nrow=1)
plot_grid(estradiol,gaba,labels=c("D","E"))
dev.new()
plot_grid(gbm_upset, label_size = 20, labels="F")
plot_grid(gbm_shotgun_plot, label_size = 20, labels="G")

plot_grid(gbm_upset, gbm_shotgun_plot, labels=c("F","G"),label_size = 20,rel_widths = c(1.25,1))

### Upset Plots ---

file_paths <- c(
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GBM-Maaslin2-SITE/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "CS-Facility-Analysis/OMIXER-RPM Results/CS_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "UCLA_V_SPF_Analysis/OMIXER-RPM/WT_Val_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/Hum_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/SPF_GBM/GBM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
  "Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
  "Donors-Analysis/differential_GBM_site/GBM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"
)


cohort_prefixes <- c("UCLA_O_SPF", 
                     "UCLA_O_SPF",
                     "CS_SPF",
                     "CS_SPF",
                     "UCLA_V_SPF",
                     "HUM_Gavage",
                     "HUM_Gavage",
                     "SPF_Gavage",
                     "SPF_Gavage",
                     "HUM_V_Gavage",
                     "HUM_V_Gavage")

all_significant_gbm <- process_results_for_upset_plot(file_paths = file_paths,
                                          cohort_prefixes = cohort_prefixes)

all_significant_gbm_mod <- all_significant_gbm %>% select(c("feature", "Cohort")) %>% unique()

module_key <- read.csv(here("GBM_Module_Key.csv"))
anno <- module_key %>% select(c("feature", "annotation"))
all_taxa <- merge(all_significant_gbm_mod, anno, by="feature")
all_taxa <- all_taxa %>% select(c("feature", "Cohort","annotation")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

df_wide <- as.data.frame(df_wide)
df_wide <- df_wide %>% mutate(SPF_Gavage = 0)
all_datasets <- names(df_wide)[-c(1:2)]

ComplexUpset::upset(df_wide, all_datasets,width_ratio=0.1,
                                  base_annotations=list(
                                    'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
                                      scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')))+
  theme_cowplot(12)


dev.new()
ComplexUpset::upset(df_wide,
                    all_datasets, width_ratio = 0.1,
                    annotations = list(
                      'Gut-Brain Modules'=(
                        ggplot(mapping=aes(fill=annotation))
                        + geom_bar(stat='count', position='fill')
                        + scale_y_continuous(labels=scales::percent_format())
                        + scale_fill_manual(values=fill_color)
                        + ylab('GBMs')
                        + theme(legend.position="none")
                      )
                    )) 

forthelegend<-ComplexUpset::upset(df_wide,
                                  all_datasets, width_ratio = 0.1,
                                  annotations = list(
                                    'Metabolic'=(
                                      ggplot(mapping=aes(fill=annotation))
                                      + geom_bar(stat='count', position='fill')
                                      + scale_y_continuous(labels=scales::percent_format())
                                      + scale_fill_manual(values=fill_color)
                                      + ylab('Features')
                                      + theme(legend.position="right")
                                    )
                                  )) 

legend <- cowplot::get_legend(forthelegend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)

### Line Plots ---
# GABA degradation
feature_value <- "MGB019"
data_all_2 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_2 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_2 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# 17 B- estradiol degradation
feature_value <- "MGB031"
data_all_3 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_3 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_3 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# NO degradation I
feature_value <- "MGB027"
data_all_4 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_4 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_4 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Tryptophan degradation
feature_value <- "MGB049"
data_all_5 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_5 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_5 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Butyrate Synthesis I 
feature_value <- "MGB052"
data_all_6 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
data_all_6 <- process_results_files(file_paths_2, feature_value, new_value, new_coef, cohort_prefixes_2)
data_all_6 <- process_results_files(file_paths_3, feature_value, new_value, new_coef, cohort_prefixes_3)

# Plot the combined data frame
paletteer::palettes_d_names
my_palette <- c(paletteer_d("basetheme::brutal",10), paletteer_d("basetheme::clean",1))
data_all_2$Cohort <- factor(data_all_2$Cohort)
names(my_palette) <-levels(data_all_2$Cohort)

plot_data <- function(data, title,color_palette) {
  # Update 'value' column with labels
  data$value <- plyr::revalue(data$value,c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value <- factor(data$value, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  # Plotting
  plot <- ggplot(data, aes(x = value, y = coef, group = Cohort, color = Cohort)) +
    geom_line(size = 2) +
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.1) +
    labs(x = "Value", y = "Coefficient") +
    scale_color_manual(values = color_palette, name="") +
    theme_cowplot(16) +
    ggtitle(title) +
    geom_point(size = 3, shape = 21, fill = "black") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}

gaba <- plot_data(data_all_2, "GABA degradation") 

gaba_legend <- plot_data(data_all_2, "GABA degradation") + 
  theme(legend.position = "right")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(1, 1, 10, 10)) 
legend <- cowplot::get_legend(gaba_legend)
grid::grid.newpage()
grid::grid.draw(legend)

estradiol <- plot_data(data_all_3, "17-B-Estradiol degradation" )

no_degradation <- plot_data(data_all_4, "Nitric oxide degradation I (NO dioxygenase)" )

trp_degradation <- plot_data(data_all_5, "Tryptophan degradation" )

butyrate_synthesis <- plot_data(data_all_6, "Butyrate synthesis" )

plot_grid(gaba, estradiol, no_degradation,
          trp_degradation, butyrate_synthesis, 
          labels=c("A","B","C","D","E"),ncol = 2,label_size = 20)

