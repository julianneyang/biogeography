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

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_6_GMM_Coef_Plots.R")

### Upset Plot ---

lum_file_paths <- c("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                "CS_SPF/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Donors-Analysis/differential_GMM_site/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")


lum_cohort_prefixes <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_MD_Gavage",
                     "HUM_SD_Gavage",
                     "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = lum_file_paths,
                                           cohort_prefixes = lum_cohort_prefixes,
                                           filter_by = "Site")

module_key <- read.csv(here("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv"))
anno <- module_key %>% select(c("feature", "annotation"))
all_taxa <- merge(all_taxa, anno, by="feature")
id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","annotation","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)
id_df_wide <- id_df_wide %>% mutate(SPF_Gavage = 0)

all_taxa <- all_taxa %>% select(c("feature", "Cohort","annotation")) %>% unique()

df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)
df_wide <- as.data.frame(df_wide)
df_wide <- df_wide %>% mutate(SPF_Gavage = 0)

df_wide <- as.data.frame(df_wide)
all_datasets <- names(df_wide)[-(1:2)]
gmm_upset <- upset(df_wide, all_datasets,width_ratio=0.1,
                                  base_annotations=list(
                                    'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
                                      scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')),
                                 themes=list(
                                   default=theme(
                                     axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                   ),
                                   intersections_matrix=theme(
                                     axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                   )
                                 ))

id_df_wide$count_ones <- rowSums(id_df_wide[, c(4:8)])
df_filtered <- id_df_wide[id_df_wide$count_ones >= 3, ]
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
gmm_of_interest <- df_filtered$feature
names(gmm_of_interest) <-df_filtered$annotation
write_rds(gmm_of_interest, "Highlighted_GMM_Fig_6.RDS")

### Coef Plots ---

lum_file_paths <- c("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                    "CS_SPF/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                    "Donors-Analysis/differential_GMM_site/GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv",
                    "Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                    "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv")

lum_cohort_prefixes <- c("UCLA_O_SPF",
                         "CS_SPF",
                         "HUM_MD_Gavage",
                         "HUM_SD_Gavage",
                         "SPF_Gavage")

new_value <- "Distal_Colon"
new_coef <- 0

# color legend for coef plots 
my_palette <- c(paletteer_d("basetheme::brutal",6))
names(my_palette) <-c(lum_cohort_prefixes, "UCLA_V_SPF")
cols <- my_palette[names(my_palette) %in% lum_cohort_prefixes]


# combine GMM results and append Map annotation
feature_value <- gmm_of_interest[1]
data_all <- process_results_files(lum_file_paths, feature_value, new_value, new_coef, lum_cohort_prefixes)

final_df <- data_all[FALSE,]
for (i in seq_along(gmm_of_interest)) {
feature_value <- gmm_of_interest[i]
data_all <- process_results_files(lum_file_paths, feature_value, new_value, new_coef, lum_cohort_prefixes)
final_df <- rbind(final_df,data_all)
}

map <- module_key %>% select(c("feature", "annotation", "Map","Map2_ammonia", "Map3_carbon_dioxide"))
map_data_all <- merge(final_df,map,by="feature")
unique(map_data_all$Map)

data <- map_data_all
data$value <- plyr::revalue(data$value,c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
data$value <- factor(data$value, levels = c("D", "J", "I", "C", "PC", "DC"))
data$annotation <- gsub("pentose phosphate pathway \\(oxidative phase\\)", "pentose phosphate pathway", data$annotation)

# Plotting
create_plot <- function(data, map_value) {
  ggplot(data %>% filter(Map == map_value),
         aes(x = value, y = coef, group = Cohort, color = Cohort)) +
    geom_line(size = 2) +
    geom_hline(yintercept=0, linetype='dotted', col = 'black',linewidth = 2)+
    geom_errorbar(aes(ymin = coef - stderr, ymax = coef + stderr), width = 0.1) +
    labs(x = "", y = "") +
    scale_color_manual(values = cols, name = "") +
    theme_cowplot(16) +
    ggtitle("") +
    geom_point(size = 3, shape = 21, fill = "black") +
    theme(legend.position = "none") +
    theme(plot.title = element_text(hjust = 0.5))
}

# Call the function for both "disaccharides" and "monosaccharides"
disaccharides <- create_plot(data, "disaccharides") + facet_grid(Map~annotation)
monosaccharides <- create_plot(data, "monosaccharides")+ facet_grid(Map~annotation)
polysaccharides <- create_plot(data, "polysaccharides")+ facet_grid(Map~annotation)
proteolytic_fermentation <- create_plot(data, "proteolytic fermentation") + facet_grid(Map~annotation) +theme_cowplot(14) + theme(legend.position = "none")
create_plot(data, "proteolytic fermentation") + facet_grid(Map~annotation) 
lipolytic_fermentation <- create_plot(data, "lipolytic fermentation")+ facet_grid(Map~annotation)
sugar_acid <- create_plot(data, "sugar acid")+ facet_grid(Map~annotation)
metabolism <- create_plot(data, "central metabolism")+ facet_grid(Map~annotation)  +theme_cowplot(12) + theme(legend.position = "none")
cross_feeding <- create_plot(data, "cross-feeding")+ facet_grid(Map~annotation)
butyrate <- create_plot(data, "butyrate")+ facet_grid(Map~annotation)
nitrate <- create_plot(data, "nitrate reduction")+ facet_grid(Map~annotation)


carbs <- plot_grid(monosaccharides,disaccharides,polysaccharides,
          rel_widths=c(0.5,0.5,1),nrow=1,
          labels=c("A"), label_size = 20)
proteins <- plot_grid(proteolytic_fermentation,nrow=1,
          labels=c("B"), label_size = 20)

dev.new()
plot_grid(carbs, proteins, rel_heights=c(0.75,1),nrow=2)

fats <- plot_grid(lipolytic_fermentation,sugar_acid,nitrate,
          rel_widths=c(1,0.5,0.5),nrow=1,
          labels=c("C"), label_size = 20)


cross_feeding <- plot_grid(cross_feeding,butyrate,
          rel_widths=c(1,0.33),nrow=1,
          labels=c("D"), label_size = 20)

dev.new()
plot_grid(fats, cross_feeding, rel_heights=c(1,1),nrow=2)


fig6E <- plot_grid(metabolism, 
          labels=c("E"), label_size = 20)
fig6F <- plot_grid(gmm_upset, 
          labels=c("F"), label_size = 20)

dev.new()
plot_grid(fig6E, fig6F, rel_widths= c(0.75, 1),ncol=2)

plot_grid(GMM[[1]], GMM[[2]],GMM[[3]], GMM[[4]],
          GMM[[5]], GMM[[6]],ncol=3,nrow=2)

plot_grid(GMM[[7]], GMM[[8]],GMM[[9]], GMM[[10]], GMM[[11]], GMM[[12]],
           GMM[[13]], gmm_upset, nrow=3, ncol=3)
         
feature_value <- gmm_of_interest[2]
names(gmm_of_interest[2])
data_all_GMM2 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
GMM2 <- plot_data(data_all_GMM1, names(gmm_of_interest[2]))

feature_value <- gmm_of_interest[2]
data_all_GMM2 <- process_results_files(file_paths, feature_value, new_value, new_coef, cohort_prefixes)
GMM2 <- plot_data(data_all_GMM1, "Lactose and Galactose Degradation")

### Draw legend ---
L2_legend <-  sugar_acid + 
  theme(legend.position = "right") +
  #guides(fill=guide_legend(nrow=8, byrow=TRUE))+
  theme_cowplot(16)+
  theme(legend.spacing.y = unit(1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(1, 1, 1, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
grid.draw(legend)

### Shotgun GMM ---
# Define the file paths, cohort prefixes, and other parameters
shotgun_fp <- c("Shotgun/CS_SPF/GMM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/HUM_Gavage/GMM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/SPF_Gavage/GMM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                "Shotgun/UCLA_O_SPF/GMM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv")
shotgun_prefix <- c("CS SPF",
                    "HUM SD Gavage",
                    "SPF Gavage",
                    "UCLA O. SPF")


feature_value <- gmm_of_interest[1]
feature_annotation <- names(gmm_of_interest[1])
data <- process_gbm_files_shotgun(shotgun_fp, feature_value, shotgun_prefix,feature_annotation)

GMM_shotgun_df <- data[FALSE,]

for (i in seq_along(gmm_of_interest)){
  feature_value <- gmm_of_interest[i]
  feature_annotation <- names(gmm_of_interest[i])
  data <- process_gbm_files_shotgun(shotgun_fp, feature_value, shotgun_prefix,feature_annotation)
  GMM_shotgun_df <- rbind(data,GMM_shotgun_df)
}


res_plot <- GMM_shotgun_df%>% select(c("coef", "qval","Cohort","Annotation"))
res_plot$Annotation_Cohort <- paste0(res_plot$Annotation,"_",res_plot$Cohort)
res_plot <- unique(res_plot)
res_plot <- res_plot %>%
  mutate(site = ifelse(coef< 0, "Distal_Colon", "Jejunum"))

y = tapply(res_plot$coef, res_plot$Annotation, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
res_plot$Annotation= factor(as.character(res_plot$Annotation), levels = names(y))

names(my_palette) <-levels(res_plot$Cohort)
cols=c("#440154FF", "#FDE725FF")

res_plot$Cohort <- factor(res_plot$Cohort, levels=c("HUM SD Gavage","SPF Gavage","CS SPF", "UCLA O. SPF"))
res_plot$Annotation  <- gsub("pentose phosphate pathway \\(oxidative phase\\)", "pentose phosphate pathway", res_plot$Annotation)

res_plot %>%
  arrange(Annotation) %>%
  # filter(qval < 0.05, abs(coef) > 0) %>%
  ggplot2::ggplot(aes(coef, Cohort, fill = site)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black")+  
  geom_bar(stat = "identity") +
  cowplot::theme_cowplot(12) +
  theme(axis.text.y = element_text(face = "bold")) +
  scale_fill_manual(values = cols) +
  labs(x = "Effect size (Jejunum/Distal_Colon)",
       y = "",
       fill = "") +
  ggtitle("Shotgun Data") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "top")+
  facet_wrap(~Annotation)


### HUM V Gavage ---
donors_filepath <- "Donors-Analysis/differential_GMM_site/"
lumtarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

muctarget <- find_concordant_features_across_sites(paste0(donors_filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"))

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
hum_v_map_lum <- generate_GMM_heat_map_by_site(paste0(donors_filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                                lumtarget,
                                                "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                Y=Map,
                                                "Map",
                                                "Luminal",
                                                cols,
                                                bk) +
  theme_cowplot(20) +
  ggtitle("HUM V. Gavage Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

?generate_GMM_heat_map_by_site()
hum_v_mmap_lum <- generate_GMM_heat_map_by_site(paste0(donors_filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                               lumtarget,
                                               "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                               Y=metabolic_map,
                                               "metabolic_map",
                                               "Luminal",
                                               cols,
                                               bk) +
  theme_cowplot(20) +
  ggtitle("HUM V. Gavage Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

dev.new()
hum_v_mmap_lum
cols=c("#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-1, -0.5, 0, 0.5, 1, 1.5)


hum_v_map_muc <- generate_GMM_heat_map_by_site(paste0(donors_filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"),
                                               muctarget,
                                               "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                               Y=Map,
                                               "Map",
                                               "Mucosal",
                                               cols,
                                               bk) +
  theme_cowplot(20) +
  ggtitle("HUM V. Gavage Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))



### UCLA. O. SPF ---
lumtarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv")

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1)
ucla_o_map_lum <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk) +
  theme_cowplot(20) +
                  ggtitle("UCLA O. SPF Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

par(mar = rep(2, 4))
ucla_o_map_lum


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF")
bk =c(-2,-1.5, -1, -0.5, 0, 0.5, 1)
ucla_o_map_muc <- generate_GMM_heat_map_by_site("Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/GMM-Maaslin2-SITE/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                     muctarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Mucosal",
                                     cols,
                                     bk) +
  theme_cowplot(20)+
                  ggtitle("UCLA O. SPF Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


### CS SPF ---

# Luminal - 36 concordant features 
lum_target <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 37 concordant features 
muctarget <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
cs_muc_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Mucosal", 
                                             cols,
                                             bk) +
  
  theme_cowplot(20) +
                  ggtitle("CS SPF Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
cs_lum_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Luminal",
                                             cols,
                                             bk) +
  theme_cowplot(20) +
                  ggtitle("CS SPF Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
                  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))


### HUM Gavage ---

#Feed in the significant results and generate a target vector with the union of all features 
lumtarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF")
bk =c(-2, -1.5, -1, -0.5, 0)

hum_lum_GMM_map <- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                     lumtarget,
                                     "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                     Y=Map,
                                     "Map",
                                     "Luminal",
                                     cols,
                                     bk) +
  theme_cowplot(20) +
  ggtitle("HUM Gavage Lum") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
 



cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

hum_muc_GMM_map<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/Hum/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=Map,
                                    "Map",
                                    "Mucosal",
                                    cols,
                                    bk) +
  theme_cowplot(20) +
  ggtitle("HUM Gavage Muc") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

hum_muc_GMM_map

### SPF Gavage ---
muctarget <- find_concordant_features_across_sites("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

spf_muc_GMM_map<- generate_GMM_heat_map_by_site("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                    muctarget,
                                    "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                    Y=Map,
                                    "Map",
                                    "SPF Gavage Muc",
                                    cols,
                                    bk)+
  theme_cowplot(20) +
  ggtitle("SPF Gavage Muc") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

fig6h <- plot_grid(NULL,spf_muc_GMM_map, ncol=1)

### UCLA V. SPF ---
muctarget <- find_concordant_features_across_sites("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv ---
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
ucla_v_muc_GMM_map <- generate_GMM_heat_map_by_site("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             titlestring="UCLA V. SPF Muc",
                                             colorvector = cols,
                                             breakvector = bk) +
  theme_cowplot(20) +
  ggtitle("UCLA V. SPF Muc") +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "none") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

thelegend <- generate_GMM_heat_map_by_site("UCLA_V_SPF_Analysis/OMIXER-RPM/WTCohort_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                           targetvector = muctarget, 
                                           path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                           Y=Map,
                                           ystring= "Map",
                                           titlestring="UCLA V. SPF Muc",
                                           colorvector = cols,
                                           breakvector = bk) +
  theme_cowplot(16) +
  ggtitle("UCLA V. SPF Muc") +
  theme(legend.position = "top") +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))

legend <- cowplot::get_legend(thelegend)
grid::grid.newpage()
grid::grid.draw(legend)

### Assemble Figure 6 ---

dev.new(width=10,height=10)
plot_grid(ucla_o_map_lum, cs_lum_GMM_map, hum_lum_GMM_map, nrow=1, labels= c("A", "B", "F"), label_size = 20,greedy = FALSE)
dev.new(width=10,height=10)
plot_grid(ucla_o_map_muc, cs_muc_GMM_map, hum_v_map_lum, nrow=1, labels= c("C", "D", "G"), label_size = 20)
dev.new(width=10,height=10)
plot_grid(ucla_v_muc_GMM_map, hum_muc_GMM_map, hum_v_map_muc, nrow=1, labels = c("H", "I", "J"), label_size = 20)
dev.new(width=10,height=10)
plot_grid(spf_muc_GMM_map, NULL, NULL, nrow=1, labels = c("K", "", "L"), label_size = 20)
