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
library(paletteer)
library(grid)

### Taxa Barplots ---
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_Taxa_Barplots_Aggregated.R")

compare_vector <- list(c("DC", "PC"),
                       c("DC", "Cec"),
                       c("DC", "Ile"),
                       c("DC", "Jej"),
                       c("DC", "Duo"))

## L6 level: Generating a global genera key (assign_cols.RDS used to be specific to each dataset)
cs_genera <- readr::read_rds(here("CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS"))
  cs_genera <- names(cs_genera)
  print(cs_genera)
ucla_v_genera <- readr::read_rds(here("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/assign_cols.RDS"))
  ucla_v_genera <- names(ucla_v_genera)
  print(ucla_v_genera)
ucla_o_genera <- readr::read_rds(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/assign_cols.RDS"))
  ucla_o_genera <- names(ucla_o_genera)
  print(ucla_o_genera)
  #ucla_o_genera <- gsub("eae\\..","eae (", ucla_o_genera)
  #ucla_o_genera <- gsub("f\\.","f)", ucla_o_genera)
spf_gavage <- readr::read_rds(here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_assign_cols.RDS"))
  spf_gavage <- names(spf_gavage)
  print(spf_gavage)
hum_gavage <- readr::read_rds(here("Humanized-Biogeography-Analysis/taxa_barplots/HUM_assign_cols.RDS"))
  hum_gavage <- names(hum_gavage)  
  print(hum_gavage)
  #hum_gavage <- gsub("eae\\..","eae (", hum_gavage)
  #hum_gavage <- gsub("f\\.","f)", hum_gavage)
  
global_genera <- union(cs_genera, ucla_v_genera)
global_genera <- union(global_genera, ucla_o_genera)
global_genera <- union(global_genera, spf_gavage)
global_genera <- union(global_genera, hum_gavage)
length(global_genera)
df<-palettes_d_names # see palette names
add_cols2 <- paletteer_d("ggthemes::Classic_20",20)	
add_cols4 <- paletteer_d("ggthemes::calc",12)
add_cols3 <- paletteer_d("dutchmasters::little_street",6)
global_genera_cols <- c(add_cols2,add_cols3,add_cols4)
global_genera_cols <- unique(global_genera_cols)
names(global_genera_cols) <- global_genera
readr::write_rds(global_genera_cols, here("global_genera_cols.RDS"))

# UCLA Original SPF
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% ucla_o_genera]
print(genera_cols)

UCLA_o_L6_muc <- generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L6.RDS",
                                                                 "UCLA O. SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
UCLA_o_L6_muc

UCLA_o_L6_lum <-generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Luminal_L6.RDS",
                                                                 "UCLA O. SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
UCLA_o_L6_lum

phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
UCLA_o_L2_lum <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Luminal_L2.csv", "UCLA O. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")
UCLA_o_L2_muc <-generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L2.csv", "UCLA O. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")

# Draw legend
L6_legend <- Microbiome.Biogeography::generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L6.RDS","Mucosal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid::grid.newpage()
grid::grid.draw(legend)
L2_legend <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L2.csv","Mucosal Phyla", ".*p__", phyla_cols,"Site")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid::grid.newpage()
grid::grid.draw(legend)

# UCLA Validation SPF
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% ucla_v_genera]
print(genera_cols)

UCLA_v_L6_muc <- generate_L6_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.RDS",
                                        "UCLA V.SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
UCLA_v_L6_muc

lummucphyla <- read.csv("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

UCLA_v_L2_muc <-generate_L2_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "UCLA V. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")

UCLA_v_L2_muc

# SPF gavage
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% spf_gavage]
print(genera_cols)

spf_gavage_L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal_level-6.RDS",
                                        "SPF Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
spf_gavage_L6_muc

spf_gavage_L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal_level-6.RDS",
                                            "SPF Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
spf_gavage_L6_lum

lummucphyla <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal/level-2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]

spf_gavage_L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal/level-2.csv", "SPF Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")
spf_gavage_L2_muc <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal/level-2.csv", "SPF Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")


# HUM gavage
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% hum_gavage]
print(genera_cols)

hum_gavage_L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS",
                                            "HUM Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
hum_gavage_L6_muc

hum_gavage_L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-6.RDS",
                                            "HUM Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
hum_gavage_L6_lum

lummucphyla <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

hum_gavage_L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-2.csv", "HUM Gavage", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")

hum_gavage_L2_muc <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-2.csv","HUM Gavage", ".*p__", phyla_cols, "Site")+
  theme(legend.position = "none")

# CS SPF 
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% cs_genera]
print(genera_cols)

cs_L6_muc <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS",
                                            "CS SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none") 
cs_L6_muc

cs_L6_lum <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS",
                                    "CS SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")
cs_L6_lum

lummucphyla <- read.csv("CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

cs_L2_lum <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-2.csv", "CS SPF", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")
cs_L2_muc <-generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "CS SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")

# Aggregate L6 Muc and draw legend
aggregated_L6_muc <- plot_grid(UCLA_o_L6_muc, cs_L6_muc, spf_gavage_L6_muc, hum_gavage_L6_muc, UCLA_v_L6_muc, nrow=1, ncol=5)
dev.new(width=20, height=5)
aggregated_L6_muc

aggregated_L2_muc <- plot_grid(UCLA_o_L2_muc, cs_L2_muc, spf_gavage_L2_muc, hum_gavage_L2_muc, UCLA_v_L2_muc, nrow=1, ncol=5)
dev.new(width=20, height=5)
aggregated_L2_muc

dev.new(width=20, height=5)
fig_mucosal <- plot_grid(alpha_diversity_muc, interregional_muc, interregional_mc, interregional_msi, aggregated_L6_muc, aggregated_L2_muc, ncol=1, nrow=6)

dev.new(width=20, height=5)
fig_mucosal

aggregated_L6_lum <- plot_grid(UCLA_o_L6_lum, cs_L6_lum, spf_gavage_L6_lum, hum_gavage_L6_lum,NULL, NULL, nrow=1, ncol=5)
dev.new(width=20, height=5)
aggregated_L6_lum

aggregated_L2_lum <- plot_grid(UCLA_o_L2_lum, cs_L2_lum, spf_gavage_L2_lum, hum_gavage_L2_lum,NULL, NULL, nrow=1, ncol=5)
dev.new(width=20, height=5)
aggregated_L2_lum

L6_L2_lum <-  cowplot::plot_grid(aggregated_L6_lum,aggregated_L2_lum, ncol=1, nrow=2,labels=c("A","B"), label_size = 20)
L6_L2_muc <-  cowplot::plot_grid(aggregated_L6_muc,aggregated_L2_muc, ncol=1, nrow=2,labels=c("C","D"), label_size = 20)

L2_lum_muc <- cowplot::plot_grid(aggregated_L2_lum,aggregated_L2_muc, ncol=1, nrow=2,labels=c("D",""), rel_widths = c(16,20)) 
dev.new(width=20, height=5)
L2_lum_muc

L6_lum_muc <- cowplot::plot_grid(aggregated_L6_lum,aggregated_L6_muc, ncol=1, nrow=2,labels=c("C",""), rel_widths = c(16,20)) 
dev.new(width=20, height=5)
L6_lum_muc

genera_cols <- readRDS("global_genera_cols.RDS")
L6_legend <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS",
                                    "CS SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=44, byrow=TRUE))+
  theme_cowplot(12)+
  theme(legend.spacing.y = unit(0.01, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 0)) 
legend <- cowplot::get_legend(L6_legend)
grid.newpage()
grid.draw(legend)

phyla_cols <- readRDS("global_phyla_cols.RDS")
L2_legend <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-2.csv", 
                                    "HUM Gavage", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "top") +
  guides(fill=guide_legend(nrow=2, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
grid.draw(legend)
