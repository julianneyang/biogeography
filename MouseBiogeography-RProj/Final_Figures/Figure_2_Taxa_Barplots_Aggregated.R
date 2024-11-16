###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 10.20.2022
### Figure Number: 1, allegedly
### Figure Contents: Mucosal alpha, beta, and DAT aggregated across datasets
###### whining ends here ---

library(cowplot)
library(ggplot2)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)
library(paletteer)
library(grid)
library(readr)
library(rlang)
library(here)
library(seecolor)
#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

### Taxa Barplots ---
here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_2_Taxa_Barplots_Aggregated.R")

compare_vector <- list(c("DC", "PC"),
                       c("DC", "Cec"),
                       c("DC", "Ile"),
                       c("DC", "Jej"),
                       c("DC", "Duo"))

## L6 level: Generating a global genera key (assign_cols.RDS used to be specific to each dataset)
hum_v_genera <- readr::read_rds(here("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_assign_cols.RDS"))
hum_v_genera <- names(hum_v_genera)
print(hum_v_genera)
cs_genera <- readr::read_rds(here("CS_SPF/Taxa-Barplots/assign_cols.RDS"))
  cs_genera <- names(cs_genera)
  print(cs_genera)
ucla_v_genera <- readr::read_rds(here("UCLA_V_SPF_Analysis/Taxa-Barplots/assign_cols.RDS"))
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
global_genera <- union(global_genera, hum_v_genera)
length(global_genera)
df<-palettes_d_names # see palette names
add_cols2 <- paletteer_d("ggthemes::Classic_20",20)	
add_cols4 <- paletteer_d("ggthemes::calc",12)
add_cols3 <- paletteer_d("dutchmasters::little_street",11)
add_cols5 <- paletteer_d("RColorBrewer::Spectral",2)
global_genera_cols <- c(add_cols2,add_cols3,add_cols4,add_cols5)
global_genera_cols <- unique(global_genera_cols)
names(global_genera_cols) <- global_genera
seecolor::print_color(global_genera_cols)
#readr::write_rds(global_genera_cols, here("global_genera_cols.RDS"))

## L2 level: Generate a global phyla key 
hum_v_lum_phyla <- get_phyla_from_plot("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-2.csv")
hum_v_muc_phyla <- get_phyla_from_plot("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-2.csv")

cs_lum_phyla <- get_phyla_from_plot("CS_SPF/Taxa-Barplots/Luminal_level-2.csv")
cs_muc_phyla <- get_phyla_from_plot("CS_SPF/Taxa-Barplots/Mucosal_level-2.csv")
ucla_v_phyla <- get_phyla_from_plot("UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-2.csv")
ucla_o_lum_phyla <- get_phyla_from_plot("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Luminal_level-2.csv")
ucla_o_muc_phyla <- get_phyla_from_plot("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Mucosal_level-2.csv")
hum_gavage_lum_phyla <- get_phyla_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Luminal_level-2.csv")
hum_gavage_muc_phyla <- get_phyla_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Mucosal_level-2.csv")
spf_gavage_lum_phyla <- get_phyla_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-2.csv")
spf_gavage_muc_phyla <- get_phyla_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-2.csv")
global_phyla <- c(hum_v_lum_phyla, hum_v_muc_phyla,
                  cs_lum_phyla, cs_muc_phyla,
                  ucla_v_phyla,
                  ucla_o_lum_phyla,ucla_o_muc_phyla,
                  hum_gavage_lum_phyla,hum_gavage_muc_phyla,
                  spf_gavage_lum_phyla,spf_gavage_muc_phyla)
global_phyla <- unique(global_phyla)
colors<- viridis::inferno(8)
seecolor::print_color(colors)
names(colors) <- global_phyla

names(colors) <- c("Actinobacteriota","Bacteroidota","Cyanobacteria","Deferribacterota","Desulfobacterota","Firmicutes","Proteobacteria","Verrucomicrobiota")
#readr::write_rds(colors, here("global_phyla_cols.RDS"))

# Donors 
genera_cols <- readRDS("global_genera_cols.RDS")
print(genera_cols)
genera_cols <- genera_cols[names(genera_cols) %in% hum_v_genera]
print(genera_cols)

hum_v_L6_muc <- generate_L6_taxa_plots("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-6.RDS",
                                        "HUM MD Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
hum_v_L6_muc

hum_v_L6_lum <- generate_L6_taxa_plots("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-6.RDS",
                                       "HUM MD Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
hum_v_L6_lum

phyla_cols <- readr::read_rds(here("global_phyla_cols.RDS"))
hum_v_phyla <- unique(c(hum_v_lum_phyla,hum_v_muc_phyla))
phyla_cols <- phyla_cols[names(phyla_cols) %in% ucla_v_phyla]
print(phyla_cols)

hum_v_L2_muc <-generate_L2_taxa_plots("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-2.csv", "HUM MD Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
hum_v_L2_lum <-generate_L2_taxa_plots("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-2.csv", "HUM MD Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))


# UCLA Original SPF
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% ucla_o_genera]
print(genera_cols)

UCLA_o_L6_muc <- generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Mucosal_level-6.RDS",
                                                                 "UCLA O. SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
UCLA_o_L6_muc

UCLA_o_L6_lum <-generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Luminal_level-6.RDS",
                                                                 "UCLA O. SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  #labs(subtitle = "Luminal Genera")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
UCLA_o_L6_lum

phyla_cols <- readRDS(here("global_phyla_cols.RDS"))
UCLA_o_L2_lum <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Luminal_level-2.csv", "UCLA O. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
UCLA_o_L2_muc<- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots/Mucosal_level-2.csv", "UCLA O. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  ylab("Relative abundance (%)")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))

# UCLA Validation SPF
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% ucla_v_genera]
print(genera_cols)

UCLA_v_L6_muc <- generate_L6_taxa_plots("UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-6.RDS",
                                        "UCLA V.SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
UCLA_v_L6_muc

phyla_cols <- readr::read_rds(here("global_phyla_cols.RDS"))
phyla_cols <- phyla_cols[names(phyla_cols) %in% ucla_v_phyla]
print(phyla_cols)

UCLA_v_L2_muc <-generate_L2_taxa_plots("UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-2.csv", "UCLA V. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))

UCLA_v_L2_muc

# SPF gavage
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% spf_gavage]
print(genera_cols)

spf_gavage_L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.RDS",
                                        "SPF Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
spf_gavage_L6_muc

spf_gavage_L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-6.RDS",
                                            "SPF Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
spf_gavage_L6_lum


phyla_cols <- readRDS("global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% spf_gavage_lum_phyla]

spf_gavage_L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-2.csv", "SPF Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
spf_gavage_L2_muc <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-2.csv", "SPF Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))


# HUM gavage
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% hum_gavage]
print(genera_cols)

hum_gavage_L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Luminal_level-6.RDS",
                                            "HUM SD Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
hum_gavage_L6_lum

hum_gavage_L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Mucosal_level-6.RDS",
                                            "HUM SD Gavage", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
hum_gavage_L6_lum

phyla_cols <- readRDS("global_phyla_cols.RDS")
hum_gavage_phyla <- unique(c(hum_gavage_lum_phyla,hum_gavage_muc_phyla))
phyla_cols <- phyla_cols[names(phyla_cols) %in% hum_gavage_phyla]
print(phyla_cols)

hum_gavage_L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Luminal_level-2.csv", "HUM SD Gavage", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))

hum_gavage_L2_muc <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Mucosal_level-2.csv","HUM SD Gavage", ".*p__", phyla_cols, "Site")+
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))

# CS SPF 
genera_cols <- readRDS("global_genera_cols.RDS")
genera_cols <- genera_cols[names(genera_cols) %in% cs_genera]
print(genera_cols)

cs_L6_muc <- generate_L6_taxa_plots("CS_SPF/Taxa-Barplots/Mucosal_level-6.RDS",
                                            "CS SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none") +
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
cs_L6_muc

cs_L6_lum <- generate_L6_taxa_plots("CS_SPF/Taxa-Barplots/Luminal_level-6.RDS",
                                    "CS SPF", ".*g__",genera_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
cs_L6_lum

cs_phyla <- unique(c(cs_lum_phyla,cs_muc_phyla))
phyla_cols <- readRDS("global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% cs_phyla]
print(phyla_cols)

cs_L2_lum <- generate_L2_taxa_plots("CS_SPF/Taxa-Barplots/Luminal_level-2.csv", "CS SPF", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))
cs_L2_muc <-generate_L2_taxa_plots("CS_SPF/Taxa-Barplots/Mucosal_level-2.csv", "CS SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x=element_text(color="black"),
        axis.text.y=element_text(color="black"),
        axis.ticks=element_line(color="black"),
        plot.title=element_text(face="plain"))

# Aggregate L6 Muc and draw legend
aggregated_L6_muc <- plot_grid(UCLA_o_L6_muc, UCLA_v_L6_muc, cs_L6_muc, spf_gavage_L6_muc, hum_gavage_L6_muc, hum_v_L6_muc,nrow=1, ncol=6,
rel_widths = c(1,1,1,1,1,1),
labels=c("","","","",""),label_size = 20)

dev.new(width=20, height=5)
aggregated_L6_muc

aggregated_L2_muc <- plot_grid(UCLA_o_L2_muc, UCLA_v_L2_muc,cs_L2_muc, spf_gavage_L2_muc, hum_gavage_L2_muc, hum_v_L2_muc,nrow=1, ncol=6)
dev.new(width=20, height=5)
aggregated_L2_muc

dev.new(width=20, height=5)
fig_mucosal <- plot_grid(alpha_diversity_muc, interregional_muc, interregional_mc, interregional_msi, aggregated_L6_muc, aggregated_L2_muc, ncol=1, nrow=6)

dev.new(width=20, height=5)
fig_mucosal

aggregated_L6_lum <- plot_grid(UCLA_o_L6_lum, cs_L6_lum, spf_gavage_L6_lum, hum_gavage_L6_lum, hum_v_L6_lum,NULL, NULL, nrow=1, ncol=6,
                               rel_widths = c(1,1,1,1,1,1),
                               labels=c("","","","",""),label_size = 20)
dev.new(width=20, height=20)
aggregated_L6_lum

aggregated_L2_lum <- plot_grid(UCLA_o_L2_lum, cs_L2_lum, spf_gavage_L2_lum, hum_gavage_L2_lum, hum_v_L2_lum,NULL, NULL, nrow=1, ncol=6,
                               rel_widths = c(1,1,1,1,1,1),
                               labels=c("","","","",""),label_size = 20)
aggregated_L2_lum
dev.new(width=20, height=5)
plot_grid(aggregated_L6_lum,aggregated_L2_lum,nrow=2)


genera_cols <- readRDS("global_genera_cols.RDS")
dummyplot<- as.data.frame(genera_cols)
dummyplot$dummyy <- seq(1,45,1)
dummyplot$dummyx <- seq(1,90,2)
dummyplot$Genus <- row.names(dummyplot)
L6_legend <-  ggplot(dummyplot, aes(x=dummyx,y=Genus,fill=Genus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=genera_cols,name="Genus")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol=6, byrow=TRUE))+
  #theme(legend.text = element_text(
    #margin = margin(r = 10, unit = "pt")))
  theme_cowplot(14)+
  theme(legend.spacing.x = unit(1, 'cm')) 
  #theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(10, 10, 10, 10)) 
legend <- cowplot::get_legend(L6_legend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)

phyla_cols <- readRDS("global_phyla_cols.RDS")
dummyplot<- as.data.frame(phyla_cols)
dummyplot$dummyy <- seq(1,8,1)
dummyplot$dummyx <- seq(1,16,2)
dummyplot$Genus <- row.names(dummyplot)
L2_legend <- ggplot(dummyplot, aes(x = dummyx, y = Genus, fill = Genus)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = phyla_cols, name = "Phylum") +
  theme_cowplot(14) +
  theme(
    legend.position = "right",
    legend.spacing.y = unit(3, "cm"),  # Increase spacing
    legend.key.height = unit(1.5, "cm")  # Increase key height
  ) +
  guides(fill = guide_legend(nrow = 8, byrow = TRUE))
  #theme(legend.background = element_rect(fill="white", size=3, linetype="solid"), legend.margin = margin(10, 10, 100, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
dev.new(width=20, height=40)
grid.draw(legend)

