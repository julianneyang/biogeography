library(paletteer)
library(dplyr)
library(here)
library(rlang)
library(funrar)


setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library(Microbiome.Biogeography)

here::i_am("MouseBiogeography-RProj/Hum-TaxaBarplots-SPF.R")

### Wrangle Genera Names --- 

# change the csv file into RDS 
file_path <- "Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-6.RDS"))

file_path <- "Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.RDS"))

## Generate a color key using paletteer colors --
labels_lum <- get_genera_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-6.RDS")
labels_muc <- get_genera_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.RDS")

#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 20)
add_cols <- paletteer_d("awtools::mpalette",2)
assign_cols <- c(assign_cols, add_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
readr::write_rds(assign_cols,here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_assign_cols.RDS"))
