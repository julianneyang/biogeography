library(paletteer)
library(dplyr)
library(here)
library(rlang)
library(funrar)


setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library(Microbiome.Biogeography)
setwd("/home/julianne/Documents/biogeography/")


here::i_am("MouseBiogeography-RProj/Shotgun_TaxaBarplots.R")

### Wrangle Genera Names --- 

# change the csv file into RDS 
file_path <- "Shotgun/starting_files/export_groupby_Site_BioGeo_Shotgun_L2/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

colnames(input_data) <- gsub("p__","",colnames(input_data))
write.csv(input_data,here("Shotgun/taxa_barplots/Shotgun_L2.csv"))

file_path <- "Shotgun/starting_files/export_groupby_Site_BioGeo_Shotgun_L6/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))


readr::write_rds(processed_data, here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.RDS"))

## Generate a color key using paletteer colors --
labels_lum <- get_genera_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Luminal_level-6.RDS")
labels_muc <- get_genera_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/SPF_Gavage_Mucosal_level-6.RDS")

#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 20)
add_cols <- paletteer_d("awtools::mpalette",1)
assign_cols <- c(assign_cols, add_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
readr::write_rds(assign_cols,here("Humanized-Biogeography-Analysis/taxa_barplots/SPF_assign_cols.RDS"))
