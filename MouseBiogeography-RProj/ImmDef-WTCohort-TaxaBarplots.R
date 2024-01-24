library(ggplot2)
library(tidyr)
library(viridis)
library(cowplot)
library(plyr)
library(dplyr)
library(rlang)
library(funrar)
library(sjmisc)
library(RColorBrewer)
library(paletteer)
library(readr)
library(here)

here::i_am("MouseBiogeography-RProj/ImmDef-WTCohort-TaxaBarplots.R")


### Wrangle Genera Names --- 

# change the csv file into RDS 
file_path <- "UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-6.csv"
test <- readr::read_csv(here(file_path))
processed_data <- process_taxonomy_data(file_path)
print(names(processed_data))
readr::write_rds(processed_data, here("UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-6.RDS"))

## Generate a color key using paletteer colors --
labels_muc <- get_genera_from_plot("UCLA_V_SPF_Analysis/Taxa-Barplots/Mucosal_level-6.RDS")

#Find out how many taxa need to be assigned colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 14)

#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_muc
readr::write_rds(assign_cols,here("UCLA_V_SPF_Analysis/Taxa-Barplots/assign_cols.RDS"))
