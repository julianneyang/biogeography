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


### Taxa Barplots ---
assign_cols <- readRDS("CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS")
L2_muc <-generate_L2_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "Mucosal Phyla", ".*p__", "Site") +
  theme(legend.position = "none")

fig7e <- plot_grid(L2_lum,L2_muc, NULL, ncol=3,labels=c("E",""), rel_widths = c(2,2,1))

L6_muc <-generate_L6_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__",assign_cols, "Site") +
  theme(legend.position = "none")

fig7f <- plot_grid(L6_lum,L6_muc, NULL, ncol=3,labels=c("F",""), rel_widths = c(1.5,1.5,1))

#handle the genera names 
wrangle_genera_names <- function(csv_dataframe, filepathstring, rds_string){
input_data<-read.csv(csv_dataframe, row.names=1, header=TRUE)
taxa<-colnames(input_data)
colnames <- strsplit(taxa, ".f__")

family=new_list(length(colnames(input_data)))
i=1
for (i in 1:length(colnames)) {
  family[i] <- colnames[[i]][2]
  i=i+1
}

family<-unlist(family)
family <- strsplit(family, ".g__")

genus =new_list(length(colnames(input_data)))
i=1
for (i in 1:length(family)) {
  genus[i] <- family[[i]][2]
  i=i+1
}

family<-as.list(family)

i=1
for (i in 1:length(genus)) {
  if (isFALSE(genus[[i]]=="NA")) {
    genus[[i]] = genus[[i]] 
  }
  else {
    
    genus[[i]] <- paste0(family[[i]]," (f)")   
  }
  i=i+1
}
colnames(input_data) <-as.character(genus)

write_rds(input_data, paste0(filepathstring, rds_string))

}

wrangle_genera_names("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.csv", "ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/","Mucosal_level-6.RDS")

## Extract taxa from Mucosal ---
L2_lum<-readRDS("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.RDS")
L2_lum<- as.matrix(L2_lum)
L2_lum<-make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa/6
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
L2_lum <-select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_muc <- unique(L2_lum$Taxa)

## Generate a color key using paletteer colors ---
length(labels_muc)
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 14)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_muc
write_rds(assign_cols,"ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/assign_cols.RDS")

