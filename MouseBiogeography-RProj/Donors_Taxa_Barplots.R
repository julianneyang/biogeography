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
library(Microbiome.Biogeography)
library(devtools)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")


here::i_am("MouseBiogeography-RProj/Donors_Taxa_Barplots.R")
phyla_cols <- readRDS("../Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")

mouse_samples = readr::read_csv(here("Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-2.csv"))

A017_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-2.csv", 
                                 titlestring = "A017", 
                                 greppattern = ".*p__", 
                                 graphby = "Site",
                                 fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B") +
  theme(plot.margin = margin(r = -2))

A018_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A018_level-2.csv", 
                                                               titlestring = "A018", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")

A041_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A041_level-2.csv", 
                                                               titlestring = "A041", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")

A043_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A043_level-2.csv", 
                                                               titlestring = "A043", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")

A045_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A045_level-2.csv", 
                                                               titlestring = "A045", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")

A046_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A046_level-2.csv", 
                                                               titlestring = "A046", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")

A047_L2_lum <- Microbiome.Biogeography::generate_L2_taxa_plots(input_data = "../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A047_level-2.csv", 
                                                               titlestring = "A047", 
                                                               greppattern = ".*p__", 
                                                               graphby = "Site",
                                                               fillvector = phyla_cols) +
  scale_fill_viridis_d(option="B")



dev.new(width=12, height=10)
plot_grid(L2_lum, L2_muc, align="hv")


## Donor Feces -- 

feces <- readr::read_csv(here("Donors-Analysis/taxa_barplots/Hoomins_level-2.csv"))

generate_L2_human_donor_taxa_barplot <- function(dataframe, donor_id){
  
  L2_lum<-dataframe
  L2_lum <- L2_lum %>% filter(Donor_ID==donor_id)
  L2_lum <- as.data.frame(t(L2_lum))
  taxa<-row.names(L2_lum)
  colnames(L2_lum)<- L2_lum[1,]

  row.names(L2_lum) <- as.character(taxa)
  L2_lum$phyla <- row.names(L2_lum)

  L2_lum <- L2_lum %>%
    filter(startsWith(phyla, "k__Bacteria"))
  L2_lum <- L2_lum %>% select(-c("phyla"))
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  
  taxa <- taxa[grep("^k__", taxa)]
  taxa<-gsub(".*p__","",taxa )
  L2_lum$Taxa <-taxa
  L2_lum<- pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="SampleID")
  L2_lum$Value <- L2_lum$Value * 100


 plot<- ggplot(data=L2_lum, aes(x=SampleID, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_viridis_d(option="B")+
    theme(legend.position = "right")+
    theme_cowplot(12) +
    ylab("") +
    xlab("")+
    labs(fill="")+
    ggtitle(paste({{donor_id}},"Human Feces")) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
    return(plot)
}

A017 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A017")
A018 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A018")
A041 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A041")
A043 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A043")
A045 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A045")
A046 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A046")
A047 <- generate_human_donor_taxa_barplot(dataframe = feces,donor_id = "A047")

## Combine --
?align_plots()
plot_grid( A017_L2_lum, A017, rel_widths = c(1,0.33))

### Mouse Samples L6 plots ---

generate_L6_taxa_plots_donors <- function(filepath, titlestring,greppattern, fillvector, graphby){
  
input_data <- readr::read_csv(here(filepath))
input_data<-as.data.frame(input_data)
row.names(input_data) <- input_data[,1]
input_data <- input_data[,-1]

#input_data <- read.csv("../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-6.csv", header=TRUE,row.names=1)
taxa<-colnames(input_data)
colnames <- strsplit(taxa, ".o__")

order=rlang::new_list(length(colnames(input_data)))
i=1
for (i in 1:length(colnames)) {
  order[i] <- colnames[[i]][2]
  i=i+1
}

order<-unlist(order)
order <- strsplit(order, ".f__")

family =rlang::new_list(length(colnames(input_data)))
i=1
for (i in 1:length(order)) {
  family[i] <- order[[i]][2]
  i=i+1
}

order<-as.list(order)

i=1
for (i in 1:length(family)) {
  if (isFALSE(family[[i]]==".g__")) {
    family[[i]] = family[[i]] 
  }
  else {
    family[[i]] <- paste0(order[[i]]," (o)")   
    family[[i]] <- family[[i]][1]
  }
  i=i+1
}


family<-unlist(family)
family <- strsplit(family, ".g__")

genus =rlang::new_list(length(colnames(input_data)))
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


input_data <- input_data[,!grepl("Mitochondria", colnames(input_data))] 
input_data <- input_data[,!grepl("Chloroplast", colnames(input_data))] 

original_string <- filepath
new_filepath<- sub("\\.csv$", ".RDS", original_string)

readr::write_rds(input_data, here(new_filepath))
#write_rds(input_data, "Donors-Analysis/Mucosal_level-6.RDS")

  titlestring<-c(titlestring)
  L2_lum<-input_data
  
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
  taxa<-gsub(greppattern,"",taxa )
  #taxa<-gsub(".*g__","",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  }
  
  cols <- fillvector
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="right") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

genera_cols <- readr::read_rds(here("global_genera_cols.RDS"))

A017_L6_lum <- generate_L6_taxa_plots_donors(filepath = "Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-6.csv", 
                                        titlestring = "A017", 
                                        greppattern = ".*g__", 
                                        graphby = "Site",
                                        fillvector = assign_cols) +
  theme(legend.position = "right")

### Human Feces L6 Plots ---
generate_L6_human_feces_taxa_plots <- function(donor_id,filepath, titlestring,greppattern, fillvector, graphby){
  #input_data <- readr::read_csv(here("Donors-Analysis/taxa_barplots/Hoomins_level-6.csv"))
 
 input_data <- readr::read_csv(here(filepath))
 input_data <- input_data %>% filter(Donor_ID==donor_id)
 
  input_data<-as.data.frame(input_data)
  row.names(input_data) <- input_data[,1]
  input_data <- input_data[,-1]
  
  #input_data <- read.csv("../Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-6.csv", header=TRUE,row.names=1)
  taxa<-colnames(input_data)
  colnames <- strsplit(taxa, ".o__")
  
  order=rlang::new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(colnames)) {
    order[i] <- colnames[[i]][2]
    i=i+1
  }
  
  order<-unlist(order)
  order <- strsplit(order, ".f__")
  
  family =rlang::new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(order)) {
    family[i] <- order[[i]][2]
    i=i+1
  }
  
  order<-as.list(order)
  
  i=1
  for (i in 1:length(family)) {
    if (isFALSE(family[[i]]==".g__")) {
      family[[i]] = family[[i]] 
    }
    else {
      family[[i]] <- paste0(order[[i]]," (o)")   
      family[[i]] <- family[[i]][1]
    }
    i=i+1
  }
  
  
  family<-unlist(family)
  family <- strsplit(family, ".g__")
  
  genus =rlang::new_list(length(colnames(input_data)))
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
  
  input_data <- input_data[,!grepl("Mitochondria", colnames(input_data))] 
  input_data <- input_data[,!grepl("Chloroplast", colnames(input_data))] 
  
  original_string <- filepath
  new_filepath<- sub("\\.csv$", ".RDS", original_string)
  
  readr::write_rds(input_data, here(new_filepath))
  #write_rds(input_data, "Donors-Analysis/Mucosal_level-6.RDS")
  
  titlestring<-c(titlestring)
  L2_lum<-input_data
  names(L2_lum)
  L2_lum <- L2_lum %>%
    select(-starts_with("NA"))
  
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
  taxa<-gsub(greppattern,"",taxa )
  #taxa<-gsub(".*g__","",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = factor(L2_lum$Site)
  }
 
  cols <- fillvector
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(paste({{donor_id}},"Human Feces")) +
    theme(legend.position="right") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

generate_L6_human_feces_taxa_plots(filepath = "Donors-Analysis/taxa_barplots/Hoomins_level-6.csv",
                                   titlestring = "Human Feces",
                                   greppattern = ".*g__", 
                                   graphby="Site",
                                   fillvector = fecal_cols,
                                   donor_id = "A017")

### Aggregated taxa barplots ---
file_path <- "Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-6.RDS"))

file_path <- "Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-6.RDS"))

## Generate a color key using paletteer colors --
labels_lum <- get_genera_from_plot("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Luminal_level-6.RDS")
labels_muc <- get_genera_from_plot("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_Mucosal_level-6.RDS")

#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 14)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
readr::write_rds(assign_cols,here("Donors-Analysis/taxa_barplots/aggregated_barplots/Mice_assign_cols.RDS"))
