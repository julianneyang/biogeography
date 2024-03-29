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


?generate_L6_taxa_plots()
here::i_am("MouseBiogeography-RProj/CS_Facility_TaxaBarplots.R")

generate_L2_taxa_plots <- function(input_data, titlestring,greppattern, graphby){
  titlestring<-c(titlestring)
  L2_lum<-read.csv(input_data)
  L2_lum <- as.data.frame(t(L2_lum))
  colnames(L2_lum)<- L2_lum[1,]
  L2_lum <- L2_lum[-1,]
  taxa<-row.names(L2_lum)
  L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
  L2_lum <- as.data.frame(prop.table(L2_lum,2))
  taxa<-gsub(greppattern,"",taxa )
  
  L2_lum$Taxa <-taxa
  L2_lum<- pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  if({{graphby}} == "Site"){
    L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  }
  else if ({{graphby}}=="Type"){ 
    L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
    }
  L2_lum$Value <- L2_lum$Value * 100
  
  
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_viridis_d(option="B")+
    theme(legend.position = "right")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="")+
    ggtitle({{titlestring}}) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
  
}

L2_lum <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-2.csv", "CS Facility Luminal", ".*p__", "Site") 
L2_muc <-generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "CS Facility Mucosal", ".*p__", "Site")

dev.new(width=12, height=10)
plot_grid(L2_lum, L2_muc, align="hv")

L2_col <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Colon_level-2.csv", "CS Facility Colon", ".*p__", "Type")
L2_SI <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/SI_level-2.csv", "CS Facility SI", ".*p__", "Type")

dev.new(width=12, height=10)
plot_grid(L2_SI, L2_col, align="hv")

### Wrangle Genera Names --- 
process_taxonomy_data <- function(file_path) {
  input_data <- readr::read_csv(here(file_path))
  input_data <- as.data.frame(input_data)
  row.names(input_data)<- input_data[,1]
  input_data <- input_data[,-1]
                              
  
  taxa<-colnames(input_data)
  taxa <- gsub(";",".",taxa)
  colnames <- strsplit(taxa, ".o__")
  
  order=new_list(length(colnames(input_data)))
  i=1
  for (i in 1:length(colnames)) {
    order[i] <- colnames[[i]][2]
    i=i+1
  }
  
  order<-unlist(order)
  order <- strsplit(order, ".f__")
  
  family =new_list(length(colnames(input_data)))
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
  
  input_data <- input_data[,!grepl("Mitochondria", colnames(input_data))] 
  input_data <- input_data[,!grepl("Chloroplast", colnames(input_data))] 
  return(input_data)
}

# change the csv file into RDS 
file_path <- "CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS"))

file_path <- "CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.csv"
processed_data <- process_taxonomy_data(file_path)
readr::write_rds(processed_data, here("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS"))

## Generate a color key using paletteer colors --

get_genera_from_plot <- function(filepath){
L2_lum<-readr::read_rds(here(filepath))
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
return(unique(L2_lum$Taxa))
}
labels_lum <- get_genera_from_plot("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS")
labels_muc <- get_genera_from_plot("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS")

#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 20)
add_cols <- paletteer_d("dutchmasters::milkmaid",3)
assign_cols <- c(assign_cols,add_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
readr::write_rds(assign_cols,here("CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS"))

dev.new(width=12, height=10)
plot_grid(L6_lum,L6_muc, align="hv")
plot_grid(L2_lum,L6_lum,L2_muc,L6_muc, nrow=1, ncols=4,align="hv",labels = c("F","G", "H", "I"))
L6_col <-generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Colon_level-6.RDS","Colon ( > 0.1% Relative Abundance)", ".*g__",assign_cols)
L6_SI <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/SI_level-6.RDS", "SI ( > 0.1% Relative Abundance)", ".*g__",assign_cols)
dev.new(width=12, height=10)
plot_grid(L6_col, L6_SI, align="hv")
plot_grid(L2_SI, L6_SI, L2_col, L6_col, nrow=1, ncols=2,align="hv",labels = c("F","G", "H", "I"))

assign_cols <- readRDS("CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS")
print(assign_cols)


generate_L6_taxa_plots <- function(path_to_RDS, titlestring,greppattern, fillvector, graphby){
  titlestring<-c(titlestring)
  L2_lum<-readRDS(path_to_RDS)
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
    theme(legend.position = "none")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

lum_rds_filepath <- print(here("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS"))
L6_lum <- generate_L6_taxa_plots(lum_rds_filepath, "Luminal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") 
muc_rds_filepath <- print(here("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS"))
L6_muc <-generate_L6_taxa_plots(muc_rds_filepath,"Mucosal ( > 0.1% Relative Abundance)", ".*g__",assign_cols, "Site")

