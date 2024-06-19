library(ggplot2) #yes
library(cowplot) #yes
#library(plyr)
library(dplyr) #yes
library(rlang) #yes
library(funrar) #yes
#library(sjmisc)
library(RColorBrewer)
library(paletteer)
library(readr)
library(here) #yes
library(Microbiome.Biogeography)
library(devtools)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")
here::i_am("MouseBiogeography-RProj/Figure_S_Donors_Taxa_Barplots.R")

## Make Donors taxa barplots -- 
phyla_cols <- readRDS("global_phyla_cols.RDS")

## Get names of donors that I need to make plots for
mouse_ASV <- read.delim(here("Donors-Analysis/starting_files/Donors-Mice-1xPrev0.15-ComBat-ASV.tsv"),row.names=1)
names(mouse_ASV) <- gsub("X","",names(mouse_ASV))
metadata <- read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"))
metadata$SampleID <- gsub("-",".",metadata$SampleID)
metadata <- metadata %>%
  filter(SampleID %in% names(mouse_ASV)) %>%
  filter(MouseID!= "U2") %>% 
  filter(Original_Human_Stool=="N")
metadata$Donor_ID <- factor(metadata$Donor_ID)
donors <- levels(metadata$Donor_ID)

## Loop over each donor and make plot of taxa by donor --
generate_plots <- function(donor) {
  input_file <- paste0("Donors-Analysis/taxa_barplots/plot_by_donor/", donor, "_Lum_level-2.csv")
  plot <- generate_L2_taxa_plots(input_data = input_file, 
                                 titlestring = paste0("(MD) ",donor), 
                                 greppattern = ".*p__", 
                                 graphby = "Site",
                                 fillvector = phyla_cols) +
    theme(plot.margin = margin(r = -2)) +
    theme_cowplot(12)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none")
  return(plot)
}

# Generate plots for each item
plots <- lapply(donors, generate_plots)
names(plots) <- paste(donors, "L2_lum", sep = "_")

# Access plots using list indexing
for (donor in donors) {
  assign(paste0(donor, "_L2_lum"), plots[[paste0(donor, "_L2_lum")]])
}

## Donor Feces -- 

feces <- readr::read_csv(here("Donors-Analysis/taxa_barplots/Hoomins_level-2.csv"))
feces$Donor_ID<-factor(feces$Donor_ID)
human_donors <- levels(feces$Donor_ID)

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
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="SampleID")
  L2_lum$Value <- L2_lum$Value * 100
  L2_lum$SampleID <- "F"

 plot<- ggplot(data=L2_lum, aes(x=SampleID, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=phyla_cols)+
    theme_cowplot(12) +
    ylab("") +
    theme(axis.text.y=element_blank()) +
    xlab("")+
    labs(fill="")+
    ggtitle(paste("Human")) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
    return(plot)
}

generate_L2_SD_human_donor_taxa_barplot <- function(dataframe){
  L2_lum <- A001
  L2_lum<-dataframe
  #L2_lum <- L2_lum %>% filter(Donor_ID==donor_id)
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
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="SampleID")
  L2_lum$Value <- L2_lum$Value * 100
  L2_lum$SampleID <- "F"
  
  plot<- ggplot(data=L2_lum, aes(x=SampleID, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_manual(values=phyla_cols)+
    theme_cowplot(12) +
    ylab("") +
    theme(axis.text.y=element_blank()) +
    xlab("")+
    labs(fill="")+
    ggtitle(paste("Human")) +
    theme(legend.position="none") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
  return(plot)
}



# Access plots using list indexing
generate_plots <- function(donor_ids, dataframe) {
  plots <- list()
  for (donor_id in donor_ids) {
    plots[[donor_id]] <- generate_L2_human_donor_taxa_barplot(dataframe = dataframe, donor_id = donor_id)
  }
  return(plots)
}

# Call the function to generate plots
plots <- generate_plots(donor_ids = human_donors, dataframe = feces)
names(plots) <- paste(human_donors, "feces", sep = "_")


# Access plots using list indexing
for (donor in human_donors) {
  assign(paste0(donor, "_feces"), plots[[paste0(donor, "_feces")]])
}

### Add HUM SD Donor ---
A001 <- read.csv("Humanized-Biogeography-Analysis/SD_Donor/taxa_barplots/HUM_SD_Donor_level-2.csv",row.names=1)
A001_feces <-  generate_L2_SD_human_donor_taxa_barplot(dataframe = A001)

A001_L2_lum<- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Luminal_level-2.csv", "(SD) A001", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")
  

### Make final figure ---
dev.new()
plot_grid(A001_L2_lum, A001_feces,
          A017_L2_lum, A017_feces,
          A018_L2_lum, A018_feces,
          A041_L2_lum, A041_feces,
          A043_L2_lum, A043_feces,
          A045_L2_lum, A045_feces,
          A046_L2_lum, A046_feces,
          A047_L2_lum, A047_feces,
          A050_L2_lum, A050_feces,
          A053_L2_lum, A053_feces,
          A054_L2_lum, A054_feces,
          A070_L2_lum, A070_feces,
          A078_L2_lum, A078_feces,
          A082_L2_lum, A082_feces,
          rel_widths = c(1,0.25,
                         1,0.25,
                         1, 0.25,
                         1, 0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1, 0.25,
                         1,0.25))
### Mouse Samples L6 plots ---

generate_L6_taxa_plots_donors <- function(filepath, titlestring,greppattern, fillvector, graphby){
#input_data <- read.csv(here("Donors-Analysis/taxa_barplots/plot_by_donor/A017_Lum_level-6.csv"))
input_data <- read.csv(here(filepath))
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
print(new_filepath)
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
    L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="C","Ileum"="I", "Jejunum"="J", "Duodenum"= "D"))
    L2_lum$Site = factor(L2_lum$Site, levels=c("D", "J", "I", "C", "PC", "DC"))
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

new_genera_legend <- readr::read_rds(here("Donors-Analysis/taxa_barplots/donors_genera_cols.RDS"))

## Loop over each donor and make plot of taxa by donor --
process_donor_data <- function(donor){
  file_path <- paste0("Donors-Analysis/taxa_barplots/plot_by_donor/", donor, "_Lum_level-6.csv")
  processed_data <- process_taxonomy_data(file_path)
  readr::write_rds(processed_data, here(paste0("Donors-Analysis/taxa_barplots/plot_by_donor/", donor, "_Lum_level-6.RDS")))
}

lapply(donors, process_donor_data)

generate_plots <- function(donor) {
  input_file <- paste0("Donors-Analysis/taxa_barplots/plot_by_donor/", donor, "_Lum_level-6.RDS")
  plot <- generate_L6_taxa_plots(path_to_RDS = input_file, 
                                 titlestring = paste0("(MD) ",donor), 
                                 greppattern = ".*g__", 
                                 graphby = "Site",
                                 fillvector = new_genera_legend) +
    theme(plot.margin = margin(r = -2)) +
    theme_cowplot(12)+
    theme(plot.title = element_text(hjust = 0.5))+
    theme(legend.position = "none")
  return(plot)
}

# Generate plots for each item
plots <- lapply(donors, generate_plots)
names(plots) <- paste(donors, "L6_lum", sep = "_")

# Access plots using list indexing
for (donor in donors) {
  assign(paste0(donor, "_L6_lum"), plots[[paste0(donor, "_L6_lum")]])
}


### Human Feces L6 Plots ---
generate_L6_human_feces_taxa_plots <- function(donor,filepath, titlestring,greppattern, fillvector, graphby){
 #input_data <- readr::read_csv(here("Donors-Analysis/taxa_barplots/Hoomins_level-6.csv"))
 
 input_data <- readr::read_csv(here(filepath))
 input_data <- as.data.frame(input_data)
 input_data <- input_data %>% filter(Donor_ID==donor)
 input_data <- input_data[,-c(99:116)]

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
 
  original_string <- filepath
  new_filepath<- sub("\\.csv$", "", original_string)
  new_filepath <- paste0(new_filepath,donor,".RDS")
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
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = factor(L2_lum$Site)
  }
 
  cols <- fillvector
  L2_lum$Site <- "F"
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = cols)+
    ylab("") +
    xlab("")+
    labs(fill="") +
    ggtitle(paste("Human")) +
    theme_cowplot(12) +
    theme(legend.position="none") +
    theme(axis.text.y=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}

generate_L6_SD_human_feces_taxa_plots <- function(filepath, titlestring,greppattern, fillvector, graphby){
  #input_data <- readr::read_csv(here("Donors-Analysis/taxa_barplots/Hoomins_level-6.csv"))
  
  input_data <- readr::read_csv(here("Humanized-Biogeography-Analysis/SD_Donor/taxa_barplots/HUM_SD_Donor_level-6.csv"))
  input_data <- as.data.frame(input_data)
  #input_data <- input_data %>% filter(Donor_ID==donor)
  #input_data <- input_data[,-c(99:116)]
  
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
  
  original_string <- filepath
  new_filepath<- sub("\\.csv$", "", original_string)
  new_filepath <- paste0(new_filepath,donor,".RDS")
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
  L2_lum<- tidyr::pivot_longer(L2_lum, -c(Taxa), values_to ="Value", names_to ="Site")
  L2_lum$Value <- L2_lum$Value * 100
  
  if({{graphby}} == "Site"){
    L2_lum$Site = factor(L2_lum$Site)
  }
  
  cols <- fillvector
  L2_lum$Site <- "F"
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    #scale_fill_paletteer_d("ggsci::category20_d3")+
    scale_fill_manual(values = genera_cols)+
    ylab("") +
    xlab("")+
    labs(fill="") +
    ggtitle(paste("Human")) +
    theme_cowplot(12) +
    theme(legend.position="none") +
    theme(axis.text.y=element_blank()) +
    theme(plot.title = element_text(hjust = 0.5))+
    #guides(fill=guide_legend(nrow=8, byrow=TRUE)) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  
}


# Access plots using list indexing
plots=list()
for(donor in human_donors){
plots[[donor]]<- generate_L6_human_feces_taxa_plots(donor = donor,
                                                   filepath = "Donors-Analysis/taxa_barplots/Hoomins_level-6.csv",
                                                           titlestring = "Human",
                                                           greppattern = ".*g__", 
                                                           graphby="Site",
                                                           fillvector = new_genera_legend)
}


# Call the function to generate plots
names(plots)
names(plots) <- paste(human_donors, "feces_L6", sep = "_")

# Access plots using list indexing
for (donor in human_donors) {
  assign(paste0(donor, "_feces_L6"), plots[[paste0(donor, "_feces_L6")]])
}

### Add HUM SD Gavage ---

A001_genera <- get_genera_from_plot("Humanized-Biogeography-Analysis/SD_Donor/taxa_barplots/HUM_SD_Donor_level-6.RDS")

A001_feces_L6 <- generate_L6_SD_human_feces_taxa_plots(filepath="Humanized-Biogeography-Analysis/SD_Donor/taxa_barplots/HUM_SD_Donor_level-6.csv",
                                                       titlestring = "Human",
                                                       greppattern = ".*g__", 
                                                       graphby="Site",
                                                       fillvector = new_genera_legend)



A001_L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/taxa_barplots/HUM_Gavage_Luminal_level-6.RDS",
                                            "(SD) A001", ".*g__",new_genera_legend, "Site") +
  theme(plot.margin = margin(r = -2)) +
  theme_cowplot(12)+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "none")

### Assemble Final Figure ---
dev.new()
plot_grid(A001_L6_lum, A001_feces_L6,
  A017_L6_lum, A017_feces_L6,
          A018_L6_lum, A018_feces_L6,
          A041_L6_lum, A041_feces_L6,
          A043_L6_lum, A043_feces_L6,
          A045_L6_lum, A045_feces_L6,
          A046_L6_lum, A046_feces_L6,
          A047_L6_lum, A047_feces_L6,
          A050_L6_lum, A050_feces_L6,
          A053_L6_lum, A053_feces_L6,
          A054_L6_lum, A054_feces_L6,
          A070_L6_lum, A070_feces_L6,
          A078_L6_lum, A078_feces_L6,
          A082_L6_lum, A082_feces_L6,
          rel_widths = c(1,0.25,
                         1,0.25,
                         1, 0.25,
                         1, 0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1,0.25,
                         1, 0.25,
                         1,0.25))

### Generate color legend ---
labels_all <- list()
for (donor_id in human_donors) {
  # Generate file path for the RDS file
  file_path <- paste0("Donors-Analysis/taxa_barplots/Hoomins_level-6", donor_id, ".RDS")
  
  # Get genera labels from the RDS file
  labels <- get_genera_from_plot(file_path)
  
  # Store the labels in the list
  labels_all[[donor_id]] <- labels
}
all_hum_labels <- unique(Reduce(union, labels_all))


# Loop through each donor ID
labels_all <- list()
for (donor_id in human_donors) {
  # Generate file path for the RDS file
  file_path <- paste0("Donors-Analysis/taxa_barplots/plot_by_donor/", donor_id, "_Lum_level-6.RDS")
  
  # Get genera labels from the RDS file
  labels <- get_genera_from_plot(file_path)
  
  # Store the labels in the list
  labels_all[[donor_id]] <- labels
}

# Combine labels from all donors
all_mouse_labels <- unique(Reduce(union, labels_all))

# combine labels from human and mice 
donors_legend <- unique(union(all_hum_labels,all_mouse_labels))
A001_feces_genera <- get_genera_from_plot("Humanized-Biogeography-Analysis/SD_Donor/taxa_barplots/HUM_SD_Donor_level-6.RDS")
A001_lum_genera <- get_genera_from_plot("Humanized-Biogeography-Analysis/taxa_barplots/Humanized_Luminal_level-6.RDS")

donors_legend <- unique(c(donors_legend, A001_feces_genera,A001_lum_genera))
names(genera_cols)


existingcols <- intersect(donors_legend,names(genera_cols))
old_legend <- genera_cols[names(genera_cols) %in% existingcols]
seecolor::print_color(old_legend,type="r")
missingcols <- setdiff(donors_legend,names(genera_cols))

cols1 <- paletteer_d("colorBlindness::Blue2DarkOrange18Steps", 18)
cols2 <- paletteer_d("colorBlindness::Green2Magenta16Steps", 12)
add_cols <- unique(c(cols1,cols2))
names(add_cols) <- missingcols

new_names <- union(names(old_legend),names(add_cols))
new_genera_legend<- union(old_legend,add_cols)
names(new_genera_legend)<- new_names

readr::write_rds(new_genera_legend, here("Donors-Analysis/taxa_barplots/donors_genera_cols.RDS"))

### Draw legend ---
genera_cols <- readRDS("Donors-Analysis/taxa_barplots/donors_genera_cols.RDS")
dummyplot<- as.data.frame(genera_cols)
dummyplot$dummyy <- seq(1,52,1)
dummyplot$dummyx <- seq(1,104,2)
dummyplot$Genus <- row.names(dummyplot)
L6_legend <-  ggplot(dummyplot, aes(x=dummyx,y=Genus,fill=Genus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=new_genera_legend,name="Genus Legend")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol=6, byrow=TRUE))+
  theme_cowplot(12)+
  theme(legend.spacing.y = unit(0.01, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(2, 11, 0, 0)) 
legend <- cowplot::get_legend(L6_legend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)

phyla_cols <- readRDS("global_phyla_cols.RDS")
dummyplot<- as.data.frame(phyla_cols)
dummyplot$dummyy <- seq(1,8,1)
dummyplot$dummyx <- seq(1,16,2)
dummyplot$Genus <- row.names(dummyplot)
L2_legend <-  ggplot(dummyplot, aes(x=dummyx,y=Genus,fill=Genus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=phyla_cols,name="Phylum Legend")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=4, byrow=TRUE))+
  theme_cowplot(16)+
  theme(legend.spacing.y = unit(1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=3, linetype="solid"), legend.margin = margin(10, 10, 100, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)

### Make Aggregated taxa barplots for main figure ---
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
