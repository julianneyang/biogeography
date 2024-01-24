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
## Extract taxa from Luminal ---
L2_lum<-readr::read_rds(here("Donors-Analysis/taxa_barplots/Mouse_Luminal/mouse_A017_level-6.RDS"))
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
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_lum <- unique(L2_lum$Taxa)
fecal_cols <- paletteer_d("ggsci::category20_d3", 14)
names(fecal_cols) <- labels_lum

## Extract taxa from Hoomins ---
L2_lum<-readr::read_rds(here("Donors-Analysis/taxa_barplots/Hoomins_level-6.RDS"))

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

## Generate a color key for ---
#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
length(labels_all)
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 13)
add_cols <- paletteer_d("dutchmasters::milkmaid",1)
assign_cols <- c(assign_cols,add_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
write_rds(assign_cols,"CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS")

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

## Merge Luminal and Mucosal together to generate Type Barplots ---
lum <- readRDS("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS")
muc <- readRDS("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS")

names(lum)[17]
names(muc)[17]

lum$Site_Type <- paste0(row.names(lum),"_","L")
muc$Site_Type <- paste0(row.names(muc),"_","M")

names(lum)==names(muc)
divergent <- c(setdiff(names(lum), names(muc)), setdiff(names(muc), names(lum)))

lum <- select(lum,-c(divergent))
muc <- select(muc, -c(divergent))

lummuc <- rbind(lum,muc)
row.names(lummuc)<-lummuc$Site_Type
lummuc <- select(lummuc, -Site_Type)

vsi <- c("Duodenum_L","Duodenum_M","Jejunum_L", "Jejunum_M","Ileum_L","Ileum_M")
vcol <-c("Cecum_L", "Cecum_M", "Proximal_Colon_L", "Proximal_Colon_M", "Distal_Colon_L","Distal_Colon_M")

lummucsi <- lummuc[row.names(lummuc) %in% vsi,]
lummuccol <- lummuc[row.names(lummuc) %in% vcol,]

readr::write_rds(lummucsi, "CS-Facility-Analysis/Taxa-Barplots/SI_LumMuc_L6.RDS")
readr::write_rds(lummuccol, "CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L6.RDS")

## Repeat at phylum level ---
lum <- read.csv("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-2.csv", row.names = 1, header = TRUE)
muc <-read.csv("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-2.csv", row.names = 1, header = TRUE)
lum$Site_Type <- paste0(row.names(lum),"_","L")
muc$Site_Type <- paste0(row.names(muc),"_","M")

names(lum)==names(muc)

lummuc <- rbind(lum,muc)
row.names(lummuc)<-lummuc$Site_Type
lummuc <- select(lummuc, -Site_Type)

vsi <- c("Duodenum_L","Duodenum_M","Jejunum_L", "Jejunum_M","Ileum_L","Ileum_M")
vcol <-c("Cecum_L", "Cecum_M", "Proximal_Colon_L", "Proximal_Colon_M", "Distal_Colon_L","Distal_Colon_M")

lummucsi <- lummuc[row.names(lummuc) %in% vsi,]
lummuccol <- lummuc[row.names(lummuc) %in% vcol,]

write.csv(lummucsi,"CS-Facility-Analysis/Taxa-Barplots/SI_LumMuc_L2.csv")
write.csv(lummuccol,"CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv")
