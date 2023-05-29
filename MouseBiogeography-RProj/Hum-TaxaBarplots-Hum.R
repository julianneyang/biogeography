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
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/")

generate_taxa_plots <- function(input_data, titlestring,greppattern){
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
  L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  #L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  L2_lum$Value <- L2_lum$Value * 100
  #L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  #L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_viridis_d(option="B")+
    theme(legend.position = "right")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="")+
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=4, byrow=TRUE))
  
  
}
L2_lum <- generate_taxa_plots("Luminal/level-2.csv", "SPF Luminal", ".*p__") 
L2_muc <-generate_taxa_plots("Mucosal/level-2.csv","SPF Mucosal", ".*p__")

dev.new(width=12, height=10)
plot_grid(L2_lum, L2_muc, align="hv")

L2_col <- generate_taxa_plots("Colon/level-2.csv", "Colon", ".*p__")
L2_SI <- generate_taxa_plots("SI/level-2.csv", "SI", ".*p__")

dev.new(width=12, height=10)
plot_grid(L2_SI, L2_col, align="hv")

#handle the genera names 
input_data <- read.csv("Luminal/level-6.csv", header=TRUE, row.names=1)
input_data <- read.csv("Mucosal/level-6.csv", header=TRUE, row.names=1)
input_data <- read.csv("Colon/level-6.csv", header=TRUE, row.names=1)
input_data <- read.csv("SI/level-6.csv", header=TRUE, row.names=1)


taxa<-colnames(input_data)
colnames <- strsplit(taxa, ".f__")

family=new_list(77)
i=1
for (i in 1:length(colnames)) {
  family[i] <- colnames[[i]][2]
  i=i+1
}

family<-unlist(family)
family <- strsplit(family, ".g__")

genus =new_list(77)
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
write_rds(input_data, "Luminal/level-6.RDS")
write_rds(input_data, "Mucosal/level-6.RDS")
write_rds(input_data, "Colon/level-6.RDS")
write_rds(input_data, "SI/level-6.RDS")

generate_L6_taxa_plots <- function(input_data, titlestring,greppattern){
  titlestring<-c(titlestring)
  L2_lum<-read.csv(input_data, row.names=1,header=TRUE)
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
  L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  #L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  
  
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    #scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d(palette="dutchmasters::milkmaid") +
    #scale_fill_paletteer_d("tvthemes::rickAndMorty")+
    scale_fill_paletteer_d("ggsci::category20_d3")+
    theme(legend.position = "none")+
    theme_cowplot(12) +
    ylab("% Relative Abundance") +
    xlab("")+
    labs(fill="") +
    ggtitle(titlestring) +
    theme(legend.position="top") +
    theme(plot.title = element_text(hjust = 0.5))+
    guides(fill=guide_legend(nrow=8, byrow=TRUE))
  
  
}

## Extract taxa from Luminal ---
L2_lum<-readRDS("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-6.RDS")
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

## Extract taxa from Mucosal ---
L2_lum<-readRDS("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS")
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
#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
length(labels_all)
#Generate that many colors 
assign_cols <- paletteer_d("ggsci::category20_d3", 20)
add_cols <- paletteer_d("dutchmasters::milkmaid",1)
assign_cols <- c(assign_cols,add_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
write_rds(assign_cols,"Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/assign_cols.RDS")

L6_lum_hum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-6.RDS", "Luminal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") 
L6_muc_hum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS","Mucosal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site")
dev.new(width=12, height=10)
plot_grid(L6_lum_hum,L6_muc_hum, align="hv")

## Merge Luminal and Mucosal together to generate Type Barplots ---
lum <- readRDS("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-6.RDS")
muc <- readRDS("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS")
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

readr::write_rds(lummucsi, "Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L6.RDS")
readr::write_rds(lummuccol, "Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Col_LumMuc_L6.RDS")

#Repeat at phylum level
lum <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-2.csv", row.names = 1, header = TRUE)
muc <-read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-2.csv", row.names = 1, header = TRUE)

names(lum)==names(muc)
lum$Site_Type <- paste0(row.names(lum),"_","L")
muc$Site_Type <- paste0(row.names(muc),"_","M")


lummuc <- rbind(lum,muc)
row.names(lummuc)<-lummuc$Site_Type
lummuc <- select(lummuc, -Site_Type)

vsi <- c("Duodenum_L","Duodenum_M","Jejunum_L", "Jejunum_M","Ileum_L","Ileum_M")
vcol <-c("Cecum_L", "Cecum_M", "Proximal_Colon_L", "Proximal_Colon_M", "Distal_Colon_L","Distal_Colon_M")

lummucsi <- lummuc[row.names(lummuc) %in% vsi,]
lummuccol <- lummuc[row.names(lummuc) %in% vcol,]

write.csv(lummucsi,"Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L2.csv")
write.csv(lummuccol,"Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Col_LumMuc_L2.csv")

