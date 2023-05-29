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
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/")

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
#L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
#L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
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
L2_lum <- generate_taxa_plots("Luminal_L2.csv", "Luminal", ".*p__") 
L2_muc <-generate_taxa_plots("Mucosal_L2.csv","Mucosal", ".*p__")
L2_col <- generate_taxa_plots("Colon_L2.csv", "Colon", ".*p__")
L2_SI <- generate_taxa_plots("SI_L2.csv", "SI", ".*p__")


dev.new(width=12, height=10)
plot_grid(L2_SI, L2_col, align="hv")

#handle the genera names 
input_data <- read.csv("Luminal_L6.csv", header=TRUE, row.names=1)
input_data <- read.csv("Mucosal_L6.csv", header=TRUE, row.names=1)
input_data <- read.csv("Colon_L6.csv", header=TRUE, row.names=1)
input_data <- read.csv("SI_L6.csv", header=TRUE, row.names=1)


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
write.csv(input_data, "Luminal_L6.csv")
write.csv(input_data, "Mucosal_L6.csv")
write.csv(input_data, "Colon_L6.csv")
write.csv(input_data, "SI_L6.csv")

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
#  L2_lum$Site = revalue(L2_lum$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
#  L2_lum$Site = factor(L2_lum$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  L2_lum$Site = revalue(L2_lum$Site, c("Luminal"="Lum", "Mucosal" = "Muc"))
  
  
  ggplot(data=L2_lum, aes(x=Site, y=Value, fill=Taxa)) +
    geom_bar(stat="identity")+
    scale_fill_paletteer_d(palette="colorBlindness::SteppedSequential5Steps") +
    #scale_fill_paletteer_d("ggthemes::manyeys")+
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

L6_lum <- generate_L6_taxa_plots("Luminal_L6.csv", "Luminal ( > 0.1% Relative Abundance)", ".*g__") 
L6_muc <-generate_L6_taxa_plots("Mucosal_L6.csv","Mucosal ( > 0.1% Relative Abundance)", ".*g__")
dev.new(width=12, height=10)
plot_grid(L2_lum,L6_lum,L2_muc,L6_muc, nrow=1, ncols=4,align="hv",labels = c("F","G", "H", "I"))
L6_col <-generate_L6_taxa_plots("Colon_L6.csv","Colon ( > 0.1% Relative Abundance)", ".*g__")
L6_SI <- generate_L6_taxa_plots("SI_L6.csv", "SI ( > 0.1% Relative Abundance)", ".*g__")
dev.new(width=12, height=10)
plot_grid(L2_SI, L6_SI, L2_col, L6_col, nrow=1, ncols=2,align="hv",labels = c("F","G", "H", "I"))

## Extract taxa from Luminal ---
temp<-read.csv("Luminal_L6.csv", row.names=1)
readr::write_rds(temp, "Luminal_L6.RDS")

temp<-read.csv("Mucosal_L6.csv", row.names=1)
readr::write_rds(temp, "Mucosal_L6.RDS")

L2_lum<-readRDS("Luminal_L6.RDS")
L2_lum<- as.matrix(L2_lum)
L2_lum<-funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa/6
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
L2_lum <- dplyr::select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- dplyr::select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_lum <- unique(L2_lum$Taxa)

## Extract taxa from Mucosal ---
L2_lum<-readRDS("Mucosal_L6.RDS")
L2_lum<- as.matrix(L2_lum)
L2_lum<- funrar::make_relative(L2_lum)
L2_lum<-as.data.frame(t(L2_lum))
toptaxa<- rowMeans(L2_lum)
L2_lum$averageRA <-toptaxa/6
L2_lum <- L2_lum %>% mutate(keeptaxa = ifelse(averageRA >0.001, row.names(L2_lum), "Other"))
L2_lum <- dplyr::select(L2_lum,-averageRA)

taxa<-L2_lum$keeptaxa
L2_lum <- dplyr::select(L2_lum,-keeptaxa)
L2_lum <- as.matrix(sapply(L2_lum,as.numeric))
L2_lum <- as.data.frame(prop.table(L2_lum,2))
taxa<-gsub(".*g__","",taxa )

L2_lum$Taxa <-taxa
labels_muc <- unique(L2_lum$Taxa)

## Generate a color key using paletteer colors ---
#Find out how many taxa need to be assigned colors 
labels_all <- union(labels_lum, labels_muc) 
#Generate that many colors 
library(paletteer)
assign_cols <- paletteer_d("ggsci::category20_d3", 16)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=labels_all
readr::write_rds(assign_cols,"assign_cols.RDS")

## Merge Luminal and Mucosal together to generate Type Barplots ---
lum <- readRDS("Luminal_L6.RDS")
muc <- readRDS("Mucosal_L6.RDS")
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

readr::write_rds(lummucsi, "SI_LumMuc_L6.RDS")
readr::write_rds(lummuccol, "Col_LumMuc_L6.RDS")

#Repeat at phylum level
lum <- read.csv("Luminal_L2.csv", row.names = 1, header = TRUE)
muc <-read.csv("Mucosal_L2.csv", row.names = 1, header = TRUE)

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

write.csv(lummucsi,"SI_LumMuc_L2.csv")
write.csv(lummuccol,"Col_LumMuc_L2.csv")

phyla<-gsub(".*p__","",names(lum) )
library(viridis)
library(seecolor)
assign_cols <-viridis(n = 8, option="B")
seecolor::print_color(assign_cols)
#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)=phyla
print(assign_cols)
readr::write_rds(assign_cols,"global_phyla_cols.RDS")
