library(paletteer)
library(dplyr)
library(here)
library(rlang)
library(funrar)
library(ggplot2)
library(cowplot)

setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library(Microbiome.Biogeography)
setwd("/home/julianne/Documents/biogeography/")


here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S1_Shotgun_TaxaBarplots.R")

### Wrangle Genera Names --- 

## Wrangle the dataframes into something that can be used with existing functions -- 
file_path <- "Shotgun/taxa_barplots/export_groupby_Site_CS_SPF_BioGeo_Shotgun_L2/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

colnames(input_data) <- gsub("p__","",colnames(input_data))
write.csv(input_data,here("Shotgun/taxa_barplots/CS_Shotgun_L2.csv"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_UCLA_O_SPF_BioGeo_Shotgun_L2/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

colnames(input_data) <- gsub("p__","",colnames(input_data))
write.csv(input_data,here("Shotgun/taxa_barplots/UCLA_Shotgun_L2.csv"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_SPF_Gavage_BioGeo_Shotgun_L2/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

colnames(input_data) <- gsub("p__","",colnames(input_data))
write.csv(input_data,here("Shotgun/taxa_barplots/SPF_Shotgun_L2.csv"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_HUM_Gavage_BioGeo_Shotgun_L2/feature-table.tsv"
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

colnames(input_data) <- gsub("p__","",colnames(input_data))
write.csv(input_data,here("Shotgun/taxa_barplots/HUM_Shotgun_L2.csv"))

shotgun_process_taxonomy_data <- function(file_path){
input_data <- read.delim(here(file_path),sep="\t",header=TRUE,row.names = 1)
input_data <- as.data.frame(t(input_data))

taxa<-colnames(input_data)
taxa <- gsub("\\|",".",taxa)
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
  if (isFALSE(startsWith(genus[[i]],"GGB"))) {
    genus[[i]] = genus[[i]] 
  }
  else {
    genus[[i]] <- paste0(family[[i]][1]," (f)")   
  }
  i=i+1
}



colnames(input_data) <-as.character(genus)
colnames(input_data)
input_data <-  as.data.frame(t(apply(input_data,1, function(x) tapply(x,colnames(input_data),sum))))

input_data <- input_data[,!grepl("Mitochondria", colnames(input_data))] 
input_data <- input_data[,!grepl("Chloroplast", colnames(input_data))] 
return(input_data)
}

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_UCLA_O_SPF_BioGeo_Shotgun_L6/feature-table.tsv"
ucla_o_L6_barplot <- shotgun_process_taxonomy_data(file_path)
readr::write_rds(ucla_o_L6_barplot, here("Shotgun/taxa_barplots/UCLA_O_SPF_Shotgun_level-6.RDS"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_CS_SPF_BioGeo_Shotgun_L6/feature-table.tsv"
cs_L6_barplot <- shotgun_process_taxonomy_data(file_path)
readr::write_rds(cs_L6_barplot, here("Shotgun/taxa_barplots/CS_SPF_Shotgun_level-6.RDS"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_HUM_Gavage_BioGeo_Shotgun_L6/feature-table.tsv"
hum_L6_barplot <- shotgun_process_taxonomy_data(file_path)
readr::write_rds(hum_L6_barplot, here("Shotgun/taxa_barplots/HUM_Gavage_Shotgun_level-6.RDS"))

file_path <- "Shotgun/taxa_barplots/export_groupby_Site_SPF_Gavage_BioGeo_Shotgun_L6/feature-table.tsv"
spf_L6_barplot <- shotgun_process_taxonomy_data(file_path)
readr::write_rds(spf_L6_barplot, here("Shotgun/taxa_barplots/SPF_Gavage_Shotgun_level-6.RDS"))

## Generate a global genera color key using paletteer colors --
ucla_labels_lum <- get_genera_from_plot("Shotgun/taxa_barplots/UCLA_O_SPF_Shotgun_level-6.RDS")
cs_labels_lum <- get_genera_from_plot("Shotgun/taxa_barplots/CS_SPF_Shotgun_level-6.RDS")
hum_labels_lum <- get_genera_from_plot("Shotgun/taxa_barplots/HUM_Gavage_Shotgun_level-6.RDS")
spf_labels_lum <- get_genera_from_plot("Shotgun/taxa_barplots/SPF_Gavage_Shotgun_level-6.RDS")

all_labels <- unique(c(ucla_labels_lum,cs_labels_lum,hum_labels_lum,spf_labels_lum))
global_genera_cols <- readr::read_rds(here("global_genera_cols.RDS"))

#Use same colors for old legend and new colors for new Genera
genera_cols <- global_genera_cols[names(global_genera_cols) %in% all_labels]
other_genera <- all_labels[!all_labels %in% names(global_genera_cols)]

assign_cols <- paletteer_d("basetheme::royal", 10)
add_cols <- paletteer_d("basetheme::clean",10)
add_cols2 <- paletteer_d("calecopal::lupinus",5)
add_cols3 <- paletteer_d("lisa::C_M_Coolidge",5)
add_cols4 <- paletteer_d("lisa::FridaKahlo",5)
add_cols5 <- paletteer_d("beyonce::X6",6)
add_cols6 <- paletteer_d("beyonce::X39",2)
assign_cols <- unique(c(assign_cols, add_cols,
                 add_cols2, add_cols3,
                 add_cols4,add_cols5,
                 add_cols6))

#Match taxa to colors and then use in scale_fill_manual
names(assign_cols)<-other_genera
length(unique(c(genera_cols,other_genera)))
shotgun_genera <- c(genera_cols, assign_cols)
readr::write_rds(shotgun_genera,here("Shotgun/taxa_barplots/global_shotgun_genera_cols.RDS"))

## Generate phyla color key using paletteer colors --
hum_phyla <- get_phyla_from_plot("Shotgun/taxa_barplots/HUM_Shotgun_L2.csv")
cs_phyla <- get_phyla_from_plot("Shotgun/taxa_barplots/CS_Shotgun_L2.csv")
spf_phyla <- get_phyla_from_plot("Shotgun/taxa_barplots/SPF_Shotgun_L2.csv")
ucla_phyla <- get_phyla_from_plot("Shotgun/taxa_barplots/UCLA_Shotgun_L2.csv")

global_phyla <- c(hum_phyla, cs_phyla,
                  spf_phyla, ucla_phyla)
global_phyla <- unique(global_phyla)
colors<- unique(c(viridis::inferno(8),viridis::mako(4)))
colors <- c(colors[1:8], colors[10:12])
seecolor::print_color(colors)
names(colors) <- global_phyla
names(colors) <- c("Bacteria_unclassified","Bacteroidetes","Cyanobacteria",
                   "Deferribacteres","Candidatus_Saccharibacteria","Firmicutes",
                   "Proteobacteria","Verrucomicrobia","Candidatus_Kryptonia",
                   "Candidatus_Melainabacteria","Tenericutes")
readr::write_rds(colors, here("Shotgun/taxa_barplots/global_phyla_cols.RDS"))


## Draw plots --
phyla_cols<-colors
hum_shotgun_L2 <-generate_L2_taxa_plots("Shotgun/taxa_barplots/HUM_Shotgun_L2.csv", "HUM SD Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())

ucla_shotgun_L2 <-generate_L2_taxa_plots("Shotgun/taxa_barplots/UCLA_Shotgun_L2.csv", "UCLA O. SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")
  #theme(axis.text.y = element_blank())

cs_shotgun_L2 <-generate_L2_taxa_plots("Shotgun/taxa_barplots/CS_Shotgun_L2.csv", "CS SPF", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())

spf_shotgun_L2 <-generate_L2_taxa_plots("Shotgun/taxa_barplots/SPF_Shotgun_L2.csv", "SPF Gavage", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())


shotgun_genera <- readr::read_rds(here("Shotgun/taxa_barplots/global_shotgun_genera_cols.RDS"))


ucla_o_shotgun_L6 <- generate_L6_taxa_plots("Shotgun/taxa_barplots/UCLA_O_SPF_Shotgun_level-6.RDS",
                                       "UCLA O. SPF", ".*g__",shotgun_genera, "Site") +
  theme(legend.position = "none")
 # theme(axis.text.y = element_blank())
ucla_o_shotgun_L6

cs_shotgun_L6 <- generate_L6_taxa_plots("Shotgun/taxa_barplots/CS_SPF_Shotgun_level-6.RDS",
                                            "CS SPF", ".*g__",shotgun_genera, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())
cs_shotgun_L6

hum_shotgun_L6 <- generate_L6_taxa_plots("Shotgun/taxa_barplots/HUM_Gavage_Shotgun_level-6.RDS",
                                        "HUM SD Gavage", ".*g__",shotgun_genera, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())
hum_shotgun_L6

spf_shotgun_L6 <- generate_L6_taxa_plots("Shotgun/taxa_barplots/SPF_Gavage_Shotgun_level-6.RDS",
                                         "SPF Gavage", ".*g__",shotgun_genera, "Site") +
  theme(legend.position = "none")+
  theme(axis.text.y = element_blank())
spf_shotgun_L6

plot_grid(ucla_shotgun_L2,cs_shotgun_L2,spf_shotgun_L2,hum_shotgun_L2, nrow=1, labels=c("A"), label_size = 16)

plot_grid(ucla_o_shotgun_L6,cs_shotgun_L6,spf_shotgun_L6,hum_shotgun_L6, nrow=1, labels=c("B"), label_size=16)

## Plot Color Legends --
phyla_cols <- readRDS("Shotgun/taxa_barplots/global_phyla_cols.RDS")
dummyplot<- as.data.frame(phyla_cols)
dummyplot$dummyy <- seq(1,11,1)
dummyplot$dummyx <- seq(1,22,2)
dummyplot$Genus <- row.names(dummyplot)
L2_legend <-  ggplot(dummyplot, aes(x=dummyx,y=Genus,fill=Genus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=phyla_cols,name="Phylum Legend")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=2))+
  theme_cowplot(12)+
  theme(legend.text = element_text(
    margin = margin(r = 20, unit = "pt")))+
  #theme(legend.spacing.y = unit(0.2, 'cm'),legend.spacing.x = unit(0.2,'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=3, linetype="solid"), legend.margin = margin(10, 10, 10, 10)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)

genera_cols <- readRDS("Shotgun/taxa_barplots/global_shotgun_genera_cols.RDS")
dummyplot<- as.data.frame(genera_cols)
dummyplot$dummyy <- seq(1,62,1)
dummyplot$dummyx <- seq(1,124,2)
dummyplot$Genus <- row.names(dummyplot)
L6_legend <-  ggplot(dummyplot, aes(x=dummyx,y=Genus,fill=Genus))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values=genera_cols,name="Genus Legend")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(ncol=5, byrow=TRUE))+
  theme_cowplot(12)+
  theme(legend.text = element_text(
    margin = margin(r = 20, unit = "pt")))+
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(10, 10, 10, 10)) 
legend <- cowplot::get_legend(L6_legend)
grid.newpage()
dev.new(width=20, height=5)
grid.draw(legend)
