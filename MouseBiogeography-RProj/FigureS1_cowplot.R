##Figure S1 - WT Validation Cohort  
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
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(plyr)
library(grid)
library(funrar)
library(here)

here::i_am("MouseBiogeography-RProj/FigureS1_cowplot.R")

### Alpha- Diversity 
generate_adiv_plots <- function(input_data, input_metadata, X, Y, facetvariable, min, max){
  #read in files
  metadata<- read.csv(input_metadata)
  data<-read.csv(input_data)
  metadata$SampleID <-gsub("-",".",metadata$SampleID)
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #Shorten site names 
  data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
  data$Type= revalue(data$Type, c("Luminal" = "Lum", "Mucosal" = "Muc"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  data$Type <- factor(data$Type, levels = c("Lum","Muc"))
  
  shannon <- ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{facetvariable}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=2,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max) +
    theme(legend.position = "none")
  
}

#Aggregate Plot by Sites Mucosal Data
otus<- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, observed_otus, Site,0,450)+
  ggtitle ("Mucosal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e <- generate_adiv_plots("ImmDef-Mouse-Biogeography-Analysis/alpha_diversity_WTCohort.csv", "ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv", Site, pielou_e, Site,0,1.2)+
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


### Beta-Diversity ---
generate_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
  data<-read.csv(ordination_file, header=FALSE)
  metadata <- read.table(metadata, sep="\t", header=TRUE)
  
  #store PC1 and Pc2
  PC1<-data[5,1]
  PC1 <-round(as.numeric(PC1)*100, digits=1)
  PC2<-data[5,2]
  PC2 <-round(as.numeric(PC2)*100, digits=1)
  PC1 <-as.character(PC1)
  str_PC1<-paste("PC 1 (", PC1,"%)")
  str_PC2<-paste("PC 2 (", PC2, "%)")
  
  #edit dataframe
  data<-data[,1:4]
  data <- slice(data, 1:(n() - 4))     # Apply slice & n functions
  data<-data[-c(1,2,3,4,5,6,7,8,9),]
  
  #rename columns
  names(data)[names(data) == "V1"] <- "SampleID" 
  names(data)[names(data)=="V2"] <- "PC1" 
  names(data)[names(data)=="V3"] <- "PC2"
  names(data)[names(data)=="V4"] <- "PC3"
  # data$SampleID<-gsub(".","",data$SampleID)
  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #declare factors
  data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))
  #data$Microbiota <-factor(data$Microbiota, levels=c("Humanized", "Cedars_SPF"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  colorvariable <- as.factor(colorvariable)
  if(colorvariable =="Site"){
    p<- ggplot(data, aes(x=PC1, y=PC2, colour=Site)) + 
      #geom_point(size=3) + 
      geom_point(aes(fill=Site), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      #scale_color_viridis_d()+
      xlab(str_PC1) +
      ylab(str_PC2) +
      theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    #coord_fixed(ratio=1/2)+
    #labs(title= paste0({{title}}, " RPCA")) 
    p
  }
  else if (colorvariable =="Type"){
    data$Type = revalue(data$Type, c("Luminal"="Lum", "Mucosal"="Muc"))
    p<- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
      #geom_point(size=3) + 
      geom_point(aes(fill=Type), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      #scale_color_viridis_d()+
      xlab(str_PC1) +
      ylab(str_PC2) +
      theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    #coord_fixed(ratio=1/2)+
    #labs(title= paste0({{title}}, " RPCA")) 
    p
  }
  else if (colorvariable =="Site_General"){
    p<- ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
      geom_point(aes(fill=Site_General), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      xlab(str_PC1) +
      ylab(str_PC2) +
      theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    p
  }
}

Colon_cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
SI_cols<- c("Duo" = "red", "Jej" = "gold", "Ile" = "green")
all_cols <-  c("Duo" = "red", "Jej" = "gold", "Ile" = "green","Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Microbiota_cols <-c("Humanized"="purple", "Cedars_SPF" = "turquoise")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")

mucosal <-generate_pcoA_plots("ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal.csv",
                         "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.tsv.txt", "Mucosal", "Site_General", cols_general) +
  labs(title="Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mc <-generate_pcoA_plots("ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal_Colon.csv",
                    "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.tsv.txt", "Colon- Mucosal", "Site", Colon_cols) +
  labs(title="Colon- Mucosal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
msi <-generate_pcoA_plots("ImmDef-Mouse-Biogeography-Analysis/RPCA/Pre-Combat/WT Cohort Site RPCA - Mucosal_SI.csv", 
                    "ImmDef-Mouse-Biogeography-Analysis/RPCA/WTCohort-Metadata.tsv.txt", "SI- Mucosal", "Site",SI_cols) +
  labs(title="SI- Mucosal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

### Taxa Barplots ---
lummucphyla <- read.csv("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

assign_cols <- readRDS("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/assign_cols.RDS")
L2_muc <-generate_L2_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "Mucosal Phyla", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")

figS1c<- plot_grid(L2_muc, NULL, ncol=2, rel_widths = c(2,1))

L6_muc <-generate_L6_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__",assign_cols, "Site") +
  theme(legend.position = "none")

figS1d <- plot_grid(L6_muc, NULL, ncol=2, rel_widths = c(1.5,1))

# Grab legend 
L6_muc_legend <-generate_L6_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__",assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_muc_legend)
  grid.newpage()
  grid.draw(legend)
L2_muc_legend <-generate_L2_taxa_plots("ImmDef-Mouse-Biogeography-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "Mucosal Phyla", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_muc_legend)
grid.newpage()
grid.draw(legend)


### Compile multi-panel figure ---

figS1a<-plot_grid(otus, pielou_e, nrow=2, labels=c("A", ""))

muc <-plot_grid(mucosal, msi,mc, nrow=3, align="hv")
figS1b <- plot_grid(muc,labels = c("B"))

figS1cd <-plot_grid(figS1c, figS1d, labels=c("C",""),nrow=2)
dev.new(width=15, height=10)
plot_grid(figS1a,figS1b,figS1cd, align="hv", ncol=3)