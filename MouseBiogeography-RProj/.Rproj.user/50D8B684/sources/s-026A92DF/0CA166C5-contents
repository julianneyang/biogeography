##Figure 6- Cedars SPF Gavage mice
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
library(gridExtra)
library(funrar)
library(here)
library(Microbiome.Biogeography)

?generate_GMM_heat_map_by_site
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Figure6_cowplot.R")

### Alpha Diversity ---
data<-read.csv("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv", header=TRUE, row.names=1)

generate_adiv_plots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  data$Type <- factor(data$Type, levels = c("Lum","Muc"))
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max)+
    theme(legend.position = "none")
  
}

#Merge metadata with alpha diversity 
data<-read.csv("Humanized-Biogeography-Analysis/alpha_diversity_Humanized.csv", header=TRUE, row.names=1)
metadata<- read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv")
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate
data$Microbiota
data <- data %>% filter(Microbiota=="Cedars_SPF")

#Shorten site names 
names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Type= revalue(data$Type, c("Luminal" = "Lum", "Mucosal" = "Muc"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
sapply(data,levels)

luminaldata<-filter(data, Type=="Lum")
mucosaldata<-filter(data, Type =="Muc")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

#Aggregate Plot by Sites Luminal Data
otus<- generate_adiv_plots(luminaldata, Site, observed_otus, 0, 450) + ggtitle ("Luminal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

pielou_e<- generate_adiv_plots(luminaldata, Site, pielou_e,0,1) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Sites Mucosal Data
otus_m <- generate_adiv_plots(mucosaldata, Site, observed_otus, 0, 450) + ggtitle ("Mucosal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_m <- generate_adiv_plots(mucosaldata, Site, pielou_e,0,1) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Type Colon Data by Site
otus_col <- generate_adiv_plots(colondata, Type, observed_otus, 0, 450)+ facet_grid(~Site) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_col <- generate_adiv_plots(colondata, Type, pielou_e, 0, 1)+ facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Type SI Data by site
otus_SI <- generate_adiv_plots(SIdata, Type, observed_otus, 0, 450) + facet_grid(~Site) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_SI <- generate_adiv_plots(SIdata, Type, pielou_e, 0, 1) + facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

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
  data$Microbiota <-factor(data$Microbiota, levels=c("Humanized", "Cedars_SPF"))
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
      theme(legend.position="top",legend.justification = "center")+ 
    #coord_fixed(ratio=1/2)+
      labs(title={{title}}) 
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

luminal <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - Lum.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Site_General",cols_general)+
  labs(title="Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mucosal <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - Muc.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Site_General",cols_general)+
  labs(title="Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LC.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Colon- Luminal", "Site",Colon_cols)+
  labs(title="Colon- Luminal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lsi <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - LSI.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "SI- Luminal", "Site",SI_cols) +
  labs(title="SI- Luminal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MC.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Colon- Mucosal", "Site",Colon_cols) +
  labs(title="Colon- Mucosal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
msi <-generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Site/Source RPCA - SPF - MSI.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "SI- Mucosal", "Site",SI_cols) +
  labs(title="SI- Mucosal") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

colon <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Colon.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  labs(title="Colon") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
SI <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - SI.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  labs(title="SI") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
DC <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Distal_Colon.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
PC <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Proximal_Colon.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
cec <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Cecum.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
ile <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Ileum.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
jej <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Jejunum.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
duo <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/Type/Source RPCA - SPF - Duodenum.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

### Taxa Barplots ---
assign_cols <- read_rds("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/assign_cols.RDS")

L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal_level-6.RDS", "Luminal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")
L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")
fig6f <- plot_grid(L6_lum,L6_muc, NULL, ncol=3,labels=c("F",""), rel_widths = c(1.5,1.5,1))

lummucphyla <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal/level-2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)
L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Luminal/level-2.csv", "Luminal Phyla", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")
L2_muc <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal/level-2.csv","Mucosal Phyla", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")
fig6e <- plot_grid(L2_lum,L2_muc, NULL, ncol=3,labels=c("E",""), rel_widths = c(2,2,1))

L2_col <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Col_LumMuc_L2.csv", "Colon Phyla", ".*p__", phyla_cols, "SiteTypeColon")+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L2_SI <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/SI_LumMuc_L2.csv", "SI Phyla", ".*p__", phyla_cols,"SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig6k <- plot_grid(L2_col,L2_SI, NULL, ncol=3,labels=c("M",""), rel_widths = c(2,2,1))

L6_col <-generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Col_LumMuc_L6.RDS","Colon Genera", ".*g__", assign_cols, "SiteTypeColon")+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L6_SI <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/SI_LumMuc_L6.RDS", "SI Genera", ".*g__", assign_cols, "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig6L <- plot_grid(L6_col,L6_SI, NULL, ncol=3,labels=c("K",""), rel_widths = c(1.5,1.5,1))


# Grab legend 
L6_muc_legend <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=25, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 0)) 
legend <- cowplot::get_legend(L6_muc_legend)
grid.newpage()
grid.draw(legend)
L2_muc_legend <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/SPF/barplots/Mucosal/level-2.csv","Mucosal Phyla",  ".*p__",phyla_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_muc_legend)
grid.newpage()
grid.draw(legend)

### Compile multi-panel figure: Longitudinal ---
fig6ab<-plot_grid(otus, otus_m, pielou_e,pielou_m, nrow=2, labels=c("A", "B", "",""))

left <- plot_grid(fig6ab,fig6e, nrow=2,rel_heights = c(2,1))

lum <-plot_grid(luminal,lsi,lc, nrow=3, align="hv")
muc <-plot_grid(mucosal, msi,mc, nrow=3, align="hv")
fig6c <- plot_grid(lum,labels = c("C"))
fig6d <- plot_grid(muc,labels = c("D"))
fig6cd <- plot_grid(fig6c,fig6d,ncol=2, align="hv")

dev.new(width=15, height=10)
right <- plot_grid(fig6cd,fig6f, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right, align = "hv")

### Compile multi-panel figure: Cross-sectional ---
fig6h <- plot_grid(otus_col, pielou_e_col, nrow=2)
fig6g<- plot_grid(otus_SI, pielou_e_SI, nrow=2) 

fig6gh <-plot_grid(otus_SI, otus_col, pielou_e_SI,pielou_e_col, SI, colon, nrow=3, labels=c("G", "H", "","", "I","J"))
dev.new(width=15, height=10)
fig6gh

fig6j <- plot_grid (cec,PC,DC, ncol=1, align ="hv")
fig6i<-plot_grid(duo,jej,ile, ncol=1, align ="hv")

fig6ij <- plot_grid(fig6i,fig6j,ncol=2, align="hv", labels=c("K", "L"))

left <- plot_grid(fig6gh,fig6k, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig6ij,fig6L, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)
fig6ghijkl <- plot_grid(left,right)

full_fig6<-plot_grid(fig6abcdef, fig6ghijkl, nrow=2, align="hv")
dev.new(width=15, height=10)
full_fig6