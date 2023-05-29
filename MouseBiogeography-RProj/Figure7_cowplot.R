##Figure 7 - CS Facility data 
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
library(ggsignif)
library(ggbeeswarm)
library(ggpubr)
library(plyr)
library(grid)
library(funrar)

here::i_am("MouseBiogeography-RProj/Figure7_cowplot.R")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
### Alpha Diversity ---

data <- readRDS("CS-Facility-Analysis/alpha_diversity_CS_Facility.RDS")
metadata<- read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv")
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

#filter data set 
luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")

generate_adiv_plots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  
  #Shorten site names 
  data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
  data$Type= revalue(data$Type, c("Luminal" = "Lum", "Mucosal" = "Muc"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
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
#Aggregate Plot by Sites Luminal Data
otus_l <- generate_adiv_plots(luminaldata, Site, observed_otus, 0, 600)+
  ggtitle ("Luminal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_l <- generate_adiv_plots(luminaldata, Site, pielou_e, 0, 1)+
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Sites Mucosal Data
otus_m<- generate_adiv_plots(mucosaldata, Site, observed_otus,0,600)+
  ggtitle ("Mucosal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_m<- generate_adiv_plots(mucosaldata, Site, pielou_e, 0, 1)+
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


#Aggregate Plot by Type Colon Data by Site
otus_col <- generate_adiv_plots(colondata, Type, observed_otus, 0, 600)+ facet_grid(~Site) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_col <- generate_adiv_plots(colondata, Type, pielou_e, 0, 1)+ facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Type SI Data by site
otus_SI <- generate_adiv_plots(SIdata, Type, observed_otus, 0, 600) + facet_grid(~Site) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_SI <- generate_adiv_plots(SIdata, Type, pielou_e, 0, 1) + facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)


## Beta-Diversity ---

generate_pcoA_plots <- function(ordination_file, metadata, title, colorvariable,colorvector){
  data<-read.csv(ordination_file, header=FALSE)
  metadata <- read.csv(metadata, header=TRUE)
  #metadata$SampleID <- metadata$X.SampleID
  
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

  #append metadata
  intermediate<- (merge(data, metadata, by = 'SampleID'))
  data<- intermediate
  
  #declare factors
  data$Site_General<-factor(data$Site_General, levels=c("SI", "Colon"))
  data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
  
  if(colorvariable =="Site"){
    p<- ggplot(data, aes(x=PC1, y=PC2, colour=Site)) + 
      geom_point(aes(fill=Site), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      xlab(str_PC1) +
      ylab(str_PC2) +
      theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
    p
  }
  else if (colorvariable =="Type"){
    data$Type = revalue(data$Type, c("Luminal"="Lum", "Mucosal"="Muc"))
    p<- ggplot(data, aes(x=PC1, y=PC2, colour=Type)) + 
      geom_point(aes(fill=Type), colour="black", pch=21, size=3) +
      scale_fill_manual(name="", values={{colorvector}}) +
      xlab(str_PC1) +
      ylab(str_PC2) +
      theme_cowplot(12)+
      theme(legend.position="top",legend.justification = "center") 
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
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
SI_cols <- c("Duo" = "firebrick", "Jej"="gold", "Ile" = "forestgreen")
Colon_cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")

luminal <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumCol", "Site_General",cols_general)+
  labs(title="Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mucosal <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumCol", "Site_General",cols_general)+
  labs(title="Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lc <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumCol", "Site",Colon_cols)+
  labs(title="Colon- Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lsi <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Luminal_SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility LumSI", "Site",SI_cols) +
  labs(title="SI- Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mc <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility MucCol", "Site",Colon_cols)+
  labs(title="Colon- Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
msi <-generate_pcoA_plots("CS-Facility-Analysis/RPCA/Site_RPCA/CS_Facility_Site_RPCA - Mucosal_SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility MucSI", "Site",SI_cols) +
  labs(title="SI- Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

colon <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility DC", "Type", Type_cols) +
  labs(title="Colon") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
SI <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - SI.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility PC", "Type", Type_cols) +
  labs(title="SI") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
DC <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Distal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility DC", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
PC <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Proximal_Colon.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility PC", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
cec <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Cecum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Cec", "Type", Type_cols) +
  facet_grid(Site~.) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
ile <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Ileum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Ile", "Type", Type_cols)+
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
jej <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Jejunum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Jej", "Type", Type_cols)+
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
duo <- generate_pcoA_plots("CS-Facility-Analysis/RPCA/Type_RPCA/CS_Facility_Beta_Diversity - Duodenum.csv", "CS-Facility-Analysis/CS_Facility_Metadata.csv", "CS_Facility Duo", "Type", Type_cols)+
  facet_grid(Site~.) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

### Taxa Barplots ---
lummucphyla <- read.csv("CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

assign_cols <- readRDS("CS-Facility-Analysis/Taxa-Barplots/assign_cols.RDS")
L2_lum <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-2.csv", "Luminal Phyla", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")
L2_muc <-generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-2.csv", "Mucosal Phyla", ".*p__", phyla_cols, "Site") +
  theme(legend.position = "none")

fig7e <- plot_grid(L2_lum,L2_muc, NULL, ncol=3,labels=c("E",""), rel_widths = c(2,2,1))

L6_lum <- Microbiome.Biogeography::generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Luminal_level-6.RDS", "Luminal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")
L6_muc <-Microbiome.Biogeography::generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS","Mucosal Genera", ".*g__",assign_cols, "Site") +
  theme(legend.position = "none")

fig7f <- plot_grid(L6_lum,L6_muc, NULL, ncol=3,labels=c("F",""), rel_widths = c(1.5,1.5,1))

lummucphyla <- read.csv("CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)
L6_col <-generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L6.RDS","Colon Genera", ".*g__",assign_cols, "SiteTypeColon") +
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L6_SI <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/SI_LumMuc_L6.RDS", "SI Genera", ".*g__",assign_cols, "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig7L <- plot_grid(L6_SI,L6_col, NULL, ncol=3,labels=c("L",""), rel_widths = c(1.5,1.5,1))

L2_col <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv", "Colon Phyla", ".*p__", phyla_cols, "SiteTypeColon") +
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L2_SI <- generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/SI_LumMuc_L2.csv", "SI Phyla", ".*p__", phyla_cols,"SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

fig7k <- plot_grid(L2_SI,L2_col, NULL, ncol=3,labels=c("K",""), rel_widths = c(2,2,1))

# Format legends 

L6_legend <- generate_L6_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-6.RDS","Mucosal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid.newpage()
grid.draw(legend)
L2_legend <-generate_L2_taxa_plots("CS-Facility-Analysis/Taxa-Barplots/Mucosal_level-2.csv","Mucosal Phyla", ".*p__", phyla_cols,"Site")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid.newpage()
grid.draw(legend)

### Compile multi-panel figure ---

#Longitudinal
fig7ab<-plot_grid(otus_l, otus_m, pielou_e_l,pielou_e_m, nrow=2, labels=c("A", "B", "",""))

lum <-plot_grid(luminal,lsi,lc, nrow=3, align="hv")
muc <-plot_grid(mucosal, msi,mc, nrow=3, align="hv")
fig7cd<-plot_grid(lum, muc, ncol=2,  align = 'hv',labels = c("C","D"))


left <- plot_grid(fig7ab,fig7e, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig7cd,fig7f, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)

fig7abcdef <- plot_grid(left,right)

#Cross-Sectional
fig7gh<-plot_grid(otus_SI, otus_col, pielou_e_SI,pielou_e_col, SI, colon, nrow=3, labels=c("A", "B", "","", "C","D"))

fig7j <- plot_grid (cec,PC,DC, ncol=1, align ="hv")
fig7i<-plot_grid(duo,jej,ile, ncol=1, align ="hv")
fig7ij <- plot_grid(fig7i,fig7j,ncol=2, align="hv", labels=c("I", "J"))

left <- plot_grid(fig7gh,fig7k, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig7ij,fig7L, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)

fig7ghijkl<- plot_grid(left,right)

dev.new(width=15, height=10)
plot_grid(fig7abcdef, fig7ghijkl, align="hv", nrow=2)