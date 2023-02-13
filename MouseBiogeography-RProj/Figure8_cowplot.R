### generate Multi-panel figures for the Manuscript 
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(ggbeeswarm)
library(ggpubr)
library(plyr)
library(grid)
library(funrar)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")

### Figure 8ab and 8gh---

here::i_am("MouseBiogeography-RProj/Figure8_cowplot.R")
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
data <- data %>% filter(Microbiota=="Humanized")

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
metadata <- read.table("Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", sep="\t", header=TRUE)
names(metadata)
Colon_cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
SI_cols<- c("Duo" = "red", "Jej" = "gold", "Ile" = "green")
all_cols <-  c("Duo" = "red", "Jej" = "gold", "Ile" = "green","Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
cols_general <- c("SI" = "#F8766D","Colon" ="#00BFC4")
Microbiota_cols <-c("Humanized"="purple", "Cedars_SPF" = "turquoise")
Type_cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")

luminal <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - Luminal.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Site_General",cols_general)+
  labs(title="Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mucosal <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - Mucosal.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Site_General",cols_general)+
  labs(title="Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LC.csv", 
                          "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized LumCol", "Site",Colon_cols) +
  labs(title="Colon- Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
lsi <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - LSI.csv", 
                           "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized LumSI", "Site",SI_cols)+
  labs(title="SI- Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
mc <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - MC.csv", 
                          "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized MucCol", "Site",Colon_cols) +
  labs(title="Colon- Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
msi <-generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Site/Source RPCA -Humanized - MSI.csv", 
                          "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized MucSI", "Site",SI_cols) +
  labs(title="SI- Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))




colon <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Colon.csv", 
                             "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized Colon", "Type", Type_cols) +
  labs(title="Colon") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
SI <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - SI.csv", 
                          "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "Humanized SI", "Type", Type_cols) +
  labs(title="SI") +
  facet_grid(Site_General~.)+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

DC <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - DC.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
PC <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - PC.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
cec <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Cec.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
ile <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Ile.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
jej <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Jej.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))
duo <- generate_pcoA_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/Type/Source RPCA -Humanized - Duo.csv", "Humanized-Biogeography-Analysis/Humanized Metadata.tsv.txt", "", "Type", Type_cols) +
  facet_grid(Site~.) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))



### Taxa Barplots ---
lummucphyla <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L2.csv",header=TRUE,row.names=1)
lummucphyla <- gsub(".*p__","",names(lummucphyla))
print(lummucphyla)
phyla_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/global_phyla_cols.RDS")
phyla_cols <- phyla_cols[names(phyla_cols) %in% lummucphyla]
print(phyla_cols)

assign_cols <- read_rds("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/assign_cols.RDS")
print(assign_cols)
seecolor::print_color(assign_cols)

L6_lum <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-6.RDS", "Luminal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")
L6_muc <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS","Mucosal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")

fig8f <- plot_grid(L6_lum,L6_muc, NULL, ncol=3,labels=c("F",""), rel_widths = c(1.5,1.5,1))

L2_lum <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Luminal/level-2.csv", "Luminal Phyla", ".*p__", phyla_cols,"Site") +
  theme(legend.position = "none")
L2_muc <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-2.csv","Mucosal Phyla", ".*p__", phyla_cols, "Site")+
  theme(legend.position = "none")

fig8e <- plot_grid(L2_lum,L2_muc, NULL, ncol=3,labels=c("E",""), rel_widths = c(2,2,1))

L2_col <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Col_LumMuc_L2.csv", "Colon Phyla", ".*p__", phyla_cols, "SiteTypeColon")+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L2_SI <- generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L2.csv", "SI Phyla", ".*p__", phyla_cols, "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig8k <- plot_grid(L2_SI,L2_col, NULL, ncol=3,labels=c("K",""), rel_widths = c(2,2,1))

L6_col <-generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Col_LumMuc_L6.RDS","Colon Genera", ".*g__", assign_cols, "SiteTypeColon")+
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L6_SI <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/SI_LumMuc_L6.RDS", "SI Genera", ".*g__", assign_cols, "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig8L <- plot_grid(L6_SI,L6_col, NULL, ncol=3,labels=c("L",""), rel_widths = c(1.5,1.5,1))

# Grab legend 
fig8f_legend <- generate_L6_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Mucosal/level-6.RDS","Mucosal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(fig8f_legend)
grid.newpage()
grid.draw(legend)
L2_muc_legend <-generate_L2_taxa_plots("Humanized-Biogeography-Analysis/Source RPCA/Hum/barplots/Col_LumMuc_L2.csv", "Colon Phyla", ".*p__", phyla_cols, "SiteTypeColon")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_muc_legend)
grid.newpage()
grid.draw(legend)



### Compile multi-panel figure ---
#Longitudinal
fig8ab<- plot_grid(otus, otus_m, pielou_e,pielou_m, cols=2, rows=2, align = 'hv',labels = c("A","B","",""))

lum <-plot_grid(luminal,lsi,lc, nrow=3, align="hv")
muc <-plot_grid(mucosal, msi,mc, nrow=3, align="hv")
fig8cd<-plot_grid(lum, muc, ncol=2,  align = 'hv',labels = c("C","D"))


left <- plot_grid(fig8ab,fig8e, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig8cd,fig8f, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)

#Cross-Sectional
fig8gh<-plot_grid(otus_SI, otus_col, pielou_e_SI,pielou_e_col, SI, colon, nrow=3, labels=c("G", "H", "","", "I","J"))

fig8j <- plot_grid (cec,PC,DC, ncol=1, align ="hv")
fig8i<-plot_grid(duo,jej,ile, ncol=1, align ="hv")
fig8ij <- plot_grid(fig7i,fig7j,ncol=2, align="hv", labels=c("K", "L"))

left <- plot_grid(fig8gh,fig8k, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig8ij,fig8L, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)
#compile fig 8ab 
dev.new(width=15, height=10)
#fig8abcd <- plot_grid(fig8ab,fig8cd)
plot_grid(fig8ab, fig8cd,fig8e, fig8f, align="hv",ncol=2,nrow=2)
## Make Fig 8 GHIJ ---
dev.new(width=15, height=10)
plot_grid(fig8gh,fig8ij,ncol=2, align="hv")
# GHIJKL 
dev.new(width=15, height=10)
plot_grid(otus_col, otus_SI, pielou_e_col,pielou_e_SI, nrow=2, labels=c("G", "H", "",""))
fig8gh<-plot_grid(otus_col, otus_SI, pielou_e_col,pielou_e_SI, nrow=2, labels=c("G", "H", "",""))

dev.new(width=15, height=10)
plot_grid(fig8gh,fig8k, nrow=2,rel_heights = c(2,1))
left <- plot_grid(fig8gh,fig8k, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(fig8ij,fig8L, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig8ij,fig8L, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)
fig8ghijkl <- plot_grid(left,right)


fig8j <- plot_grid (cec,PC,DC, ncol=1, align ="hv")
fig8i<-plot_grid(duo,jej,ile, ncol=1, align ="hv")

fig8ij <- plot_grid(fig8i,fig8j,ncol=2, align="hv", labels=c("I", "J"))


dev.new(width=15, height=10)
plot_grid(otus_col, pielou_e_col, nrow=2)
fig8h <- plot_grid(otus_col, pielou_e_col, nrow=2)

dev.new(width=15, height=10)
fig8g<- plot_grid(otus_SI, pielou_e_SI, nrow=2) 

fig8cd<-plot_grid(lum, muc, ncol=2,  align = 'hv',labels = c("C","D"))
dev.new(width=15, height=10)
fig8cd


lum <-plot_grid(lsi,lc, nrow=2, align="hv")
muc <-plot_grid(msi,mc, nrow=2, align="hv")

###Heatmaps: Figs 1fg
#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maaslin2 Type Genus Level/")
duodenum<-read.table("L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
a<-duodenum_significant$feature
jejunum<-read.table("L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
b<-jejunum_significant$feature
ileum<-read.table("L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
ileum_significant<-filter(ileum, metadata=="Type" & value=="Mucosal" &qval<0.05)
c<-ileum_significant$feature
cecum<-read.table("L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
d<-cecum_significant$feature  
pc<-read.table("L6-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
e<-pc_significant$feature  
DC<-read.table("L6-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID/significant_results.tsv", header=TRUE)
DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

#Query the target vector against all_results.tsv for each site
duodenum<-read.table("L6-LumRef-CLR-Duodenum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
duodenum_all$Site<- "Duodenum"
jejunum<-read.table("L6-LumRef-CLR-Jejunum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
jejunum_all$Site<- "Jejunum"
ileum<-read.table("L6-LumRef-CLR-Ileum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
ileum_all<-ileum_all[match(target,ileum_all$feature),]
ileum_all$Site<- "Ileum"
cecum<-read.table("L6-LumRef-CLR-Cecum-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
cecum_all<-cecum_all[match(target,cecum_all$feature),]
cecum_all$Site<- "Cecum"
pc<-read.table("L6-LumRef-CLR-ProximalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
pc_all<-pc_all[match(target,pc_all$feature),]
pc_all$Site<- "Proximal_Colon"
DC<-read.table("L6-LumRef-CLR-DistalColon-ComBat-SeqRunLineSexType-1-MsID/all_results.tsv", header=TRUE)
DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
DC_all<-DC_all[match(target,DC_all$feature),]
DC_all$Site<- "Distal_Colon"

duojej<-rbind(duodenum_all,jejunum_all)
ilecec<-rbind(ileum_all, cecum_all)
pcdc<-rbind(pc_all,DC_all)
duojejilecec<-rbind(duojej,ilecec)
duojejilececpcdc<-rbind(duojejilecec,pcdc)

#write.csv(duojejilececpcdc, "Genus-TYPE-Heatmap.csv") 
#from here make sure all NA rows are filled with feature name corresponding to NA via copy paste
#remove across six sites all GMM that failed to converge
duojejilececpcdc<-read.csv("Genus-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("genus_taxonomy.csv", header=TRUE)
data<- (merge(gmm_heatmap_final, annotation, by = 'feature'))
data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")

qval<-data$qval
asterisk<-c("")
for (item in qval){
  if (item < 0.05){
    asterisk<-c(asterisk,"*")
  }
  else {
    asterisk<-c(asterisk,"")
  }
}
asterisk<-asterisk[-1]
data$asterisk<-asterisk
data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 
y = tapply(data$coef_d, data$Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Genus= factor(as.character(data$Genus), levels = names(y))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
ggplotdata<-data
cols=viridis(8)

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

#assign colors to phyla
phyla<-ggplotdata$Phylum
ggplotdata$Phylum <- factor(ggplotdata$Phylum)
tick_colors<-c("")
i=1
phyla[1]
for (item in phyla){
  if (phyla[i]=="Proteobacteria"){
    tick_colors<-c(tick_colors,"firebrick")
  }
  else if (phyla[i]=="Actinobacteria") {
    tick_colors<-c(tick_colors,"purple")
  }
  else if (phyla[i]=="Bacteroidetes"){
    tick_colors<-c(tick_colors, "blue")
  }
  else if (phyla[i]=="Firmicutes"){
    tick_colors<-c(tick_colors, "forestgreen")
  }
  else if (phyla[i]=="Verrucomicrobia"){
    tick_colors<-c(tick_colors, "cyan")
  }
  i=i+1
}
tick_colors<-tick_colors[-1]
ggplotdata$tick_colors <- tick_colors
write.csv(ggplotdata, ("Phyla_Colors_MucvsLum.csv"))

### make graph after reordering tick colors 
tick_colors<-read.csv("Fig 3c Muc vs Lum - Phyla_Colors_MucvsLum.csv", header=TRUE)
tick_colors<-tick_colors$tick_colors
tick_colors<-rev(tick_colors)
g1 <- ggplot(ggplotdata, aes(x = Site, y=Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="right",legend.justification = "center") +
  xlab("")+
  ylab("") +
  theme(axis.text.y = element_text(colour = tick_colors))+
  guides(fill=guide_colourbar(title="",label=TRUE, barheight = 15))
dev.new(width=15, height=10)
f3d <- g1 +ggtitle("Mucosal vs. Luminal") +theme(plot.title = element_text(hjust = 0.5))

plot_grid(f3d, labels="C")

