### generate Multi-panel figures for the Manuscript 
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(plyr)
library(ggpubr)
library(tidyr)
library(dplyr)
library(gridExtra)

remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

### Alpha Diversity ---

here::i_am("MouseBiogeography-RProj/Figure1_cowplot.R")

data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv", header=TRUE, row.names=1)

generate_adiv_plots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max)+
    theme(legend.position = "none")
  
}

#merge metadata with alpha diversity
data$SampleID <- gsub("-",".",data$X)
write.csv(data,"Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv")
metadata<- read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.tsv", header=TRUE)
write.csv(metadata,"Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.csv")
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= factor(data$Type, levels=c("Luminal", "Mucosal"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
sapply(data,levels)

luminaldata<-data %>% dplyr::filter(Type =="Luminal")
mucosaldata<-data %>% dplyr::filter(Type =="Mucosal")
colondata<-data %>% dplyr::filter(Site_General =="Colon")
SIdata<-data %>% dplyr::filter(Site_General =="SI")

#Luminal Data
otus<- generate_adiv_plots(luminaldata, Site, observed_otus, 0, 725) + ggtitle ("Luminal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e<- generate_adiv_plots(luminaldata, Site, pielou_e,0,1.1) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
#Mucosal Data

otus_m <- generate_adiv_plots(mucosaldata, Site, observed_otus, 0, 725) + ggtitle ("Mucosal") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox", vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_m <- generate_adiv_plots(mucosaldata, Site, pielou_e,0,1.1) +
  stat_compare_means(comparisons = list(c("DC", "PC"),
                                        c("DC", "Cec"),
                                        c("DC", "Ile"),
                                        c("DC", "Jej"),
                                        c("DC", "Duo")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

### Beta Diversity ---
#Luminal
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_Luminal_PcoA.csv", row.names=1)
data$Site <- factor(data$Site.1)
data$Site_General <- factor(data$Site_General)
cols <- c("SI" = "#F8766D","Colon" ="#00BFC4")
luminal <- ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
  geom_point(aes(fill=Site_General), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="",values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
luminal <- luminal + labs(title="Luminal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

#Mucosal
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_Mucosal_PcoA.csv", row.names=1)
data$Site <- factor(data$Site.1)
data$Site_General <- factor(data$Site_General)
cols <- c("SI" = "#F8766D","Colon" ="#00BFC4")
mucosal <- ggplot(data, aes(x=PC1, y=PC2, colour=Site_General)) + 
  geom_point(aes(fill=Site_General), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="",values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
mucosal <- mucosal + labs(title="Mucosal") + 
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  theme(plot.title = element_text(hjust = 0.5))

#Luminal Colon
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_LumCol_PcoA.csv", row.names=1)
data$Site = factor(data$Site.1, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site.1, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec"))
cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
#dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Site)) + 
  geom_point(aes(fill=Site), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
lumcol <- p + labs(title="Colon- Luminal") + theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +theme(plot.title = element_text(hjust = 0.5))

#Luminal SI
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_LumSI_PcoA.csv", row.names=1)
data$Site = factor(data$Site.1, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site.1, c("Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
data$Site = factor(data$Site, levels=c("Duo", "Jej", "Ile"))
cols <- c("Duo" = "firebrick", "Jej"="gold", "Ile" = "forestgreen")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Site)) + 
  geom_point(aes(fill=Site), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
silum <- p + labs(title="SI- Luminal") + theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +theme(plot.title = element_text(hjust = 0.5))

#Mucosal Colon
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_MucCol_PcoA.csv", row.names=1)
data$Site = factor(data$Site.1, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec"))
cols <- c("Cec" = "cyan", "PC" = "blue", "DC" = "magenta")
#dev.new(width=12, height=10)
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Site)) + 
  geom_point(aes(fill=Site), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
muccol <- p + labs(title="Colon- Mucosal") + theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +theme(plot.title = element_text(hjust = 0.5))

#Mucosal SI 
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA-PCoA/RPCA for all Sites - RPCA_MucSI_PcoA.csv", row.names=1)
data$Site = factor(data$Site.1, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site.1, c("Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
data$Site = factor(data$Site, levels=c("Duo", "Jej", "Ile"))
cols <- c("Duo" = "firebrick", "Jej"="gold", "Ile" = "forestgreen")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Site)) + 
  geom_point(aes(fill=Site), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") 
mucsi <- p + labs(title="SI- Mucosal") + theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +theme(plot.title = element_text(hjust = 0.5))

### Taxa Barplots ---

assign_cols <- readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/assign_cols.RDS")
??Microbiome.Biogeography
L2_lum <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Luminal_L2.csv", "Luminal Phyla", ".*p__", "Site") +
  theme(legend.position = "none")
L2_muc <-generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L2.csv", "Mucosal Phyla", ".*p__", "Site") +
  theme(legend.position = "none")

L6_lum <- Microbiome.Biogeography::generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Luminal_L6.RDS", "Luminal Genera", ".*g__", assign_cols, "Site") +
  theme(legend.position = "none")
L6_muc <- Microbiome.Biogeography::generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L6.RDS","Mucosal Genera", ".*g__",assign_cols, "Site") +
  theme(legend.position = "none")

# Draw legend
L6_legend <- Microbiome.Biogeography::generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L6.RDS","Mucosal ( > 0.1% Relative Abundance)", ".*g__", assign_cols, "Site") +
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L6_legend)
grid::grid.newpage()
grid::grid.draw(legend)
L2_legend <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Mucosal_L2.csv","Mucosal Phyla", ".*p__", "Site")+
  theme(legend.position = "right") +
  guides(fill=guide_legend(nrow=22, byrow=TRUE))+
  theme(legend.spacing.y = unit(0.1, 'cm')) +
  theme(legend.background = element_rect(fill="lightblue", size=1, linetype="solid"), legend.margin = margin(0, 11, 0, 1)) 
legend <- cowplot::get_legend(L2_legend)
grid::grid.newpage()
grid::grid.draw(legend)


###Heatmaps: Figs 1fg
#Feed in the significant results and generate a target vector with the union of all features 
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/")
luminal<-read.table("L6-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv", header=TRUE)
luminal<-read.table("L6-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv", header=TRUE)

duodenum_significant<-filter(luminal, metadata=="Site" & value=="Duodenum" &qval<0.05)
a<-duodenum_significant$feature
jejunum_significant<-filter(luminal, metadata=="Site" & value=="Jejunum" &qval<0.05)
b<-jejunum_significant$feature
ileum_significant<-filter(luminal, metadata=="Site" & value=="Ileum" &qval<0.05)
c<-ileum_significant$feature
cecum_significant<-filter(luminal, metadata=="Site" & value=="Cecum" &qval<0.05)
d<-cecum_significant$feature  
pc_significant<-filter(luminal, metadata=="Site" & value=="Proximal_Colon" &qval<0.05)
e<-pc_significant$feature  
DC_significant<-filter(luminal, metadata=="Site" & value=="Distal_Colon" &qval<0.05)
f<-DC_significant$feature  

joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

#Query the target vector against all_results.tsv 
luminal<-read.table("L6-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)
luminal<-read.table("L6-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv", header=TRUE)

luminal_all<-filter(luminal, metadata=="Site")
#length(luminal_all$value[luminal_all$value=="Duodenum"])
data<-luminal_all[luminal_all$feature %in% target, ]

length(target)
#make an empty dataframe to store the reference variable 
y <- data.frame(matrix(NA,nrow=length(target),ncol=9))
#Assign x, a string vector, to y as its column names:
x <- c(colnames(data))
colnames(y) <- x
y$feature<-target
y$coef <- 0
y$value <- "Distal_Colon"
y$metadata <-"Site"
y$qval<-100

site_heatmap<-rbind(data,y)

site_heatmap$feature <- gsub("X","",as.character(site_heatmap$feature))
#write.csv(site_heatmap,"SITE Genus Heatmap.csv")

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("genus_Luminal_taxonomy.csv", header=TRUE)
annotation <- read.csv("genus_Mucosal_taxonomy.csv", header=TRUE)
annotation$feature<-annotation$X
data<- (merge(site_heatmap, annotation, by = 'feature'))
data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")

qval<-data$qval
asterisk<-c("")
for (item in qval){
  if (item < 0.05){
    asterisk<-c(asterisk,"*")
  }
  else if (item=="NA"){}
  else {
    asterisk<-c(asterisk,"")
  }
}
asterisk<-asterisk[-1]
data$asterisk<-asterisk
data$value<-factor(data$value, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
data$coef_d[data$coef_d < (-2)] <- (-2)
summary(data$coef_d) 
  y = tapply(data$coef_d, data$Genus, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
  y = sort(y, FALSE)   #switch to TRUE to reverse direction
  data$Genus= factor(as.character(data$Genus), levels = names(y))
  data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
ggplotdata<-data
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

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
write.csv(ggplotdata, ("Phyla_Colors_Mucosal.csv"))

g1 <- ggplot(ggplotdata, aes(x = value, y=Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="top",legend.justification = "center") +
  xlab("")+
  ylab("") +
  guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))
dev.new(width=15, height=10) 
muc <- g1 +ggtitle("Mucosal") +theme(plot.title = element_text(hjust = 0.5))
lum <- g1 +ggtitle("Luminal") +theme(plot.title = element_text(hjust = 0.5))

fig1fg <- plot_grid(lum, muc, align= "hv", labels = c("F","G"))
dev.new(width=15, height=10)
fig1fg
fig1bcde

### Generating a Multi-Panel figure ---

dev.new(width=15, height=20)
plot_grid(otus, otus_m, pielou_e,pielou_m, cols=2, rows=2, align = 'hv',labels = c("B","C","",""))
fig1bc <- plot_grid(otus, otus_m, pielou_e,pielou_m, cols=2, rows=2, align = 'hv',labels = c("B","C","",""))

lum <-plot_grid(silum,lumcol, nrow=2, align="hv")
muc <-plot_grid(mucsi,muccol, nrow=2, align="hv")
fig1de<-plot_grid(lum, muc, ncol=2,  align = 'hv',labels = c("D","E"))
dev.new(width=12, height=10)
fig1de
dev.new(width=12, height=10)
plot_grid(fig1bc,fig1de, align="hv")

#Longitudinal
fig1bc<-plot_grid(otus, otus_m, pielou_e,pielou_m, nrow=2, labels=c("B", "C", "",""))

warnings()
lum <-plot_grid(luminal, silum,lumcol, nrow=3, ncol=1, align="hv")
muc <-plot_grid(mucosal, mucsi,muccol, nrow=3, ncol=1,align="hv")
fig1fg<-plot_grid(lum, muc, ncol=2,labels = c("D","E"))

fig1h <- plot_grid(L2_lum,L2_muc, NULL, ncol=3,labels=c("F",""), rel_widths = c(2,2,1))
fig1i <- plot_grid(L6_lum,L6_muc, NULL, ncol=3,labels=c("G",""), rel_widths = c(1.5,1.5,1))

left <- plot_grid(fig1bc, fig1h, nrow=2, rel_heights = c(2,1))
right <- plot_grid(fig1fg,fig1i, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)

fig7abcdef <- plot_grid(left,right)

