### generate Multi-panel figures for the Manuscript 
library(cowplot)
library(ggplot2)
library(RColorBrewer)
library(ggsignif)
library(ggbeeswarm)
library(ggpubr)
library(plyr)

## Figure 1ab
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Figure2_cowplot.R")
data<-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/alpha_diversity_Regional.csv", header=TRUE, row.names=1)

generate_adiv_plots <- function(input_data, X, Y, min,max){
  data<-as.data.frame(input_data)
  #Ensure correct ordering of levels 
  data$Site_General <- factor(data$Site_General, levels = c("SI", "Colon"))
  data$Site <- factor(data$Site, levels = c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  ggplot(data=data,aes(x={{X}},y={{Y}}, fill={{X}})) + 
    geom_violin(alpha=0.25,position=position_dodge(width=.75),size=1,color="black",draw_quantiles=c(0.5))+
    #geom_boxplot(alpha=0.25)+ 
    #geom_quasirandom(alpha=0.1)+
    scale_fill_viridis_d()+
    geom_point(size=1,position=position_jitter(width=0.25),alpha=1)+
    theme_cowplot(12) +
    ylim(min,max)+
    theme(legend.position = "none")
  
}

#merge metadata with alpha diversity
data$SampleID <- gsub("-",".",data$X)
metadata<- read.table("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/alpha_Regional-ASV-table_d11238/Regional-Combat-Metadata.tsv", header=TRUE)
intermediate<- (merge(data, metadata, by = 'SampleID'))
data<- intermediate

names(data)
data$Sequencing_Run= factor(data$Sequencing_Run)
data$Type= revalue(data$Type, c("Luminal" = "Lum", "Mucosal" ="Muc"))
data$Site_General = factor(data$Site_General, levels=c("Colon", "SI"))
data$Site = factor(data$Site, levels= c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec", "Ileum"="Ile", "Jejunum" = "Jej", "Duodenum"= "Duo"))
sapply(data,levels)

luminaldata<-filter(data, Type=="Luminal")
mucosaldata<-filter(data, Type =="Mucosal")
colondata<-filter(data, Site_General =="Colon")
SIdata<-filter(data, Site_General =="SI")
 
#Aggregate Plot by Type Colon Data by Site
otus_col <- generate_adiv_plots(colondata, Type, observed_otus, 0, 600)+ facet_grid(~Site) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_col <- generate_adiv_plots(colondata, Type, pielou_e, 0, 1)+ facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

#Aggregate Plot by Type SI Data by site
otus_si <- generate_adiv_plots(SIdata, Type, observed_otus, 0, 600) + facet_grid(~Site) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)
pielou_e_si <- generate_adiv_plots(SIdata, Type, pielou_e, 0, 1) + facet_grid(~Site) +
  stat_compare_means(comparisons = list(c("Lum", "Muc")),method="wilcox",vjust=0.5,label="p.signif",step.increase=0.08, hide.ns = TRUE)

## Figure 1de
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Duodenum PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Jejunum PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Ileum PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Cecum PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Proximal_Colon PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Distal_Colon PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Colon PcoA.csv", header = TRUE)
data <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/SI PcoA.csv", header = TRUE)

metadata <- read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/RPCA_LuminalvsMucosal_PcoA/Regional-Metadata-All.csv", header=TRUE)
intermediate<- (merge(data, metadata, by = 'Description'))
data<- intermediate
data$Type= revalue(data$Type, c("Luminal" = "Lum", "Mucosal" ="Muc"))

#Colon

cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))
colon <- p + labs(title="Colon")+theme(plot.title = element_text(hjust = 0.5))


#SI 

cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))
SI <- p + labs(title="SI") +theme(plot.title = element_text(hjust = 0.5))

data$Site = revalue(data$Site, c("Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))
duo <- p + facet_grid(Site~.) + ggtitle ("SI") +theme(plot.title = element_text(hjust = 0.5)) 

data$Site = revalue(data$Site, c("Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="none") 
jej <- p + facet_grid(Site~.) 

data$Site = revalue(data$Site, c("Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="none") 
ile <- p + facet_grid(Site~.) 

dev.new(width=15, height=10)
fig2c<-plot_grid (duo,jej,ile, ncol=1, align ="hv")

dev.new(width=15, height=10)
plot_grid(test,test2,test3, ncol=1, align ="hv")

dev.new(width=5, height=5)
test<-duo + coord_fixed(1/1)
test2<-jej+coord_fixed(1/1)
test3<- ile+coord_fixed(1/1)
#Colon
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec"))

cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="top",legend.justification = "center") +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid"))
cec <- p + facet_grid(Site~.) + ggtitle ("Colon") +theme(plot.title = element_text(hjust = 0.5)) 

data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec"))
cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="none") 
PC <- p + facet_grid(Site~.) 

data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon"="PC", "Cecum"= "Cec"))
cols<-c("Lum"="#481567FF", "Muc" = "#3CBB75FF")
p <- ggplot(data, aes(x=PC1, y=PC2, fill=Type)) + 
  geom_point(aes(fill=Type), colour="black", pch=21, size=3) + xlab("PC1") + ylab("PC2") +
  scale_fill_manual(name="Type", values=cols) +
  theme_cowplot(12) +
  #coord_fixed(ratio=1/1)+
  theme(legend.position="none") 
DC <- p + facet_grid(Site~.)



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

### Taxa Barplots ---
assign_cols<-readRDS("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/assign_cols.RDS")
L6_col <-generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Col_LumMuc_L6.RDS","Colon Genera", ".*g__",assign_cols, "SiteTypeColon") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L6_SI <- generate_L6_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/SI_LumMuc_L6.RDS", "SI Genera", ".*g__",assign_cols, "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig2h <- plot_grid(L6_SI,L6_col, NULL, ncol=3,labels=c("H",""), rel_widths = c(1.5,1.5,1))

L2_col <-generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/Col_LumMuc_L2.csv","Colon Phyla", ".*p__", "SiteTypeColon") +
  theme(legend.position = "none")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
L2_SI <- generate_L2_taxa_plots("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Taxa-Barplots/SI_LumMuc_L2.csv", "SI Phyla", ".*p__", "SiteTypeSI") +
  theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
fig2g <- plot_grid(L2_SI,L2_col, NULL, ncol=3,labels=c("G",""), rel_widths = c(2,2,1))

### Generate multi-panel figure ---
fig2b <- plot_grid(otus_col, pielou_e_col, nrow=2)
fig2a <- plot_grid(otus, pielou_e, nrow=2) 

fig2ab <-plot_grid(fig2a, fig2b,ncol=2, labels=c("A", "B"))
dev.new(width=15, height=10)
fig2ab

dev.new(width=15, height=10)
fig2d <- plot_grid (cec,PC,DC, ncol=1, align ="hv")

dev.new(width=15, height=10)
fig2cd<-plot_grid(fig2c,fig2d,ncol=2, labels=c("C", "D"))

dev.new(width=15, height=10)
plot_grid(fig2ab,fig2cd,ncol=2, align="hv")

#Cross-Sectional
fig2abcd<-plot_grid(otus_si, otus_col, pielou_e_si,pielou_e_col, SI, colon, nrow=3, labels=c("A", "B", "","", "C","D"))

fig2e <- plot_grid (cec,PC,DC, ncol=1, align ="hv")
fig2f<-plot_grid(duo,jej,ile, ncol=1, align ="hv")
fig2ef <- plot_grid(fig2f,fig2e,ncol=2, align="hv", labels=c("E", "F"))

left <- plot_grid(fig2abcd,fig2g, nrow=2,rel_heights = c(2,1))
right <- plot_grid(fig2ef,fig2h, nrow=2,rel_heights = c(2,1))

dev.new(width=15, height=10)
plot_grid(left,right)

