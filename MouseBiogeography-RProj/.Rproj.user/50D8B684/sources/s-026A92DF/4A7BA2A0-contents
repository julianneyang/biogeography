library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/Hum-Maaslin2-SITE-Genus-SPF.R")
input_data <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/Luminal_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/Mucosal_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
sapply(df_input_metadata,levels)

#Luminal Site
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Site
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Luminal Site_General
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Site_General
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Heatmap-----

luminal<-read.table("Humanized-Biogeography-Analysis/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv", header=TRUE)
luminal<-read.table("Humanized-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv", header=TRUE)

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

immdefsignificant<-target
immdefsignificant<-as.data.frame(immdefsignificant)
immdefsignificant$feature <- immdefsignificant[,1]

annotation <- read.csv("Humanized-Biogeography-Analysis/Genus_Taxonomy.csv", header=TRUE)
annotation$feature<-annotation$taxonomy
annotation$feature<-gsub("; ",".",annotation$feature)
annotation$feature<-gsub(".s__.*","",annotation$feature)
annotation$feature<-gsub("-",".",annotation$feature)
annotation$feature<-gsub("/",".",annotation$feature)

tempdf<- (merge(immdefsignificant, annotation, by = 'feature'))
genera<-unique(tempdf$Genus)
readr::write_rds(genera, "Hum_SPF_mucosal_genera.RDS")

luminal<-read.table("Humanized-Biogeography-Analysis/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv", header=TRUE)
luminal<-read.table("Humanized-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv", header=TRUE)
luminal<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv", header=TRUE)
data <- luminal %>% filter(metadata=="Site_General" & qval<0.05)
luminal<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/Maaslin2_Site_L6/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID/all_results.tsv", header=TRUE)
data <- luminal %>% filter(metadata=="Site_General" & qval<0.05)


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
annotation <- read.csv("Humanized-Biogeography-Analysis/Genus_Taxonomy.csv", header=TRUE)
annotation$feature<-annotation$taxonomy
annotation$feature<-gsub("; ",".",annotation$feature)
annotation$feature<-gsub(".s__.*","",annotation$feature)
annotation$feature<-gsub("-",".",annotation$feature)
annotation$feature<-gsub("/",".",annotation$feature)

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

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(ggplotdata, aes(x = value, y=Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="top",legend.justification = "center") +
  xlab("")+
  ylab("") +
  guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))+
  ggtitle("Mucosal - SPF Mice")

##Make Venn diagram
library(ggvenn)
csgenera <- readRDS("CS_Facility_Mucosal_Genera.RDS")
genera <- readRDS("Hum_SPF_mucosal_genera.RDS")

sanitycheck <- c(setdiff(genera, regionalgenera), setdiff(regionalgenera, genera))
sanitycheck <- sort(sanitycheck)
sanitycheck

A=list()
A=list("CS Mucosal" = csgenera, "CS Gavage, Mucosal" = genera)
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggvenn(A,
       fill_color=c("red", "purple"),
       text_size = 8,
       set_name_size = 8)

?ggvenn()
