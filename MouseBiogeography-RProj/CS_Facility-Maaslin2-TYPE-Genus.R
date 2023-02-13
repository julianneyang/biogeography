library(Maaslin2)
library(funrar)
library(dplyr)

input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Distal_Colon_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Proximal_Colon_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Cecum_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Ileum_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Jejunum_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Duodenum_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - Colon_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file
input_data <- read.csv("CS-Facility-Analysis/Type_L6/Maaslin2_L6 SITE and TYPE - SI_L6.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
sapply(df_input_metadata,levels)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

sapply(df_input_metadata,levels)

#Distal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE, plot_scatter = FALSE)

#Proximal_Colon 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE )

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID", fixed_effects = c("Sequencing_Run","Sex", "Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), random_effects = c("MouseID"),normalization="clr", transform ="none", plot_heatmap = FALSE,plot_scatter = FALSE)

## Make a heatmap ---
duodenum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
a<-duodenum_significant$feature
jejunum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
b<-jejunum_significant$feature
ileum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
ileum_significant<-filter(ileum, metadata=="Type" & value=="Mucosal" &qval<0.05)
c<-ileum_significant$feature
cecum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
d<-cecum_significant$feature  
pc<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
e<-pc_significant$feature  
DC<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

#Query the target vector against all_results.tsv for each site
duodenum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
duodenum_all$Site<- "Duodenum"
jejunum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
jejunum_all$Site<- "Jejunum"
ileum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
ileum_all<-ileum_all[match(target,ileum_all$feature),]
ileum_all$Site<- "Ileum"
cecum<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
cecum_all<-cecum_all[match(target,cecum_all$feature),]
cecum_all$Site<- "Cecum"
pc<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
pc_all<-pc_all[match(target,pc_all$feature),]
pc_all$Site<- "Proximal_Colon"
DC<-read.table("CS-Facility-Analysis/Type_L6/L6-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
DC_all<-DC_all[match(target,DC_all$feature),]
DC_all$Site<- "Distal_Colon"

duojej<-rbind(duodenum_all,jejunum_all)
ilecec<-rbind(ileum_all, cecum_all)
pcdc<-rbind(pc_all,DC_all)
duojejilecec<-rbind(duojej,ilecec)
duojejilececpcdc<-rbind(duojejilecec,pcdc)

#write.csv(duojejilececpcdc, "CS-Facility-Analysis/Type_L6/Genus-TYPE-Heatmap.csv") 
#from here make sure all NA rows are filled with feature name corresponding to NA via copy paste
#remove across six sites all GMM that failed to converge
duojejilececpcdc<-read.csv("CS-Facility-Analysis/Type_L6/Genus-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
write.csv(offtarget, "CS-Facility-Analysis/Type_L6/deleted_taxa.csv")
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

#for dropout taxa viz 
gmm_heatmap_final<-subset(gmm_heatmap,  gmm_heatmap[,3] %in% offtarget )
#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("CS-Facility-Analysis/Type_L6/genus_taxonomy.csv", header=TRUE)
annotation$feature<-gsub("; ",".",annotation$feature)
annotation$feature<-gsub(".s__.*","",annotation$feature)
annotation$feature<-gsub("/",".",annotation$feature)
data<- (merge(gmm_heatmap_final, annotation, by = 'feature'))
data$Family_Genus<-paste(data$Family,data$Genus,sep=" : ")
data$Phylum_Genus<-paste(data$Phylum,data$Genus,sep=" : ")

qval<-data$qval
qval <- replace_na(qval, value =100)
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
data$Site = plyr::revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
data$Site = factor(data$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
ggplotdata<-data
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(ggplotdata, aes(x = Site, y=Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  cowplot::theme_cowplot(12) +
  #theme(legend.position="top",legend.justification = "center") +
  theme(legend.position="right") +
  xlab("")+
  ylab("") +
  #guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))
  guides(fill=guide_colourbar(title="",label=TRUE,barheight = 15))

  #ggtitle("Mucosal vs Luminal - CS Facility")
