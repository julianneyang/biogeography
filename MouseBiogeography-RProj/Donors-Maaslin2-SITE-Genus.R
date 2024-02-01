library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(gplots)

here::i_am("MouseBiogeography-RProj/Donors-Maaslin2-SITE-Genus.R")
input_data <- readr::read_delim(here("Donors-Analysis/site_subsets/export_L6_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), delim="\t") # choose filtered non rarefied csv file
input_data <- subset(input_data, !grepl("Mitochondria|Chloroplast", OTU.ID))

input_data <- as.data.frame(input_data)
row.names(input_data)<-input_data$OTU.ID
df_input_data <- select(input_data, -c("taxonomy","OTU.ID"))


input_metadata <-readr::read_delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),delim="\t") #mapping file
input_metadata <- as.data.frame(input_metadata)
row.names(input_metadata) <- input_metadata$SampleID
input_metadata <- input_metadata %>% select(-c("SampleID"))

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Donor_ID<- factor(df_input_metadata$Donor_ID)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
#df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Mucosal", "Luminal"))
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))


#Luminal Site
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site_General,Colon'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site,Distal_Colon'))

## Mucosal --
input_data <- readr::read_delim(here("Donors-Analysis/site_subsets/export_L6_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV/feature-table.tsv"), delim="\t") # choose filtered non rarefied csv file

input_data <- as.data.frame(input_data)
row.names(input_data)<-input_data$OTU.ID
df_input_data <- select(input_data, -c("taxonomy","OTU.ID"))

input_metadata <-readr::read_delim(here("Donors-Analysis/Donors_Metadata.tsv"),delim="\t") #mapping file
input_metadata <- as.data.frame(input_metadata)
row.names(input_metadata) <- input_metadata$SampleID
input_metadata <- input_metadata %>% select(-c("SampleID"))

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Donor_ID<- factor(df_input_metadata$Donor_ID)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("SI","Colon"))
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
#df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Mucosal", "Luminal"))
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site_General,Colon'))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID","DonorID"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE,
                    min_prevalence=0.15,
                    reference=c('Sequencing_Run,Jan_2017','Site,Distal_Colon'))

#Heatmap-----

luminal<-readr::read_delim(here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"), delim="\t")
luminal<-readr::read_delim(here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/significant_results.tsv"), delim="\t")

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

df<-target
df<-as.data.frame(df)
df$feature <- df[,1]

df$Phylum <- gsub(".*\\.p__", "", df$feature)
df$Phylum <- gsub("\\.c__.*", "", df$Phylum)
df$Order <- gsub(".*\\.o__", "", df$feature)
df$Order <- gsub("\\.f__.*", "", df$Order)
df$Order <- paste0(df$Order, " (o)")
df$Family <- gsub(".*\\.f__", "", df$feature)
df$Family <- gsub("\\.g__.*", "", df$Family)
df$Family<- paste0(df$Family, " (f)")
df$Genus <- gsub(".*\\.g__", "", df$feature)

df <- df %>%
  mutate(annotation = ifelse(Genus!="", Genus,
                             ifelse(Family!=" (f)", Family, Order)))

readr::write_csv(df,here("Donors-Analysis/differential_genera_site/Genus_Luminal_taxonomy.csv"))
readr::write_rds(df, here("Donors-Analysis/differential_genera_site/Genus_Mucosal_taxonomy.RDS"))

luminal<-readr::read_delim(here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"), delim="\t")
luminal<-readr::read_delim(here("Donors-Analysis/differential_genera_site/L6-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID/all_results.tsv"), delim="\t")

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
annotation <- readr::read_csv("Donors-Analysis/differential_genera_site/Genus_Luminal_taxonomy.csv")
annotation <- readr::read_rds("Donors-Analysis/differential_genera_site/Genus_Mucosal_taxonomy.RDS")

data<- (merge(site_heatmap, annotation, by = 'feature'))
data$Phylum_Genus<-paste(data$Phylum,data$annotation,sep=" : ")
data<- data %>% filter(!annotation=="")

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
ggplot_data<-unique(data)


cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")

bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggplot(ggplot_data, aes(x = value, y=Phylum_Genus)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="top",legend.justification = "center") +
  xlab("")+
  ylab("") +
  guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15))+
  ggtitle("HUM V. Gavage Luminal")

##Make Heatmap with Heatmap2
#construct heatmap using heatmap2 with dendrogram
?pivot_wider
data_long<-pivot_wider(ggplot_data, id_cols=Phylum_Genus, names_from = value, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,Duo,Jej, Ile,Cec,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)

#construct heatmap using heatmap2 with hierarchical clustering
library(dendextend) #Creating a color palette & color breaks

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
matrix.data
distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

# Luminal Cluster Colors
cols_branches <- c("purple", "forestgreen", "royalblue", "firebrick") # Set the colors of branches
# Mucosal Cluster Colors
cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple") # Set the colors of branches


dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
nrow(matrix.data)
data_long_qval<-pivot_wider(ggplot_data, id_cols=annotation, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval <- select(data_long_qval,c("Duo","Ile", "Jej","Cec","PC","DC"))

for(i in 1:ncol(data_long_qval)){       # for-loop over columns
  v<-data_long_qval %>% pull(i)
  for(ctr in 1:length(v)){
    if (data_long_qval[ctr,i]<0.05){
      data_long_qval[ctr,i]<- 100
    }
    else {
      data_long_qval[ctr,i]<- 0
    }
  }
}


for(i in 1:ncol(data_long_qval)){ 
  v<-data_long_qval %>% pull(i)
  data_long_qval[,i]<-as.character(v)
}


for(i in 1:ncol(data_long_qval)){       # for-loop over columns
  v<-data_long_qval %>% pull(i)
  for(ctr in 1:length(v)){
    if (data_long_qval[ctr,i]=="100"){
      data_long_qval[ctr,i]<- "*"
    }
    else {
      data_long_qval[ctr,i]<- ""
    }
  }
}

asterisk_matrix<-as.matrix.data.frame(data_long_qval)


#row.names(matrix.data) <- c()
dev.new(width=15, height=10)
heatmap.2(matrix.data,
          Colv= FALSE,
          breaks=bk,
          Rowv = as.dendrogram(hcluster),
          symkey=FALSE,
          dendrogram="none",
          scale="none",
          key.xlab="coef",
          density.info="none",
          trace="none",
          cexCol = 1,
          cexRow = 1,
          margins=c(5,25),
          col=coul,
          keysize=0.5,
          RowSideColors = col_labels,
          cellnote=asterisk_matrix,
          notecol="black",
          labCol = c("Duo", "Jej","Ile","Cec","PC","DC"),
          colRow=col_labels,
          srtCol=0)
