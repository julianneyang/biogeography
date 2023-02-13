library(Maaslin2)
library(funrar)
library(dplyr)
library(ggplot2)
library(cowplot)
library(plyr)
library(here)
library(lessR)
library(ggplot2)
library(tidyr)
library(gplots)
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/ImmDef-WTCohort-Maaslin2-Heatmap.R")
?generate_GMM_heat_map_by_site()

input_data <- read.csv("ImmDef-Mouse-Biogeography-Analysis/Maaslin2_Site_Differences/WT-Cohort-L6-min10k.csv", header=TRUE, row.names=1) # choose filtered non rarefied csv file

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("ImmDef-Mouse-Biogeography-Analysis/Full-Metadata.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata

#Mucosal Site
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "ImmDef-Mouse-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-SeqRunSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), random_effects = c("MouseID_Original"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal Site_General
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "ImmDef-Mouse-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-SeqRunSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), random_effects = c("MouseID_Original"),normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Heatmap-----

luminal<-read.table("ImmDef-Mouse-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv", header=TRUE)


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
annotation <- read.csv("ImmDef-Mouse-Biogeography-Analysis/genus_Mucosal_taxonomy.csv", header=TRUE)
annotation$feature<-annotation$X.OTU.ID
annotation$feature<-gsub("/",".",annotation$feature)
annotation$feature<-gsub("-",".",annotation$feature)
tempdf<- (merge(immdefsignificant, annotation, by = 'feature'))

immdefgenera<-tempdf$Genus
readr::write_rds(immdefgenera,"immdefgenera.RDS")

luminal<-read.table("ImmDef-Mouse-Biogeography-Analysis/L6-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv", header=TRUE)

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

site_heatmap$feature <- as.character(site_heatmap$feature)
testing <- site_heatmap$feature
#write.csv(site_heatmap,"SITE Genus Heatmap.csv")

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("ImmDef-Mouse-Biogeography-Analysis/genus_Mucosal_taxonomy.csv", header=TRUE)
annotation$feature<-annotation$X.OTU.ID
annotation$feature<-gsub("/",".",annotation$feature)
annotation$feature<-gsub("-",".",annotation$feature)

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
cols=viridis(10)

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
  ggtitle("Mucosal- WT Cohort")

#construct heatmap using heatmap2 with dendrogram
data_long<-pivot_wider(data, id_cols=Phylum_Genus, names_from = value, values_from =coef_d)
data_long_final<-data_long[,-1]
data_long_final<-select(data_long_final,Duo,Jej, Ile,Cec,PC,DC)
row.names(data_long_final)= data_long$Phylum_Genus
matrix.data<- as.matrix.data.frame(data_long_final)
library(RColorBrewer)
coul = colorRampPalette(brewer.pal(8, "Blues"))(5)#, the number (25) represents the number of shades of the gradient
dev.new(width=15, height=10)
heatmap.2(matrix.data, colv= NA, rowv =TRUE, dendrogram ="row", scale="row",density.info="none", trace="none", cexCol = 1, cexRow = 1, margins=c(9,20), col=coul, keysize=1.5, )
#looks like 6 clusters could work

#construct heatmap using heatmap2 with hierarchical clustering
library(dendextend) #Creating a color palette & color breaks

coul=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

distance= dist(matrix.data, method ="euclidean")  
hcluster = hclust(distance, method ="ward.D")

cols_branches <- c("purple", "firebrick", "forestgreen", "royalblue") # Set the colors of branches
#cols_branches <- c("firebrick", "royalblue", "forestgreen", "purple", "navy","black") # Set the colors of branches
#cols_branches <-c("purple", "pink", "grey", "navy")
#cols_branches <-c("cyan","purple", "pink", "grey", "navy")
dend1<-as.dendrogram(hcluster)
dend1 <- color_branches(dend1, k = 4, col = cols_branches)  


col_labels <- cols_branches[cutree(dend1, k = 4)] # sync with num clusters
col_labels <- get_leaves_branches_col(dend1)
col_labels <- col_labels[order(order.dendrogram(dend1))]

# dendrogram tuning from: https://stackoverflow.com/questions/29265536/how-to-color-the-branches-and-tick-labels-in-the-heatmap-2
data_long_qval<-pivot_wider(data, id_cols=feature, names_from = value, values_from =qval)
data_long_qval<-data_long_qval[,-1]

data_long_qval<-pivot_wider(data, id_cols=feature, names_from = value, values_from =qval)
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
dev.new(width=15, height=10)
heatmap.2(matrix.data,
          Colv= FALSE,
          breaks=bk,
          Rowv = as.dendrogram(hcluster),
          symkey=FALSE,
          dendrogram="row",
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

##Make Venn diagram
library(ggvenn)   
here()
regionalgenera <- readRDS("regionalgenera.RDS")
immdefgenera <- readRDS("immdefgenera.RDS")
c(setdiff(regionalgenera, immdefgenera), setdiff(immdefgenera, regionalgenera))
A=list()
A=list("Original" = regionalgenera, "Validation" = immdefgenera)
dev.new(width=15, height=10)  # can adjust window size of the plot output this way
ggvenn(A)
overlap <-intersect(regionalgenera, immdefgenera)
