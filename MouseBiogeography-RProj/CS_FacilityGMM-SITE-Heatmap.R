###Purpose: Aggregate all significant results from each of 6 intestinal sites into one vector; then query this vector against "all results" output from each of six sites 
library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)
library(plyr)
library(Microbiome.Biogeography)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
here::i_am("MouseBiogeography-RProj/CS_FacilityGMM-SITE-Heatmap.R")
###for SITE:DC vs all data
find_concordant_features_across_sites <- function(filepath_to_significant_results_tsv) {
  luminal<-read.table(filepath_to_significant_results_tsv, header=TRUE)
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
  return(target)
}

find_discordant_features_across_sites <- function(filepath_to_significant_results_tsv) {
  luminal<-read.table(filepath_to_significant_results_tsv, header=TRUE)
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
  disab<- c(setdiff(a, b), setdiff(b, a))
  discd<- c(setdiff(c, d), setdiff(d, c))
  disef<- c(setdiff(e, f), setdiff(f, e))
  disabcd <- c(setdiff(disab, discd), setdiff(discd, disab))
  target<-c(setdiff(disef, disabcd), setdiff(disabcd, disef))
  return(target)
}

# Luminal - 36 concordant features 
lum_target <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

# Mucosal - 37 concordant features 
muctarget <- find_concordant_features_across_sites("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv")

### Query the target vector against all_results.tsv ---

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

muc_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                              targetvector = muctarget, 
                              path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                              Y=metabolic_map,
                              ystring= "metabolic_map",
                              "Mucosal",
                              cols,
                              bk)

dev.new(width=15, height=10)
muc_GMM_map

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
muc_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Muc-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = muctarget, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Mucosal", 
                                             cols,
                                             bk)
dev.new(width=15, height=10)
muc_GMM_map

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
lum_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=metabolic_map,
                                             ystring= "metabolic_map",
                                             "Luminal",
                                             cols,
                                             bk)

dev.new(width=15, height=10)
lum_GMM_map

cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5)
lum_GMM_map <- generate_GMM_heat_map_by_site("CS-Facility-Analysis/OMIXER-RPM Results/CS_GMM/GMM-DCvsAll-CLR-Lum-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                             targetvector = lum_target, 
                                             path_to_Module_Key = "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                             Y=Map,
                                             ystring= "Map",
                                             "Luminal",
                                             cols,
                                             bk)

dev.new(width=15, height=10)
lum_GMM_map

# Y- Choose from feature_annotations, hierachy_annotations, metabolic_map

generate_GMM_heat_map_by_site <- function(path_to_all_results_tsv, targetvector, path_to_Module_Key, Y, ystring,titlestring){
  
    luminal<-read.table(path_to_all_results_tsv, header=TRUE)
    target <- targetvector
    luminal_all<-filter(luminal, metadata=="Site")
    data<-luminal_all[luminal_all$feature %in% target, ]
    
    # make an empty dataframe to store the reference variable and assign x, a string vector, to y as its column names:
    y <- data.frame(matrix(NA,nrow=length(target),ncol=9))
    x <- c(colnames(data))
    colnames(y) <- x
    y$feature<-target
    y$coef <- 0
    y$value <- "Distal_Colon"
    y$metadata <-"Site"
    y$qval<-100
    
    site_heatmap<-rbind(data,y)
    
  
    #construct the heatmap using ggplot
    annotation <- read.csv(path_to_Module_Key, header=TRUE)
    
    data<- (merge(site_heatmap, annotation, by = 'feature'))
    data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
    data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
    data$metabolic_map<-paste(data$Map,data$annotation,sep=" : ")

    #write.csv(data,"Luminal-DCvsall-ggplot-Heatmap.csv")
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
    data$value<-factor(data$value, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
    data <- data %>% mutate(coef_d= ifelse(coef>2, 2, coef))
    data$coef_d[data$coef_d < (-2)] <- (-2)
    data$value = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
    data$value = factor(data$value, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
    
  if(ystring=="hierachy_annotations"){
    y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
  }
  else if(ystring=="metabolic_map"){
    y = tapply(data$coef_d, data$metabolic_map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$metabolic_map= factor(as.character(data$metabolic_map), levels = names(y))
  }
  else if(ystring=="feature_annotations"){
    y = tapply(data$coef_d, data$feature_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$feature_annotations= factor(as.character(data$feature_annotations), levels = names(y))
  }
   
    ggplotdata<-data
    cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
    bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
    
    g1 <- ggplot(ggplotdata, aes(x = value, y={{Y}})) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      theme_cowplot(12) +
      theme(legend.position="top",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barwidth = 15)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))

    g1
}


