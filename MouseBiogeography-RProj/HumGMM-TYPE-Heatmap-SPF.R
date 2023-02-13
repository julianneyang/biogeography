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
library(seecolor)
library(dplyr)
library(plyr)

here::i_am("MouseBiogeography-RProj/HumGMM-TYPE-Heatmap-SPF.R") #conflicts with plyr so use the :: notation
here::here()

###for TYPE:Mucosal vs Luminal Data
#Feed in the significant results and generate a target vector with the union of all features 
duodenum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
duodenum_significant<-filter(duodenum, metadata=="Type" & value=="Mucosal" &qval<0.05)
a<-duodenum_significant$feature
jejunum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
jejunum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
b<-jejunum_significant$feature
ileum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv", header=TRUE)
ileum_significant<-filter(jejunum, metadata=="Type" & value=="Mucosal" &qval<0.05)
c<-ileum_significant$feature
cecum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",header=TRUE)
cecum_significant<-filter(cecum, metadata=="Type" & value=="Mucosal" &qval<0.05)
d<-cecum_significant$feature  
pc<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",header=TRUE)
pc_significant<-filter(pc, metadata=="Type" & value=="Mucosal" &qval<0.05)
e<-pc_significant$feature  
DC<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/significant_results.tsv",header=TRUE)
DC_significant<-filter(DC, metadata=="Type" & value=="Mucosal" &qval<0.05)
f<-DC_significant$feature  
joinab<- union(a,b)
joincd<- union(c,d)
joinef<- union(e,f)
joinabcd <- union(joinab,joincd)
target<-union(joinabcd,joinef)

#Query the target vector against all_results.tsv for each site
duodenum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
duodenum_all<-filter(duodenum, metadata=="Type" & value=="Mucosal")
duodenum_all<-duodenum_all[match(target,duodenum_all$feature),]
duodenum_all$Site<- "Duodenum"
jejunum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
jejunum_all<-filter(jejunum, metadata=="Type" & value=="Mucosal")
jejunum_all<-jejunum_all[match(target,jejunum_all$feature),]
jejunum_all$Site<- "Jejunum"
ileum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
ileum_all<-filter(ileum, metadata=="Type" & value=="Mucosal")
ileum_all<-ileum_all[match(target,ileum_all$feature),]
ileum_all$Site<- "Ileum"
cecum<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
cecum_all<-filter(cecum, metadata=="Type" & value=="Mucosal")
cecum_all<-cecum_all[match(target,cecum_all$feature),]
cecum_all$Site<- "Cecum"
pc<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
pc_all<-filter(pc, metadata=="Type" & value=="Mucosal")
pc_all<-pc_all[match(target,pc_all$feature),]
pc_all$Site<- "Proximal_Colon"
DC<-read.table("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID/all_results.tsv", header=TRUE)
DC_all<-filter(DC, metadata=="Type" & value=="Mucosal")
DC_all<-DC_all[match(target,DC_all$feature),]
DC_all$Site<- "Distal_Colon"

duojej<-rbind(duodenum_all,jejunum_all)
ilecec<-rbind(ileum_all, cecum_all)
pcdc<-rbind(pc_all,DC_all)
duojejilecec<-rbind(duojej,ilecec)
duojejilececpcdc<-rbind(duojejilecec,pcdc)

#write.csv(duojejilececpcdc, "Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-TYPE-Heatmap.csv") 
#from here make sure all NA rows are filled with feature name corresponding to NA

#remove across six sites all GMM that didn't exist across all datsets
duojejilececpcdc<-read.csv("Humanized-Biogeography-Analysis/Source RPCA/SPF/OMIXER-RPM/GMM-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

readr::write_rds(offtarget)

#construct the heatmap using ggplot
remove.packages("Microbiome.Biogeography")
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Microbiome.Biogeography/")
devtools::document()
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/")
devtools::install("Microbiome.Biogeography")
library("Microbiome.Biogeography")

library(viridis)
cols=c("#440154FF","#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF", "#4AC16DFF", "#9FDA3AFF", "#FDE725FF")
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)
?generate_GMM_heat_map_by_Type
mucvlum <- Microbiome.Biogeography::generate_GMM_heat_map_by_Type(gmm_heatmap_final, 
                                                                  "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                                                  metabolic_map, 
                                                                  "metabolic_map", 
                                                                  "Mucosal vs Luminal", 
                                                                  colorvector=cols, 
                                                                  breakvector=bk)
dev.new(width=15, height=10)
mucvlum

cols=c("#46337EFF", "#365C8DFF" ,"#277F8EFF", "#1FA187FF")
bk =c(-1.5, -1, -0.5, 0, 0.5)
mucvlum <- Microbiome.Biogeography::generate_GMM_heat_map_by_Type(gmm_heatmap_final, 
                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv", 
                                         Map, 
                                         "Map", 
                                         "Mucosal vs Luminal",
                                         cols,
                                         bk)
dev.new(width=15, height=10)
mucvlum
#construct the heatmap using ggplot
generate_GMM_heat_map_by_Type<- function(dataframe, path_to_Module_Key, Y, ystring, titlestring,colorvector,breakvector){
  data<-merge(gmm_heatmap_final,annotation,by="feature")
  
  
  annotation <- read.csv(path_to_Module_Key, header=TRUE)
  site_heatmap<-as.data.frame(dataframe)
  
  data<- (merge(site_heatmap, annotation, by = 'feature'))
  
  
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
  
  data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
  data <- data %>% mutate(coef_d= ifelse(coef>=1.5, 2, coef))
  data$coef_d[data$coef_d <= (-1.5)] <- (-2)
  data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))
  data$Site = factor(data$Site, levels=c("Duo", "Jej", "Ile", "Cec", "PC", "DC"))
  
  #Add ammonia and carbon dioxide and remove lactaldehyde degradation
  data <- data[!grepl("#N/A",data$Map), ]
  
  if("ammonia" %in% colnames(data) | "carbon dioxide" %in% colnames(data)){
    if("ammonia" %in% colnames(data)) {
      ammonia<- data %>% filter(Map2_ammonia=="ammonia") 
      ammonia$Map <- "ammonia"
      data<-rbind(data, ammonia)
    }
    if("carbon dioxide" %in% colnames(data)){
      co2 <- data %>% filter(Map3_carbon_dioxide=="carbon dioxide")
      co2$Map <- "carbon dioxide"
      data <- rbind(data, co2)
    }
  }
  
  #Generate Y values (what you want to show on Y axis)
  data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
  data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
  data$metabolic_map<-paste(data$Map,data$annotation,sep=" : ")
  
  #Generate median coef to show if aggregating rows by GMM
  data$Site_Map <- paste(data$Site,"_",data$Map)
  test<- data %>%                                    # Calculate mean by group
    group_by(Site_Map) %>%
    mutate(coef_mean = median(coef)) %>% 
    as.data.frame()
  
  #Set colors and breaks
  
  cols=c(colorvector)
  bk =c(breakvector)
  
  if(ystring=="hierachy_annotations"){
    y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))
    ggplotdata<-data
    
    
    g1 <- ggplot(ggplotdata, aes(x = Site, y={{Y}})) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      theme_cowplot(12) +
      theme(legend.position="right",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barheight = 5)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  else if(ystring=="metabolic_map"){
    y = tapply(data$coef_d, data$metabolic_map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$metabolic_map= factor(as.character(data$metabolic_map), levels = names(y))
    ggplotdata<-data
    
    g1 <- ggplot(ggplotdata, aes(x = Site, y={{Y}})) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      theme_cowplot(12) +
      theme(legend.position="right",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barheight = 5)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  else if(ystring=="feature_annotations"){
    y = tapply(data$coef_d, data$feature_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    data$feature_annotations= factor(as.character(data$feature_annotations), levels = names(y))
    ggplotdata<-data
    
    g1 <- ggplot(ggplotdata, aes(x = Site, y={{Y}})) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      theme_cowplot(12) +
      theme(legend.position="right",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barheight = 5)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
  }
  else if(ystring=="Map"){
    
    ggplotdata<-data
    
    
    #ggplotdata <- ggplotdata %>% mutate(coef_d_mean= ifelse(coef_mean>=1.5, 2, coef_mean))
    #ggplotdata$coef_d_mean[ggplotdata$coef_d_mean <= (-1.5)] <- (-2)
    
    y = tapply(ggplotdata$coef_d, data$Map, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
    y = sort(y, FALSE)   #switch to TRUE to reverse direction
    ggplotdata$Map= factor(as.character(data$Map), levels = names(y))
    
    print(summary(ggplotdata$coef_d))
    
    g1 <- ggplot(ggplotdata, aes(x = Site, y={{Y}})) + 
      geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
      geom_text(aes(label=asterisk)) +
      #scale_fill_stepsn(colors = cols) +
      scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
      cowplot::theme_cowplot(12) +
      theme(legend.position="right",legend.justification = "center") +
      xlab("")+
      ylab("") +
      guides(fill=guide_colourbar(title="",label=TRUE,barheight = 5)) +
      ggtitle({{titlestring}}) +
      theme(plot.title = element_text(hjust = 0.5))
    
    
  }
  
  return(g1)
}
