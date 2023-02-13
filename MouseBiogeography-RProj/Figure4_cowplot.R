
###Make Figure 4: Heatmap of features within modules------------ [ARCHIVE]

## Mucosal vs. Luminal ---
#remove across six sites all GMM that failed to converge
duojejilececpcdc<-read.csv("GMM-TYPE-Heatmap.csv")
gmm_heatmap<-duojejilececpcdc
discard_gmm<- gmm_heatmap[is.na(gmm_heatmap$metadata), ]
offtarget<- discard_gmm$feature
offtarget<-unique(offtarget)
gmm_heatmap_final<-subset(gmm_heatmap,  !gmm_heatmap[,3] %in% offtarget )

#construct the heatmap using ggplot
library(viridis)
annotation <- read.csv("Revised_Module_Key.csv", header=TRUE)
data<- (merge(gmm_heatmap_final, annotation, by = 'feature'))
data$Map <- factor(data$Map)
#data <- data %>% filter(Map =="central metabolism")
data$feature_annotations<-paste(data$feature,data$annotation,sep=" : ")
data$hierachy_annotations<-paste(data$Hierarchy_L2,data$annotation,sep=" : ")
data$Map_annotations<-paste(data$Map,data$annotation,sep=" : ")

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
data <- data %>% mutate(coef_d= ifelse(coef>1.0, 1.5, coef))
data$coef_d[data$coef_d < (-1.5)] <- (-2)
data$Site<-factor(data$Site, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
data$Site = revalue(data$Site, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))

summary(data$coef_d)

#don't run this unless you want a color gradient ----
y = tapply(data$coef_d, data$Map_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$Map_annotations= factor(as.character(data$Map_annotations), levels = names(y))
###-----

ggplotdata<-data
cols=viridis(8)
bk=c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

#bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)

g1 <- ggplot(ggplotdata, aes(x = Site, y=Map)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="right",legend.justification = "center") +
  xlab("")+
  ylab("") +
  #theme(axis.text.y = element_text(colour = tick_colors))+
  guides(fill=guide_colourbar(title="",label=TRUE, barheight = 15))
dev.new(width=15, height=10)
g1 +ggtitle("Mucosal vs. Luminal") +theme(plot.title = element_text(hjust = 0.5))

## DC vs all other Sites ---

data <- read.csv("Luminal-DCvsall-ggplot-Heatmap.csv", row.names=1)
data <- read.csv("Mucosal-DCvsall-ggplot-Heatmap.csv", row.names=1)

annotation <- read.csv("Revised_Module_Key.csv", header=TRUE)
data<- (merge(data, annotation, by = 'feature'))
data$Map <- factor(data$Map)


#data <- data %>% filter(Map =="central metabolism")
data$Map_annotations<-paste(data$Map,data$annotation.x,sep=" : ") #for Luminal

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
data <- data %>% mutate(coef_d= ifelse(coef>1.0, 1.5, coef))
data$coef_d[data$coef_d < (-1.5)] <- (-2)
summary(data$coef_d)
data$Site<-factor(data$value, levels=c("Duodenum", "Jejunum", "Ileum", "Cecum", "Proximal_Colon", "Distal_Colon"))
data$Site = revalue(data$value, c("Distal_Colon"="DC", "Proximal_Colon" = "PC", "Cecum" ="Cec","Ileum"="Ile", "Jejunum"="Jej", "Duodenum"= "Duo"))

y = tapply(data$coef, data$hierachy_annotations, function(y) mean(y))  # orders the genera by the highest fold change of any ASV in the genus; can change max(y) to mean(y) if you want to order genera by the average log2 fold change
y = sort(y, FALSE)   #switch to TRUE to reverse direction
data$hierachy_annotations= factor(as.character(data$hierachy_annotations), levels = names(y))

ggplotdata<-data
cols=viridis(8)
bk =c(-2, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 2)


g1 <- ggplot(ggplotdata, aes(x = Site, y=Map)) + geom_tile(aes(fill = coef_d),colour="white",size=0.25) +
  #stat_summary(fun = "mean")+
  geom_text(aes(label=asterisk)) +
  scale_fill_stepsn(breaks=bk, values = NULL, colors = cols) +
  theme_cowplot(12) +
  theme(legend.position="right",legend.justification = "center") +
  xlab("")+
  ylab("") +
  #theme(axis.text.y = element_text(colour = tick_colors))+
  guides(fill=guide_colourbar(title="",label=TRUE, barheight = 15))
dev.new(width=15, height=10)
g1 +ggtitle("Luminal") +theme(plot.title = element_text(hjust = 0.5))

