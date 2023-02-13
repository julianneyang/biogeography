BiocManager::install(c("sva"))
library("sva")
library("dplyr")
library("plyr")
##NOTE: Metadata file needs to exactly match the data file in the order of samples
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/")
newmeta=read.csv("Sequencing Run Overlap ASV and Metadata - Regional-Overlap-Metadata.csv", header=T,row.names=1)

newmeta=as.data.frame(newmeta)
filtereddata=read.csv("Sequencing Run Overlap ASV and Metadata - Regional-Overlap-ASV.csv", row.names=1)
#filtereddata= filtereddata[,-c(90)] #omit taxonomy column 
warnings()
batch=newmeta$Sequencing_Run
modcombat=model.matrix(~ Type + Site,data=newmeta)
#group=newmeta$Group


input=as.matrix(filtereddata)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  
#combat_adjusted_counts=ComBat_seq(input,batch=batch,group=group)  
#combat_adjusted_counts=ComBat_seq(input,batch=batch)  
#write.table(combat_adjusted_counts,"...  .txt",sep='\t',col.names=NA,row.names=TRUE) 
write.csv(combat_adjusted_counts,"RegionalMicrobiome-Type-Site-Overlap-ComBat-Adjusted.csv") 

?ComBat_seq
