BiocManager::install(c("sva"))
library("sva")
library("dplyr")
library("plyr")
##NOTE: Metadata file needs to exactly match the data file in the order of samples
setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/ImmDef-Mouse-Biogeography-Analysis/")
newmeta=read.csv("ImmDef-Overlap-Combat-ASV - ImmDef-Combat-Metadata.csv", header=T,row.names=1)

newmeta=as.data.frame(newmeta)
filtereddata=read.csv("ImmDef-Overlap-Combat-ASV - ImmDef-Overlap-ASV.csv", row.names=1)
#filtereddata= filtereddata[,-c(90)] #omit taxonomy column 
warnings()
batch=newmeta$Sequencing_Run
modcombat=model.matrix(~ Sex + Site,data=newmeta)
#modcombat=model.matrix(~ Group + Obesity,data=newmeta)
#group=newmeta$Group


input=as.matrix(filtereddata)
combat_adjusted_counts=ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  
#combat_adjusted_counts=ComBat_seq(input,batch=batch,group=group)  
#combat_adjusted_counts=ComBat_seq(input,batch=batch)  
#write.table(combat_adjusted_counts,"...  .txt",sep='\t',col.names=NA,row.names=TRUE) 
write.csv(combat_adjusted_counts,"ImmDef-Overlap-ComBat-Adjusted-ASV.csv") 

?ComBat_seq
