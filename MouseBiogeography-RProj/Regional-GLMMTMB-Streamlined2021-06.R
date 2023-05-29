library("glmmTMB")
library("bbmle") ## for AICtab
library("qvalue")
library("funrar")
library(dplyr)
library(tidyverse)
library(data.table)
library(janitor)


set.seed(123)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography/SubRegional-Mouse-Biogeography-Analysis/")

input_data <- read.csv(file= "MucosalColon-DiffAbundanceTestingFiles - MucosalColon-GLMMTMB Input.csv", header=TRUE, row.names=1)

Skupskydata<-as.data.frame(input_data)

#Take out only the ASV count data and store it into SkupskyASVonly or read just the ASV file 
SkupskyASVonly = Skupskydata[, -c(275:285)] #mucosalcolon
View(SkupskyASVonly)
SkupskyASVonly <- as.matrix(SkupskyASVonly)

#Double check that your variables (all the ones you will be using in your model) have the appropriate number of levels

Skupskydata$Sex <- factor(Skupskydata$Sex)
Skupskydata$Line<- factor(Skupskydata$Line)
Skupskydata$Sequencing_Run <- factor(Skupskydata$Sequencing_Run)
Skupskydata$Site <- factor(Skupskydata$Site)
#Skupskydata$Cre <- factor(Skupskydata$Cre)
Skupskydata$MouseID_Line <- factor(Skupskydata$MouseID_Line)
sapply(Skupskydata,levels)

#Store your metadata separately 
Sex <- Skupskydata$Sex
Sequencing_Run <- Skupskydata$Sequencing_Run
Line <- Skupskydata$Line
Site <- Skupskydata$Site
N <- Skupskydata[, "Counts_by_Sample"]
MouseID_Line <- Skupskydata$MouseID_Line


out2<-data.frame() 
out3<-data.frame()
out4<-data.frame() 
out5 <- list() #store full outputs here. then you can call each individual full result 
for (j in 1:ncol(SkupskyASVonly)){
  y = as.numeric(SkupskyASVonly[, j]) 
  out1 = glmmTMB(y ~ Sequencing_Run + Sex + Site + offset(log(N)) + (1|Line/MouseID_Line), family=nbinom2) #SLC Total Pairs
  output <- summary(out1) 
  outputcoefficients <- output$coefficients 
  tTable_fixed <- outputcoefficients$cond 
  tTable_zinf <- outputcoefficients$zi 
  tTable_AIC <- output$AICtab 
  out2<-rbind(out2, tTable_fixed) 
  out3<-rbind(out3, tTable_zinf) 
  out4<-rbind(out4, tTable_AIC) 
  #out5[j] <- output, if you want this script to take forever but take an individual look at the results for a taxon which caused problems
}
warnings() 
out4 = as.data.frame(out4) 
rownames(out4) = colnames(SkupskyASVonly)
#write.csv(out4, "LuminalColon-vil-nonzero_SiteSexCre_1-Litter-MsID_GLMMTMB_AICtable.csv")
#out3 = as.data.frame(out3) 
#rownames(out3) = colnames(SkupskyASVonly)
#write.csv(out3, "nonzero_SiteSexGen_1-Study-Litter_lumen_SI_GLMMTMB_ZeroInfValues.csv")

out2$bonferroni<-p.adjust(out2$`Pr(>|z|)`, method="bonferroni")
out2$BH<-p.adjust(out2$`Pr(>|z|)`, method="BH")
out2$`Pr(>|z|)`[is.na(out2$`Pr(>|z|)`)] <- 1
out2_Qvalue=qvalue(out2$`Pr(>|z|)`,pi0.method="bootstrap",robust=FALSE)
out2$Qvalues=out2_Qvalue$qvalues

#write.csv(out2, "LuminalColon-vil-nonzero_SiteSexCre_1-Litter-MsID_GLMMTMB_FixedCoeff.csv")

###################Now, Run everything to the end to get your taxonomic inputs file. Just rename the outputs.
#First, convert your ASV counts to relative abundances.

dfASVonly<-SkupskyASVonly
dfASVonly <- as.matrix(dfASVonly) #taxa are columns, samples are rows. 

df_relative_ASV <- make_relative(dfASVonly)
df_relative_ASV <- as.data.frame(df_relative_ASV)
Relative_Abundance <- summarize_all(df_relative_ASV, mean)
Relative_Abundance <- t(Relative_Abundance)
#write.csv(relative_abundances, "mucosalSI-Vil-RelativeAbundance.csv") #if you just want a vector of relative abundances

row.names(out4)==row.names(Relative_Abundance)
taxa_names <- cbind(out4,Relative_Abundance)
taxa_names['TaxonNbr']<-NA
ctr=0
for (j in 1:nrow(taxa_names)){
  taxa_names[j,7] <- ctr
  ctr=ctr+1
}
ASV <- row.names(taxa_names)
ASV <- as.data.frame(ASV)
taxa_names_final <- cbind(taxa_names,ASV)
write.csv(taxa_names_final, "AICtable-MucosalSI-SLCTotal.csv") #save your AIC table in case you need it

#Append the taxon name to the Fixed Coefficients output#
out2
out2['Relative_Abundance']<-NA
out2['TaxonNbr']<- NA
out2['ASV'] <- NA

ctr=0
for (j in 1:nrow(out2)){
  if((str_detect(row.names(out2[j,]), "Intercept"))== TRUE){
    ctr=ctr+1
    print(ctr)
    out2[j,] <- cbind(out2[j,1:7],taxa_names_final[ctr,])
  }
  out2[j,8] <- taxa_names_final[ctr,6]
  out2[j,9] <- taxa_names_final[ctr,7]
  out2[j,10] <- taxa_names_final[ctr,8]
}

df<-out2
Taxa_Model_Not_Fit <- df[rowSums(is.na(df)) > 0,]
#View(Taxa_Model_Not_Fit)
Taxa_Model_Did_Fit <- na.omit(df)

#write.csv(Taxa_Model_Not_Fit, "MucosalSI-SLC-Total-FailedTaxa.csv") #If you want to see the taxa rejected later

Taxa_Model_Did_Fit$eEstimate <- exp(Taxa_Model_Did_Fit$Estimate)
Taxa_Model_Did_Fit$log2FoldChange<- log(Taxa_Model_Did_Fit$eEstimate,2)

Cre_Taxa <- Taxa_Model_Did_Fit[grep("^Genotype", row.names(Taxa_Model_Did_Fit)),] #grep your variable of interest 
Cre_Taxa$Qvalues

Cre_Taxa_test<- Cre_Taxa[order(Cre_Taxa$Qvalues),]
Cre_Taxa_test<- filter(Cre_Taxa_test,Qvalues<0.05)

write.csv(Cre_Taxa_test, "MucosalSI-SLC-Total-taxonomicinputs.csv")
