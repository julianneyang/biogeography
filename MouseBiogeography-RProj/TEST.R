library(dplyr)
library(tidyverse)
library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)

#########Removing all taxa which have no assignment############################
df <- read.csv("RTP_dysbiosis_ASV - 8to15_unaffected.tsv", header=TRUE)
colnames(df)#lalala
View(df)
niceASV <- subset(df, select=-c(None) )
niceASV <- df[,!grepl("None", colnames(df))]
View(niceASV)
write.csv(niceASV, "Skupsky_Final_ASV.csv")
###################################################################################################

#########################WHAT TO DO AFTER RUNNING GLMMTMB#########################################

setwd("C:/Users/Jacobs Laboratory/Documents/JCYang/2021 Vil1Cre Re-analysis/2021-6 Vil1Cre Big Batch/")
#################First, Append the taxon name to the Fixed Coefficients output################################ 
df <- read.csv("vil-nonzero_SiteSexCre_1-Litter-MsID_LuminalColon_GLMMTMB_FixedCoeff.csv", row.names=1)
df <- read.csv("vil-nonzero_SiteSexCre_1-Litter-MsID_luminal_SI_GLMMTMB_FixedCoeff.csv", row.names=1)
df <- read.csv("vil-nonzero_SiteSexCre_1-Litter-MsID_mucosal-colon_GLMMTMB_FixedCoeff.csv", row.names=1)
df <- read.csv("vil-nonzero_SiteSexCre_1-Litter-MsID_mucosalSI_GLMMTMB_FixedCoeff.csv", row.names=1)

df['Taxon'] <- NA
df['TaxonNbr']<- NA
View(df)
ctr=0
#View(df)
taxa_names <- read.csv("LuminalColon-Taxonomy_Key.csv", header= TRUE) #all the taxa names. Get this from the AIC output csv file. 
taxa_names <- read.csv("LuminalSI-Taxonomy-Key.csv", header= TRUE) #all the taxa names. Get this from the AIC output csv file. 
taxa_names <- read.csv("MucosalColon-Taxonomy-Key.csv", header= TRUE) #all the taxa names. Get this from the AIC output csv file. 
taxa_names <- read.csv("MucosalSI-Taxonomy-Key.csv", header= TRUE) #all the taxa names. Get this from the AIC output csv file. 

taxa_names <- as.data.frame(taxa_names)
for (j in 1:nrow(df)){
  if((str_detect(row.names(df[j,]), "Intercept"))== TRUE){
    ctr=ctr+1
    print(ctr)
    df[j,] <- cbind(df[j,1:9],taxa_names[ctr,])
  }
  df[j,8] <- taxa_names[ctr,1]
  df[j,9] <- taxa_names[ctr,2]
}
View(df)
warnings()

Taxa_Model_Not_Fit <- df[rowSums(is.na(df)) > 0,]
#View(Taxa_Model_Not_Fit)
Taxa_Model_Did_Fit <- na.omit(df)
#View(Taxa_Model_Did_Fit)
write.csv(Taxa_Model_Did_Fit, "MucosalSI-Vil-SuccessTaxa.csv")
write.csv(Taxa_Model_Not_Fit, "MucosalSI-Vil-FailedTaxa.csv")

df_nonzero <- read.csv("nonzero_SexGen_1-Study_lumen_Colon_GLMMTMB_FixedCoeff.csv", row.names=1)

df_nonzero['Taxon'] <- NA
df_nonzero['TaxonNbr']<- NA
#View(df_nonzero)
ctr=0
taxa_names <- read.csv("All_Taxa_Names_Nbrs.csv", header= TRUE) #all the taxa names. Get this from the AIC output csv file. 
taxa_names <- as.data.frame(taxa_names)
for (j in 1:nrow(df_nonzero)){
  if((str_detect(row.names(df[j,]), "Intercept"))== TRUE){
    ctr=ctr+1
    print(ctr)
    df[j,] <- cbind(df[j,1:6],taxa_names[ctr,])
  }
  df_nonzero[j,7] <- taxa_names[ctr,1]
  df_nonzero[j,8] <- taxa_names[ctr,2]
}
View(df_nonzero)

####################Next, Keep just the taxa which correspond to a lower AIC in Zero Inflated NB model#########################
df_lowerAIC <- read.csv("zero_AIC.csv", row.names = 1)
#View(df_lowerAIC)
df_subset <- subset(df, df$TaxonNbr %in% df_lowerAIC$TaxonNbr) 
View(df_subset)

#later go and run this for your non-zero inflated model. 

df_lowerAIC_nonzero <- read.csv("nonzero_AIC.csv")
#View(df_lowerAIC_nonzero)
df_subset_nonzero <- subset(df_nonzero, df_nonzero$TaxonNbr %in% df_lowerAIC_nonzero$TaxonNbr) 
View(df_subset_nonzero)


####################Then, subset out the taxa which have NA for the std Error################################################

Taxa_Model_Not_Fit <- df_subset[rowSums(is.na(df_subset)) > 0,]
#View(Taxa_Model_Not_Fit)
Taxa_Model_Did_Fit <- na.omit(df_subset)
View(Taxa_Model_Did_Fit)


write.csv(Taxa_Model_Not_Fit, "Zero_SexGen_1-Study_1-Study_lumencolon_FailedTaxa.csv")
write.csv(Taxa_Model_Did_Fit, "Zero_SexGen_1-Study_1-Study_lumencolon__SuccessTaxa.csv")
#later go and run this for your non-zero inflated model.
Taxa_Model_Not_Fit <- df_subset_nonzero[rowSums(is.na(df_subset_nonzero)) > 0,]
View(Taxa_Model_Not_Fit)
Taxa_Model_Did_Fit <- na.omit(df_subset_nonzero)
View(Taxa_Model_Did_Fit)

write.csv(Taxa_Model_Not_Fit, "Nonzero_SexGen_1-Study_lumencolon_FailedTaxa.csv")
write.csv(Taxa_Model_Did_Fit, "Nonzero_SexGen_1-Study_lumencolon_SuccessTaxa.csv") #all succeeded
#From here, take your Success Taxa and sort by BH p-values. Then, sort all entries with BH p-values <0.05 by TaxonNbr. 
#You now have all your taxa which are significantly differentially abundant and fit to the zero inf model! now fit to the other model. 

###################Then, convert the ASV table to relative abundances. You will use this later for your taxonomic inputs plot
#First, convert your ASV counts to relative abundances.
df<- read.csv("LuminalColon-Vil-GLMMTMB-inputs.csv", row.names=1 , header = TRUE)
df<- read.csv("LuminalSI-Vil-GLMMTMB-inputs.csv", row.names=1 , header = TRUE)
df<- read.csv("MucosalColon-Vil-GLMMTMB-inputs.csv", row.names=1 , header = TRUE)
df<- read.csv("MucosalSI-Vil-GLMMTMB-inputs.csv", row.names=1 , header = TRUE)


df<- as.data.frame(df)
dfASVonly = df[, -c(451:469)] #luminalcolon vil
dfASVonly = df[,-c(334:351)] #luminalSI vil
dfASVonly = df[,-c(213:230)] #mucosalColon
dfASVonly = df[,-c(80:97)] #mucosalSI
dfASVonly <- as.matrix(dfASVonly) #taxa are columns, samples are rows. 
make_relative(dfASVonly)
df_relative_ASV <- make_relative(dfASVonly)
df_relative_ASV <- as.data.frame(df_relative_ASV)
relative_abundances <- summarize_all(df_relative_ASV, mean)
relative_abundances <- t(relative_abundances)
write.csv(relative_abundances, "mucosalSI-Vil-RelativeAbundance.csv")

#Subset out all the relative abundance counts for significantly differentially abundant taxa
significant_diff_abundant_taxa <- read.csv("RTP_dysbiosis_15plus - 15plus_Dysbiosis_SignificantTaxa.csv", row.names = 1)
significant_diff_abundant_taxa$XTaxonNbr <- paste("X", significant_diff_abundant_taxa$TaxonNbr, sep="")
significant_diff_abundant_taxa$XTaxonNbr
transposed_relative_ASV <- t(df_relative_ASV)
transposed_relative_ASV <- as.data.frame(transposed_relative_ASV)
row.names(transposed_relative_ASV)
relative_ASV_subset <- subset(transposed_relative_ASV, row.names(transposed_relative_ASV) %in% significant_diff_abundant_taxa$XTaxonNbr) 
write.csv(relative_ASV_subset, "Significant_15plus_RelativeAbundances.csv", sep =",")

#Subset even Further, if needed 
diet_taxa <- read.csv("ZINBMM - Diet_Significantly_Diff_Abundant_Taxa.csv", row.names=1)
diet_relative_ASV_subset <- subset(relative_ASV_subset, row.names(relative_ASV_subset)  %in% diet_taxa$Taxon)
write.csv(diet_relative_ASV_subset, "ASV_DietOnly_Significantly_Differentially_Abundant_Taxa.csv")

##########Get absolute abundances also, subsetted by parameter of interest (in this case, it's Diet). You only need to do this if you're planning to make a heatmap
transposed_absolute_ASV <- t(dfASVonly)
transposed_absolute_ASV <- as.data.frame(transposed_absolute_ASV)
absolute_ASV_subset <- subset(transposed_absolute_ASV, row.names(transposed_absolute_ASV) %in% significant_diff_abundant_taxa$Taxon) 
write.csv(absolute_ASV_subset, "ASV_All_Significantly_Differentially_Abundant_Taxa_absolute.csv")
diet_absolute_ASV_subset <- subset(absolute_ASV_subset, row.names(absolute_ASV_subset)  %in% diet_taxa$Taxon)
write.csv(diet_absolute_ASV_subset, "ASV_DietOnly_Significantly_Differentially_Abundant_Taxa_absolute.csv")




