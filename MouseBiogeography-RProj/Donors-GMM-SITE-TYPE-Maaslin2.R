library(Maaslin2)
library(funrar)
library(dplyr)

here::i_am("MouseBiogeography-RProj/Donors-GMM-SITE-TYPE-Maaslin2.R")
input_data <- read.delim(here("Donors-Analysis/omixer_output/Donors_GMM_modules.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
input_metadata <-read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

samples<- input_metadata %>%
  filter(SampleID %in% names(df_input_data)) %>%
  filter(Donor_ID!= "A072") %>%
  pull(SampleID)

df_input_data <- df_input_data[,samples]
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata$SampleID <- row.names(input_metadata)
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)


## Get Luminal and Mucosal datasets to test DC vs all ---
lum_samples <- df_input_metadata %>% filter(Type =="Luminal", SampleID %in% names(df_input_data)) %>% pull(SampleID)
lum_input_data <- df_input_data[, lum_samples]

muc_samples <- df_input_metadata %>% filter(Type =="Mucosal", SampleID %in% names(df_input_data)) %>% pull(SampleID)
muc_input_data <- df_input_data[, muc_samples]

## Get Site datasets to test Lum v Muc ---
input_data <- read.delim(here("Donors-Analysis/omixer_output/Donors_GMM_modules.tsv"), header=TRUE,row.names=1) 
df_input_data<-as.data.frame(input_data)
names(df_input_data)<-gsub("X","",names(df_input_data))
input_metadata <-read.delim(here("Donors-Analysis/starting_files/Donors_Metadata.tsv"),header=TRUE, row.names=1) #mapping file
row.names(input_metadata)
row.names(input_metadata)<-gsub("-",".",row.names(input_metadata))
input_metadata$SampleID <- row.names(input_metadata)

samples<- input_metadata %>%
  filter(SampleID %in% names(df_input_data)) %>%
  pull(SampleID)

df_input_data <- df_input_data[,samples]
target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata$SampleID <- row.names(input_metadata)
df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Donor_ID <- factor(df_input_metadata$Donor_ID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
sapply(df_input_metadata,levels)

si_samples <- df_input_metadata %>% filter(Site_General =="SI", SampleID %in% names(df_input_data)) %>% pull(SampleID)
si_input_data <- df_input_data[, si_samples]

col_samples <- df_input_metadata %>% filter(Site_General =="Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
col_input_data <- df_input_data[, col_samples]

duo_samples <- df_input_metadata %>% filter(Site =="Duodenum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
duo_input_data <- df_input_data[, duo_samples]

jej_samples <- df_input_metadata %>% filter(Site =="Jejunum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
jej_input_data <- df_input_data[, jej_samples]

ile_samples <- df_input_metadata %>% filter(Site =="Ileum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
ile_input_data <- df_input_data[, ile_samples]

cec_samples <- df_input_metadata %>% filter(Site =="Cecum", SampleID %in% names(df_input_data)) %>% pull(SampleID)
cec_input_data <- df_input_data[, cec_samples]

PC_samples <- df_input_metadata %>% filter(Site =="Proximal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
PC_input_data <- df_input_data[, PC_samples]

DC_samples <- df_input_metadata %>% filter(Site =="Distal_Colon", SampleID %in% names(df_input_data)) %>% pull(SampleID)
DC_input_data <- df_input_data[, DC_samples]

## Run Maaslin2 ---

#Luminal 
target <- colnames(lum_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

filepath <- "Donors-Analysis/differential_GMM_site/"
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=lum_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID", "Donor_ID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=lum_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-ColonRef-CLR-Lum-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID", "Donor_ID"),
                    reference=c("Sequencing_Run,Jan_2017","Site,Distal_Colon"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
target <- colnames(muc_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=muc_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite_General-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site_General"), 
                    random_effects = c("MouseID", "Donor_ID"),
                    reference=c("Sequencing_Run,Jan_2017","Site_General,Colon"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum", "Jejunum", "Duodenum"))
fit_data = Maaslin2(input_data=muc_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-ColonRef-CLR-Muc-ComBat-SeqRunSexSite-1-MsID-DonorID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID", "Donor_ID"),
                    reference=c("Sequencing_Run,Jan_2017","Site,Distal_Colon"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)


## Run Maaslin2 ---
filepath <- "Donors-Analysis/differential_GMM_type/"

#Distal_Colon 
target <- colnames(DC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=DC_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-DistalColon-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex","Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017","Type,Luminal"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)



#Proximal_Colon -
target <- colnames(PC_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=PC_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-ProximalColon-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex","Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum 
target <- colnames(cec_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=cec_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
target <- colnames(ile_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=ile_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum -lot of overfitting
target <- colnames(jej_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=jej_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
target <- colnames(duo_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=duo_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Colon
target <- colnames(col_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=col_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-Colon-ComBat-SeqRunSexSiteType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex","Site","Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017","Site,Distal_Colon"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)


#SI
target <- colnames(si_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

fit_data = Maaslin2(input_data=si_input_data, 
                    input_metadata=df_input_metadata, 
                    output = paste0(filepath,"GMM-LumRef-CLR-SI-ComBat-SeqRunSexSiteType-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    reference=c("Sequencing_Run,Jan_2017","Site,Distal_Colon"),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)


