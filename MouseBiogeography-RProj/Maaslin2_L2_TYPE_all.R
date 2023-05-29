library(Maaslin2)
library(funrar)
library(dplyr)
library(tidyr)

here::i_am("MouseBiogeography-RProj/Maaslin2_L2_TYPE_all.R")
here::here()

### Type Comparisons at the Phylum Level---
## UCLA O. SPF --
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Duodenum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Jejunum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Ileum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Cecum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - PC_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - DC_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - Colon_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/UCLA_O_SPF/Maaslin2 SITE and TYPE L2 - SI_L2.csv", header=TRUE, row.names=1) 

df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))

input_metadata <-read.csv("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/Maasllin2 Site Genus Level/Regional-Metadata-All.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-Colon-CLR-ComBat-SeqRunLineSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-SI-CLR-ComBat-SeqRunLineSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site","Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-Duodenum-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-Jejunum-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-Ileum-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-Cecum-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-PC-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/UCLA_O_SPF/L2-DC-CLR-ComBat-SeqRunLineSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Type"), 
                    random_effects = c("MouseID_Line"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

## CS SPF --
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - Colon_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - SI_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - Duodenum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - Jejunum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - Ileum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - Cecum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - PC_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/CS_SPF/Type_L2/Maaslin2_L2 SITE and TYPE  - DC_L2.csv", header=TRUE, row.names=1) 

df_input_data<-as.data.frame(input_data)
df_input_data <- df_input_data %>% select(-c("taxonomy"))

input_metadata <-read.csv("CS-Facility-Analysis/CS_Facility_Metadata.csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-Colon-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-SI-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/CS_SPF/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

## SPF Gavage --
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Colon.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-SI.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Duo.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Jejunum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Ileum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-Cecum.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-PC.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/SPF_Gavage/Maaslin2_Type_L2/SPF_Maaslin2_L2_SITE_and_TYPE - L2-DC.csv", header=TRUE, row.names=1) 

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))

input_metadata <-read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Colon-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-SI-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Duodenum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/SPF_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
## HUM Gavage --
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Colon_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - SI_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Duodenum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Jejunum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Ileum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - Cecum_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - PC_L2.csv", header=TRUE, row.names=1) 
input_data <- read.csv("Maaslin2_L2/HUM_Gavage/HUM_Maaslin2_SITE_and_TYPE_L2 - DC_L2.csv", header=TRUE, row.names=1) 

df_input_data <- as.data.frame(input_data)
df_input_data <- select(df_input_data, -c("taxonomy"))
input_metadata <-read.csv("Humanized-Biogeography-Analysis/Humanized Metadata - All-Humanized-Metadata (1).csv",header=TRUE, row.names=1) #mapping file

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)
names(input_metadata)

df_input_metadata<-input_metadata
df_input_metadata$MouseID <- factor(df_input_metadata$MouseID)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))

sapply(df_input_metadata,levels)

#Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Colon-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SI
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-SI-CLR-ComBat-SeqRunSexSiteType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Site","Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Duodenum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Duodenum-CLR-ComBat-SexType-1-MsID", 
                    fixed_effects = c("Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Jejunum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Jejunum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Ileum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Ileum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Cecum
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-Cecum-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#Proximal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-PC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
#Distal Colon
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, 
                    output = "Maaslin2_L2/HUM_Gavage/L2-DC-CLR-ComBat-SeqRunSexType-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Sex", "Type"), 
                    random_effects = c("MouseID"),
                    normalization="clr", transform ="none",
                    plot_heatmap = FALSE,plot_scatter = FALSE)
