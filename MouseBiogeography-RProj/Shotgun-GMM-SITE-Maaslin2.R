library(Maaslin2)
library(dplyr)

here::i_am("MouseBiogeography-RProj/Shotgun-GMM-SITE-Maaslin2.R")
input_data <- readr::read_delim(here("Shotgun/omixer/Shotgun_GMM_modules.tsv"),delim="\t")
df_input_data<-as.data.frame(input_data)
row.names(df_input_data)<-df_input_data$Module
df_input_data <- select(df_input_data, -c("Module"))
input_metadata <-readr::read_delim(here("Shotgun/BioGeo_Shotgun_Metadata.tsv"),delim="\t") #mapping file
input_metadata <- as.data.frame(input_metadata)
row.names(input_metadata) <- input_metadata$humann_sampleid

target <- colnames(df_input_data)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

input_metadata <- input_metadata
input_metadata$Sequencing_Run <- factor(input_metadata$Sequencing_Run)
input_metadata$MouseID <- factor(input_metadata$MouseID)
input_metadata$Sex <- factor(input_metadata$Sex)
sapply(input_metadata,levels)

## Get UCLA O SPF datasets to test DC vs Jej---
ucla_samples <- input_metadata %>% filter(Dataset =="UCLA_O_SPF", humann_sampleid %in% names(df_input_data)) %>% pull(humann_sampleid)
ucla_input_data <- df_input_data[, ucla_samples]

cs_samples <- input_metadata %>% filter(Dataset =="CS_SPF", humann_sampleid %in% names(df_input_data)) %>% pull(humann_sampleid)
cs_input_data <- df_input_data[, cs_samples]

hum_samples <- input_metadata %>% filter(Dataset =="HUM_Gavage", humann_sampleid %in% names(df_input_data)) %>% pull(humann_sampleid)
hum_input_data <- df_input_data[, hum_samples]

spf_samples <- input_metadata %>% filter(Dataset =="SPF_Gavage", humann_sampleid %in% names(df_input_data)) %>% pull(humann_sampleid)
spf_input_data <- df_input_data[, spf_samples]

## Run Maaslin2 ---

#UCLA 
target <- colnames(ucla_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Jejunum"))
fit_data = Maaslin2(input_data=ucla_input_data, input_metadata=df_input_metadata, 
                    output = here("Shotgun/UCLA_O_SPF/GMM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID"),
                    min_prevalence = 0.15,
                    reference= c('Line,JJWT','Site,Distal_Colon'),
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#CS SPF 
target <- colnames(cs_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Jejunum"))
fit_data = Maaslin2(input_data=cs_input_data, input_metadata=df_input_metadata, 
                    output = here("Shotgun/CS_SPF/GMM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),
                    reference= c('Site,Distal_Colon'),
                    min_prevalence = 0.15,
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#HUM Gavage
target <- colnames(hum_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Jejunum"))
fit_data = Maaslin2(input_data=hum_input_data, input_metadata=df_input_metadata, 
                    output = here("Shotgun/CS_SPF/GMM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),
                    reference= c('Site,Distal_Colon'),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    plot_heatmap = FALSE,plot_scatter = FALSE)

#SPF Gavage
target <- colnames(spf_input_data)
df_input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(df_input_metadata)

df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Jejunum"))
fit_data = Maaslin2(input_data=spf_input_data, input_metadata=df_input_metadata, 
                    output = here("Shotgun/CS_SPF/GMM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID"), 
                    fixed_effects = c("Sequencing_Run","Sex", "Site"), 
                    random_effects = c("MouseID"),
                    reference= c('Site,Distal_Colon'),
                    normalization="clr", 
                    transform ="none",
                    min_prevalence = 0.15,
                    plot_heatmap = FALSE,plot_scatter = FALSE)

