### Assign Taxonomy after running DeBlur for denoising reads: Donors ---
library(here)
library(sva)
library(dplyr)

BiocManager::install("sva")

here::i_am("MouseBiogeography-RProj/Donors_ASV_Preprocessing.R")

## Combine taxonomy assignments --
PE_taxonomy <- readr::read_csv(here("Donors-Analysis/starting_files/Donors_PE_taxonomy_assignments.csv"))
SE_taxonomy <- readr::read_csv(here("Donors-Analysis/starting_files/Donors_single_taxonomy_assignments.csv"))

Donors_combined_taxonomy <- rbind(PE_taxonomy,SE_taxonomy)
Donors_combined_taxonomy$ASV <- Donors_combined_taxonomy$...1
Donors_combined_taxonomy <- Donors_combined_taxonomy %>% select(c("ASV","taxonomy"))
Donors_combined_taxonomy <- unique(Donors_combined_taxonomy)
readr::write_csv(Donors_combined_taxonomy,here("Donors-Analysis/Donors_combined_taxonomy_assignments.csv"))

## Combine aligned DNA sequences and make a taxonomy file --
PE_aligned <- readr::read_csv(here("Donors-Analysis/Donors_PE_aligned_DNA_sequences.csv"))
SE_aligned <- readr::read_csv(here("Donors-Analysis/Donors_single_aligned_DNA_sequences.csv"))

Donors_combined_aligned <- rbind(PE_aligned,SE_aligned)
Donors_combined_aligned <- Donors_combined_aligned[,-1]
Donors_combined_aligned <- unique(Donors_combined_aligned)
Donors_combined_aligned$ASV <- Donors_combined_aligned$`fasta3[-1, ]`
readr::write_csv(Donors_combined_taxonomy,here("Donors-Analysis/Donors_combined_DNA_sequences.csv"))

## Make a taxonomy file for importing into QIIME -- 
Donors_taxonomy_for_qza <- merge(Donors_combined_taxonomy, Donors_combined_aligned, by="ASV")
Donors_taxonomy_for_qza <- Donors_taxonomy_for_qza %>% select(c("ASV","taxonomy"))
Donors_taxonomy_for_qza <- rename(Donors_taxonomy_for_qza, c("Feature ID" = "ASV"))
Donors_taxonomy_for_qza <- rename(Donors_taxonomy_for_qza, c("Taxon" = "taxonomy"))

#Custom function to calculate the length of the taxonomy string
get_taxonomy_length <- function(taxonomy) {
  return(nchar(as.character(taxonomy)))
}

# Find rows with duplicate uniqueIDs and keep the one with the longer taxonomy string
result <- Donors_taxonomy_for_qza %>%
  group_by(`Feature ID`) %>%
  filter(n_distinct(Taxon) > 1) %>%
  slice(which.max(sapply(Taxon, get_taxonomy_length))) %>%
  ungroup()
print(result$Taxon)

# Remove the entries in the OG dataframe that had duplicate uniqueIDs and then rbind the deduped uniqueIDs back
filtered_df <- anti_join(Donors_taxonomy_for_qza, result, by = "Feature ID")
Donors_taxonomy_for_qza <- rbind(filtered_df, result)

readr::write_delim(Donors_taxonomy_for_qza,here("Donors-Analysis/Donors_taxonomy.tsv"),delim="\t")

## Combine tables --
PE_table <- readr::read_delim(here("Donors-Analysis/export_Donors-paired-end-table/feature-table.tsv"),delim="\t")
SE_table <- readr::read_delim(here("Donors-Analysis/export_Donors-single-filtered-table/feature-table.tsv"),delim="\t")

row.names(PE_table) <- PE_table$OTU.ID
row.names(SE_table) <- SE_table$OTU.ID

# Merge the data frames based on ASV and assign 0 where these are lacking
merged_df <- merge(PE_table, SE_table, by = "OTU.ID", all = TRUE)
merged_df[is.na(merged_df)] <- 0

rownames(merged_df) <- merged_df$OTU.ID
merged_df <- merged_df[, -1]

# Filter samples with total counts greater than or equal to 10,000
total_counts <- colSums(merged_df)
print(total_counts)
selected_samples <- names(total_counts[total_counts >= 10000])
rejected_samples <- names(total_counts[total_counts < 10000])

#Use filtered ASV and metadata files to then split the ASV into each sequencing run, 
#apply 15% prevalence feature filtering threshold, and select only intersecting ASVs
ASV <- merged_df[,selected_samples]
metadata <- readr::read_csv(here("Donors-Analysis/Donors_Metadata.csv"))
readr::write_csv(ASV[1,],here(("Donors-Analysis/Donors_ASV_names.csv")))

metadata$SampleID <- gsub("_","-",metadata$SampleID)

target <- names(ASV)
metadata <- metadata[match(target, metadata$SampleID),]
target == metadata$SampleID

metadata <- metadata %>%
  filter(SampleID %in% names(ASV)) %>%
  filter(MouseID!= "U2") %>% 
  filter(Original_Human_Stool=="N")
names(ASV)
metadata$SampleID


#Filter by Sequencing Run 
filter_by_metadata <- function(metadata, asv, run) {
  samples <- metadata %>%
    filter(Sequencing_Run == run, SampleID %in% names(asv)) %>%
    pull(SampleID)
  
  asv_data <- asv[, samples]
  return(asv_data)
}

metadata$Sequencing_Run
srun1 <- filter_by_metadata(metadata,ASV,"April_2017")
srun2 <- filter_by_metadata(metadata,ASV,"Dec_2017")
srun3 <- filter_by_metadata(metadata,ASV,"Dec_2018")
srun4 <- filter_by_metadata(metadata,ASV,"Jan_2017")
srun5 <- filter_by_metadata(metadata,ASV,"June_2016")
srun6 <- filter_by_metadata(metadata,ASV,"May_2018")
srun7 <- filter_by_metadata(metadata,ASV,"Oct_2017")

# Prevalence filter each sequencing run 
prevalence_filter <- function(counts_df, min_sample_number){
  
  lumcol_counts<-as.data.frame(counts_df)
  t_df_input_data<-as.data.frame(t(lumcol_counts))
  
  ctr= 0
  prevalence <- vector(mode="numeric")
  
  for(i in 1:ncol(t_df_input_data)){
    v<-t_df_input_data %>% pull(i)
    for(j in 1:length(v)){
      if (t_df_input_data[j,i]>0){
        ctr=1+ctr
      }
      else {
        ctr=ctr
      }
    }
    prevalence<-append(prevalence,ctr)
    ctr=0
  }
  
  lumcol_counts$prevalence<-prevalence 
  lumcol_counts<- lumcol_counts%>% filter(prevalence>=min_sample_number) 
  lumcol_counts <- lumcol_counts %>% select(-c(prevalence))
  return(lumcol_counts)
}

0.15*14
srun1_prev <- prevalence_filter(srun1, 2)
0.15*107
srun2_prev <- prevalence_filter(srun2, 16)
0.15*130
srun3_prev <- prevalence_filter(srun3, 20)
0.15*128
srun4_prev <- prevalence_filter(srun4, 19)
0.15*12
srun5_prev <- prevalence_filter(srun5, 2)
0.15*82
srun6_prev <- prevalence_filter(srun6, 12)
0.15*70
srun7_prev <- prevalence_filter(srun7, 11)

##  Find the intersection of all features then query target vector against each dataset, then merge datasets by feature
# exclude srun1 (April 2017) and all samples where mouseID = U2, and also all human fecal samples
common_row_names <- rownames(srun2_prev)
for (df in list( srun3_prev, srun4_prev, srun6_prev, srun7_prev)) {
  common_row_names <- intersect(common_row_names, rownames(df))
}

print(common_row_names)

## Query the intersection against each dataset 
srun2<-srun2_prev[row.names(srun2_prev) %in% common_row_names,]
srun3<-srun3_prev[row.names(srun3_prev) %in% common_row_names,]
srun4<-srun4_prev[row.names(srun4_prev) %in% common_row_names,]
srun6<-srun6_prev[row.names(srun6_prev) %in% common_row_names,]
srun7<-srun7_prev[row.names(srun7_prev) %in% common_row_names,]

combined_df<- srun2
for (df in list( srun3, srun4, srun6, srun7)) {
  combined_df <- cbind(combined_df,df)
}

### REDO WITH PREV FILTERING ON WHOLE DATASET ---
### Prevalence filtering on all mouse stool samples (exclude April 2017 and MouseID U2) and do combat seq
metadata <- readr::read_csv(here("Donors-Analysis/Donors_Metadata.csv"))
readr::write_delim(metadata, here("Donors-Analysis/Donors_Metadata.tsv"),delim="\t")
target <- names(ASV)
metadata <- metadata[match(target, metadata$SampleID),]
target == metadata$SampleID

metadata <- metadata %>%
  filter(SampleID %in% names(ASV)) %>%
  filter(MouseID!= "U2") %>% 
  filter(Original_Human_Stool=="N")
names(ASV)
metadata$SampleID


#Filter ASV -
samples <- metadata %>%
    filter(SampleID %in% names(ASV)) %>%
    pull(SampleID)
  
mouse_ASV <- ASV[, samples]

0.15*531

mouse_ASV_prev <- prevalence_filter(mouse_ASV,80)

#ComBatSeq - 
batch<- as.character(metadata$Sequencing_Run)
metadata$Type <- factor(metadata$Type)
metadata$Site <- factor(metadata$Site)
modcombat=model.matrix(~Type + Site,data=metadata)

input=as.matrix(mouse_ASV_prev)
combat_adjusted_counts=sva::ComBat_seq(input,batch=batch,group=NULL,covar_mod=modcombat)  
combat_adjusted_counts <- as.data.frame(combat_adjusted_counts)

# Replace QIIME sequences with actual ASV sequences 
combat_adjusted_counts$QIIME_seqs <- rownames(combat_adjusted_counts)
taxonomy <- Donors_combined_aligned %>% select(c("QIIME_seqs","ASV"))
combat_adjusted_counts_ASV <-merge(combat_adjusted_counts, taxonomy, by="QIIME_seqs")
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>% select(-c("QIIME_seqs"))
combat_adjusted_counts_ASV <- combat_adjusted_counts_ASV %>%
  select(ASV, everything())
combat_adjusted_counts_ASV <- rename(combat_adjusted_counts_ASV, c("#OTU.ID" = "ASV"))

readr::write_delim(combat_adjusted_counts_ASV,here("Donors-Analysis/Donors-Mice-1xPrev0.15-ComBat-ASV.tsv"), delim="\t")

### Utilize the same set of ASVs from the mouse stool on the human stool ---
target <- combat_adjusted_counts_ASV$`#OTU.ID`
metadata <- readr::read_csv(here("Donors-Analysis/Donors_Metadata.csv"))

mouse_metadata <- metadata %>%
  filter(SampleID %in% names(ASV)) %>%
  filter(MouseID!= "U2") %>% 
  filter(Original_Human_Stool=="N")
metadata <- metadata %>%
  filter(SampleID %in% names(ASV)) %>%
  filter(Original_Human_Stool=="Y")%>%
  filter(Donor_ID %in% mouse_metadata$Donor_ID)

metadata$Donor_ID
donors <- unique(mouse_metadata$Donor_ID)
setdiff(donors, metadata$Donor_ID) # missing A072 sequence data

names(ASV)
metadata$SampleID

samples <- metadata %>%
  pull(SampleID)

hoomins_ASV <- ASV[, samples]
hoomins_ASV$QIIME_seqs <- rownames(hoomins_ASV)
taxonomy <- Donors_combined_aligned %>% select(c("QIIME_seqs","ASV"))
hoomins_ASV <-merge(hoomins_ASV, taxonomy, by="QIIME_seqs")
hoomins_ASV <- hoomins_ASV %>% select(-c("QIIME_seqs"))

hoomins_ASV <- hoomins_ASV %>%
  select(ASV, everything())
hoomins_ASV <- rename(hoomins_ASV, c("#OTU.ID" = "ASV"))


hoomins_ASV<-hoomins_ASV[hoomins_ASV$`#OTU.ID` %in% target,]
readr::write_delim(hoomins_ASV,here("Donors-Analysis/Donors-Hoomins-ASV.tsv"), delim="\t")
