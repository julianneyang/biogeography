library(renv)
devtools::install_github("biobakery/melonnpan")
renv::restore()
renv::snapshot()
library(melonnpan)
library(here)
library(stringr)

here::i_am("MouseBiogeography-RProj/MelonnPan.R")

## Functions --
impute.knn.obs.sel <- function(dat, K=10) { #rownames are samples, columns are features
  
  results <- list()
  cor.cutoff <- 0.2     # use only variables with cor>0.2 for distance computation
  
  da1 <- dat 
  da1list <- da2list <- rep(list(dat),length(K)) 
  
  incom.vars <- which(apply(dat,2,function(x) any(is.na(x))))
  incom.obs <- which(apply(da1,1,function(x) any(is.na(x))))
  
  Cor <- cor(da1,use="p")
  
  D2list <- lapply(incom.vars, function(j) {
    varsel <- which(abs(Cor[j,])>cor.cutoff)  
    if(length(varsel)>10) varsel <- order(abs(Cor[j,]),decreasing=T)[1:11]
    if(length(varsel)<5) varsel <- order(abs(Cor[j,]),decreasing=T)[1:6]
    D2 <- as.matrix(dist(scale(da1[,varsel])),upper=T,diag=T) 
    if(any(is.na(D2))) {
      D2a <- as.matrix(dist(scale(da1)),upper=T,diag=T)*sqrt(length(varsel)/ncol(da1)) 
      D2[is.na(D2)] <- D2a[is.na(D2)] 
    }
    diag(D2) <- NA
    D2})
  names(D2list) <- incom.vars
  
  for (i in incom.obs){
    comvars <-  complete.cases(as.numeric(da1[i,]))
    for (j in which(!comvars)) {
      D2 <- D2list[[as.character(j)]]                                 
      if(any(!is.na(D2[i,]))) {
        KNNids <- order(D2[i,],na.last=NA)
        KNNids_naomit <- KNNids[sapply(KNNids,function(ii) any(!is.na(da1[ii,j])))] 
      } else {
        KNNids  <- NULL
      }
      da1list <- lapply(1:length(da1list),function(ii) {
        k <- K[ii] 
        da <-  da1list[[ii]]
        if(!is.null(KNNids)) {
          KNNids_sel <- intersect(KNNids[1:min(k,length(KNNids))],KNNids_naomit)
        }
        if(length(KNNids_sel)<1) {
          KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]
        } else if (length(which(sapply(KNNids_sel,function(ii) !is.na(da1[ii,j])))) < floor(k/2) ){
          KNNids_sel <- KNNids_naomit[1:min(floor(k/2),length(KNNids_naomit))]} 
        
        if(any(!is.na(D2[i,])) & length(KNNids)>=1) {
          da_sel <- da[KNNids_sel,j]
          da[i,j] <- sum(da_sel*exp(-D2[i,KNNids_sel]),na.rm=T)/sum(exp(-D2[i,KNNids_sel])[!is.na(da_sel)],na.rm=T) }
        da}) 
    }
  }
  da1list <- lapply(da1list, function(da) {
    da <- apply(da,2, function(x) {
      if(any(is.na(x))) x[is.na(x)] <- mean(x,na.rm=T)
      x}) 
    da})
  
  results <- c(results,list(da1list))
  names(results)[length(results)] <- "knn.sample.euc.sel" 
  rm(da1list,da2list)  
  return(results$knn.sample.euc.sel[[1]])
}

## Read in TL1A file ---
tl1a<-read.delim(here("../melonnpan_data/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),row.names=1)
names(tl1a) <- gsub("\\.DC\\.","_",names(tl1a)) 
names(tl1a) <- gsub("\\.FP\\.","_",names(tl1a))
names(tl1a) <- gsub("_DC_","_",names(tl1a))
names(tl1a) <- gsub("_FP_","_",names(tl1a))
names(tl1a) <- gsub("merged.*TLIA","TL1A", names(tl1a))
names(tl1a) <- gsub("_S.*__kneaddata_paired_Abundance.RPKs","",names(tl1a))
names(tl1a) <- gsub("__kneaddata_paired_Abundance.RPKs","",names(tl1a))


lipid_tl1a <- read.csv(here("../melonnpan_data/Metabolomics_TL1A - TL1A_lipidomics.csv"))
names(lipid_tl1a)
lipid_tl1a[lipid_tl1a=="na"] <- NA
lipid_tl1a[402,1] <- "PE 34:1 repeat"
lipid_tl1a <- lipid_tl1a %>% filter(name!="") %>%
  column_to_rownames(var="name") 
lipid_tl1a <- lipid_tl1a %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))
lipid_tl1a <- as.data.frame(t(lipid_tl1a))
lipid_tl1a <- lipid_tl1a%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
lipid_tl1a <- as.data.frame(t(lipid_tl1a))

#lipid_tl1a <- lipid_tl1a %>%
  #mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%  # Convert everything to numeric
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values


primary_tl1a <- read.csv(here("../melonnpan_data/Metabolomics_TL1A - TL1A_primary.csv"))
primary_tl1a <- primary_tl1a %>% filter(is.na(as.numeric(name))) %>% 
  column_to_rownames(var="name") %>%
  mutate(across(everything(), ~ as.numeric(as.character(.))))  # Convert everything to numeric
primary_tl1a <- as.data.frame(t(primary_tl1a))
primary_tl1a <- primary_tl1a%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values
primary_tl1a <- as.data.frame(t(primary_tl1a))

amines_tl1a <- read.csv(here("../melonnpan_data/Metabolomics_TL1A - TL1A_amines.csv"))
amines_tl1a[36,1] <- "2 -Deoxycytidine Repeat"
amines_tl1a[389,1] <- "Pantothenic acid Repeat"
amines_tl1a<- amines_tl1a%>% filter(name!="") %>%
  column_to_rownames(var="name") %>% 
  mutate(across(everything(), ~ as.numeric(as.character(.))))  # Convert everything to numeric
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values
amines_tl1a <- as.data.frame(t(amines_tl1a))
amines_tl1a <- amines_tl1a%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
amines_tl1a <- as.data.frame((t(amines_tl1a)))

# Find common column names
common_cols <- Reduce(intersect, list(colnames(lipid_tl1a), 
                                      colnames(primary_tl1a), 
                                      colnames(amines_tl1a),
                                      colnames(tl1a)))

colnames(tl1a) %in% rownames(amines_tl1a)
?sapply()

# Convert columns with the same name to character type before binding
lipid_tl1a[common_cols] <- lapply(lipid_tl1a[common_cols], as.character)
primary_tl1a[common_cols] <- lapply(primary_tl1a[common_cols], as.character)
amines_tl1a[common_cols] <- lapply(amines_tl1a[common_cols], as.character)

# Combine dataframes with only common columns
tl1a_combined_df <- bind_rows(
  lipid_tl1a %>% select(all_of(common_cols)),
  primary_tl1a %>% select(all_of(common_cols)),
  amines_tl1a %>% select(all_of(common_cols))
)
#write.csv(tl1a_combined_df, here("../melonnpan_data/Tl1A_combined_metabolites.csv"))

tl1a_combined_relative_abundance <- tl1a_combined_df %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
#write.csv(tl1a_combined_relative_abundance, here("../melonnpan_data/TL1A_combined_relabun.csv"))

# Select only the samples with paired microbiome data 
tl1a_microbiome <- tl1a %>% select(all_of(common_cols))
tl1a_microbiome <- tl1a_microbiome %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
names(tl1a_microbiome)==names(tl1a_combined_df)

### Read in input files ---
df<-read.delim(here("../melonnpan_data/IL10_merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),row.names=1)
names(df) <- gsub("merged_.*IL10FP","IL",names(df)) 
names(df) <- gsub("_S.*__kneaddata_paired_Abundance.RPKs","",names(df))
names(df) <- gsub("_Rep2", "", names(df))
names(df) <- gsub("_S.*","",names(df))
names(df) <- gsub(".Rep2","",names(df))
names(df)

lipid <- read.csv(here("../melonnpan_data/Metabolomics_IL10_Lipid.csv"))
lipid[lipid=="na"] <- NA
lipid <- lipid %>% filter(name!="") %>%
  column_to_rownames(var="name") 
lipid <- lipid %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) #%>%  # Convert everything to numeric
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values
lipid <- as.data.frame(t(lipid))
lipid <- lipid%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
lipid <- as.data.frame(t(lipid))


primary <- read.csv(here("../melonnpan_data/Metabolomics_IL10_Primary.csv"))
primary <- primary %>% filter(is.na(as.numeric(name))) %>% 
  column_to_rownames(var="name") %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) #%>%  # Convert everything to numeric
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values
primary <- as.data.frame(t(primary))
primary <- primary%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
primary <- as.data.frame(t(primary))

amines <- read.csv(here("../melonnpan_data/Metabolomics_IL10_Amine.csv"))
amines <- amines %>% filter(name!="") %>%
  column_to_rownames(var="name") 
amines <- amines %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) #%>%  # Convert everything to numeric
  #mutate(across(everything(), ~ ifelse(is.na(.), median(., na.rm = TRUE), .)))  # Impute missing values
amines <- as.data.frame(t(amines))
amines <- amines%>%
  impute.knn.obs.sel() %>%
  as.data.frame()
amines <- as.data.frame(t(amines))


# Find common column names
common_cols <- Reduce(intersect, list(colnames(lipid), colnames(primary), colnames(amines)))

# Convert columns with the same name to character type before binding
lipid[common_cols] <- lapply(lipid[common_cols], as.character)
primary[common_cols] <- lapply(primary[common_cols], as.character)
amines[common_cols] <- lapply(amines[common_cols], as.character)

# Combine dataframes with only common columns
combined_df <- bind_rows(
  lipid %>% select(all_of(common_cols)),
  primary %>% select(all_of(common_cols)),
  amines %>% select(all_of(common_cols))
)
write.csv(combined_df, here("../melonnpan_data/IL10_combined_metabolites.csv"))


combined_relative_abundance <- combined_df %>%
  mutate(across(everything(), ~ as.numeric(as.character(.)))) %>%
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
write.csv(combined_relative_abundance, here("../melonnpan_data/IL10_combined_relabun.csv"))

# Select only the samples with paired microbiome data 
microbiome <- df %>% select(all_of(common_cols))
microbiome <- microbiome %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
names(microbiome)==names(combined_df)

### Merge metabolome data and microbiome data ---
combined_relative_abundance <- read.csv(here("../melonnpan_data/IL10_combined_relabun.csv"),row.names = 1)
tl1a_combined_relative_abundance <- read.csv(here("../melonnpan_data/TL1A_combined_relabun.csv"), row.names=1)

rownames(tl1a_combined_relative_abundance) <- rownames(combined_relative_abundance)
complete_metabolome_training_data <- cbind(tl1a_combined_relative_abundance, combined_relative_abundance)

tl1a_microbiome$KO <- row.names(tl1a_microbiome)
microbiome$KO <- row.names(microbiome)

complete_microbiome_training_data <- merge(tl1a_microbiome, microbiome, by="KO")
complete_microbiome_training_data <- complete_microbiome_training_data %>% 
  column_to_rownames("KO")

### Perform additional feature filtering ---

filter_features <- function(df, threshold = 0.0001, min_percentage = 0.10) {
  # Store the row names (features)
  #feature_names <- rownames(df)
  
  # Calculate the number of samples where the feature is above the threshold
  num_samples_above_threshold <- rowSums(df > threshold)
  
  # Find the minimum number of samples required based on the percentage
  min_samples <- ncol(df) * min_percentage
  
  # Filter out features that are below the threshold in ≥10% of samples
  filtered_df <- df[num_samples_above_threshold >= min_samples, ]
  
  # Reassign the row names to the filtered dataframe
  #rownames(filtered_df) <- feature_names[num_samples_above_threshold >= min_samples]
  
  return(filtered_df)
}


filtered_microbiome_data <- filter_features(complete_microbiome_training_data, threshold = 0.000001, min_percentage = 0.10)
final_microbiome <- as.data.frame(t(filtered_microbiome_data))
final_microbiome <- final_microbiome %>% select(-c("UNMAPPED","UNGROUPED"))
saveRDS(final_microbiome, here("../melonnpan_data/microbiome_training_data.RDS"))

filtered_metabolome_data <- filter_features(complete_metabolome_training_data, threshold = 0.000001, min_percentage = 0.10)
final_metabolome <- as.data.frame(t(filtered_metabolome_data))

# Perform variance filtering
#metabolome <- combined_relative_abundance %>%
  #rowwise() %>%
  #mutate(variance = var(c_across(everything()), na.rm = TRUE))
#summary(metabolome$variance)
#min_variance_threshold <- 2.1e-09
#metabolome <- metabolome %>%
  #filter(variance > min_variance_threshold) %>%
  #select(-variance) 

### Run Melonnpan Train ---
?melonnpan.train
# Remove columns (features) that contain any NA values
df_clean <- final_metabolome[, colSums(is.na(final_metabolome)) == 0]

melonnpan.train(metab = df_clean,
                metag = final_microbiome,
                output = here("Shotgun/melonnpan/IL10_TL1A_model/"),
                cores =4, 
                seed= 1234)






### If receiving error messages of "binary" "monomorphic" ---


# Step 1: Identify binary columns
binary_columns <- sapply(final_microbiome, function(x) length(unique(na.omit(x))) == 2)

# Print names of binary columns for reference
if (any(binary_columns)) {
  print("Binary columns identified:")
  print(names(microbiome)[binary_columns])
} else {
  print("No binary columns found.")
}


# Step 2: Remove binary columns from the metabolome dataset
microbiome_filtered <- microbiome %>%
  select(which(!binary_columns))  # Keep columns that are NOT binary

?melonnpan.train()
str(metabolome)
data_types <- sapply(metabolome, class)
# Check if all columns in the dataframe are numeric
are_all_numeric <- all(sapply(final_metabolome, is.numeric))

# Print the result
if (are_all_numeric) {
  print("All columns are numeric.")
} else {
  print("Not all columns are numeric.")
}

# Step 1: Identify monomorphic columns
monomorphic_columns <- sapply(microbiome_filtered, function(x) length(unique(na.omit(x))) == 1)

# Print names of monomorphic columns for reference
if (any(monomorphic_columns)) {
  print("Monomorphic columns identified:")
  print(names(microbiome_filtered)[monomorphic_columns])
} else {
  print("No monomorphic columns found.")
}

# Step 2: Remove monomorphic columns from the metabolome dataset
microbiome_final <- microbiome_filtered %>%
  select(which(!monomorphic_columns))  # Keep columns that are NOT monomorphic


melonnpan.train(metab = metabolome,
                metag = microbiome_final,
                output = here("Shotgun/melonnpan/"))

### Run Melonnpan using pretrained model ---
predict_metabolites <- function(df_input, weights_path, output_dir, train_metag, threshold = 0.01/100, sample_fraction = 0.10) {
  
  # Check if df_input is a file path or a dataframe
  if (is.character(df_input)) {
    # If it's a file path, read the file
    df <- read.delim(df_input, row.names = 1)
  } else if (is.data.frame(df_input)) {
    # If it's already a dataframe, use it directly
    df <- df_input
  } else {
    stop("df_input must be either a file path or a dataframe")
  }
  # Read in the normalized data and transpose it
  shotgun_dat <- as.data.frame(t(df))
  shotgun_dat <- shotgun_dat[,-1]
  
  # Read in the trained weights
  weights <- read.delim(weights_path)
  
  # Find overlapping IDs
  overlap <- intersect(names(shotgun_dat), weights$ID)
  
  # Filter data based on abundance threshold and prevalence across samples
  binary_abundance <- shotgun_dat > threshold
  proportion_samples <- colMeans(binary_abundance)
  filtered_data <- shotgun_dat[, proportion_samples >= sample_fraction]
  
  # Predict metabolites using MelonnPan
  metabolites <- melonnpan::melonnpan.predict(metag = filtered_data,
                                              output = output_dir,
                                              weight.matrix = weights_path,
                                              train.metag = train_metag)
  
  return(metabolites)
}


microbiome_train <- readRDS(here("../melonnpan_data/microbiome_training_data.RDS"))

test <- read.delim(here("Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),row.names=1)
shotgun_result <- predict_metabolites(
  df_path = here("Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),
  weights_path = here("Shotgun/melonnpan/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Shotgun/melonnpan/IL10_TL1A_model/"),
  train_metag = microbiome_train
)
shotgun_result$RTSI

microbiome <- read.delim(here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_counts.tsv"), row.names=1)
microbiome <- microbiome %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))
write.table(microbiome, here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_relabun.tsv"),sep = "\t")
test <- read.delim(here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_relabun.tsv"),row.names=1)

ucla_o_spf_result <- predict_metabolites(
  df_path = here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_relabun.tsv"),
  weights_path = here("Shotgun/melonnpan/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Regional-Mouse-Biogeography-Analysis/melonnpan/"),
  train_metag = microbiome_train
)

ucla_o_spf_result$RTSI

df<-read.delim(here("Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),row.names=1)
shotgun_dat <- as.data.frame(t(df))
shotgun_dat <- shotgun_dat[,-1]
names(shotgun_dat)

weights <- read.delim(here("Shotgun/melonnpan/MelonnPan_Trained_Weights.txt"))
weights$ID

overlap <- intersect(names(shotgun_dat), weights$ID)

# Mean relative abundance of 0.01% and prevalent in 10% of samples
threshold <- 0.01/100
sample_fraction <- 0.10
binary_abundance <- shotgun_dat > threshold
proportion_samples <- colMeans(binary_abundance)
filtered_data <- shotgun_dat[,proportion_samples>=sample_fraction]

metabolites <- melonnpan::melonnpan.predict(metag = filtered_data,
                            output = here("Shotgun/melonnpan/IL10_TL1A_model/"),
                            weight.matrix=here("Shotgun/melonnpan/IL10_TL1A_model/MelonnPan_Trained_Weights.txt"),
                            train.metag = final_microbiome)

rtsi_scores <- metabolites$RTSI

# Prevalent in 15 % of samples 
threshold <- 0.15
binary_presence <- shotgun_dat > 0
proportion_present <- colMeans(binary_presence)
filtered_dat2 <- shotgun_dat[,proportion_present>=threshold]

melonnpan::melonnpan.predict(metag = filtered_dat2,
                             output = here("Shotgun/melonnpan/"))

