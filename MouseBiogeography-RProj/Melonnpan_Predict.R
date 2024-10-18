library(melonnpan)
library(here)
library(stringr)
library(Maaslin2)

here::i_am("MouseBiogeography-RProj/Melonnpan_Predict.R")

microbiome_train <- readRDS(here("../melonnpan_data/microbiome_training_data.RDS"))

shotgun_result <- predict_metabolites(
  df_input = here("Shotgun/relab_normalized/merged_humann_genefamilies_relab_normalized_unstratified_ko.tsv"),
  weights_path = here("Shotgun/melonnpan/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Shotgun/melonnpan/IL10_TL1A_model/"),
  train_metag = microbiome_train
)
shotgun_result$RTSI

### UCLA O. SPF 




## UCLA O. SPF --
microbiome <- read.delim(here("Regional-Mouse-Biogeography-Analysis/picrust_output/UCLA_O_SPF_KO_counts.tsv"), row.names=1)
microbiome <- microbiome %>% 
  mutate(across(everything(), ~ . / sum(., na.rm = TRUE)))

input_metadata <-read.delim(here("Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/starting_files/Regional-Combat-Metadata.tsv"),header=TRUE, row.names=1) #mapping file

target <- colnames(microbiome)
input_metadata = input_metadata[match(target, row.names(input_metadata)),]
target == row.names(input_metadata)

df_input_metadata <- as.data.frame(input_metadata)
df_input_metadata$Sequencing_Run <- factor(df_input_metadata$Sequencing_Run)
df_input_metadata$Line <- factor(df_input_metadata$Line)
df_input_metadata$MouseID_Line <- factor(df_input_metadata$MouseID_Line)
df_input_metadata$Sex <- factor(df_input_metadata$Sex)
df_input_metadata$Type <- factor(df_input_metadata$Type, levels=c("Luminal", "Mucosal"))
df_input_metadata$SampleID <- row.names(df_input_metadata)
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))

samples <- df_input_metadata %>%
  filter(SampleID %in% names(microbiome)) %>%
  filter(Type=="Luminal") %>%
  pull(SampleID)

df_input_data <- microbiome[,samples]

#df_input_data <- filter_features(df_input_data)

# Predict Metabolite Compostion - 
ucla_o_spf_result <- predict_metabolites(
  df_input = df_input_data,
  weights_path = here("Shotgun/melonnpan/MelonnPan_Trained_Weights.txt"),
  output_dir = here("Regional-Mouse-Biogeography-Analysis/melonnpan/"),
  train_metag = microbiome_train
)

lum_mets <- ucla_o_spf_result$pred
lum_mets <- read.delim(here("Regional-Mouse-Biogeography-Analysis/melonnpan/MelonnPan_Predicted_Metabolites.txt"))
lum_mets <- lum_mets %>% 
  column_to_rownames("ID")
lum_mets <- as.data.frame(t(lum_mets))

# Run Maaslin2 

lum_meta <- df_input_metadata %>%
  filter(SampleID %in% names(microbiome)) %>%
  filter(Type=="Luminal") 

site_fe <- c("Sequencing_Run", "Line","Sex", "Type")
ranef <- c("MouseID_Line")
refs <- c("Sequencing_Run,Hiseq_April_Nineteen","Line,ItgCre","Site,Distal_Colon")

test_normality_metabolites <- function(df) {
  # Subset the first 50 metabolites (rows)
  metabolites_subset <- df[1:50, ]
  
  # Apply Shapiro-Wilk test to each row (metabolite) and collect p-values
  normality_results <- apply(metabolites_subset, 1, function(metabolite) {
    shapiro.test(metabolite)$p.value
  })
  
  # Create a data frame of p-values
  results_df <- data.frame(
    Metabolite = rownames(metabolites_subset),
    Shapiro_p_value = normality_results
  )
  
  return(results_df)
}
normality <- test_normality_metabolites(lum_mets)

plot_metabolite_density_patchwork <- function(df) {
  # Subset the first 50 metabolites (rows)
  metabolites_subset <- df[1:50, ]
  
  # List to store the plots
  plot_list <- list()
  
  # Loop through each metabolite (row)
  for (i in 1:5) {
    metabolite_data <- as.numeric(metabolites_subset[i, ])
    metabolite_name <- rownames(metabolites_subset)[i]
    
    # Create the density plot
    p <- ggplot(data.frame(x = metabolite_data), aes(x = x)) +
      geom_density(fill = "blue", alpha = 0.5) +
      labs(title = metabolite_name, x = "Relative abundance", y = "Density") +
      theme_minimal()
    
    # Add the plot to the list
    plot_list[[i]] <- p
  }
  
  # Use patchwork to combine the plots into a grid
  grid <- wrap_plots(plot_list, ncol = 5)
  
  return(grid)
}

# Usage example:
plot_metabolite_density_patchwork(lum_mets)


fit_data = Maaslin2(input_data=lum_mets, 
                    input_metadata=lum_meta, 
                    output = here("Regional-Mouse-Biogeography-Analysis/melonnpan/Luminal-SeqRun-Line-Sex-Site_General_1-MouseID"), 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), 
                    random_effects = c("MouseID_Line"),
                    reference =  c("Sequencing_Run,Hiseq_April_Nineteen","Line,ItgCre","Site,Distal_Colon"),
                    min_prevalence=0.15,
                    #analysis_method = "LM",
                    normalization="none", transform ="log",plot_heatmap = FALSE,plot_scatter = FALSE)

res <- read.delim(here("Regional-Mouse-Biogeography-Analysis/melonnpan/Luminal-SeqRun-Line-Sex-Site_General_1-MouseID/all_results.tsv"))
res_signif <- res %>% 
  filter(metadata=="Site_General") %>% 
  filter(qval < 0.05) %>% 
  filter(abs(coef) > 1)
res_signif$feature
ggplot(res_signif, aes(x = coef, y = feature, fill = value)) +
  geom_bar(stat = 'identity', position = 'dodge') +  # Create horizontal bars
  #scale_fill_manual(values = c('SI' = 'blue', 'colon' = 'orange')) +  # Custom colors
  labs(x = 'Coefficient', y = 'Feature', title = 'UCLA O. SPF Luminal Metabolites') +
  theme_cowplot(16)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Lum-ComBat-SeqRunLineSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID_Line"),normalization="clr",
                    min_prevalence=0.15,
                    transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

#Mucosal
df_input_metadata$Site_General <- factor(df_input_metadata$Site_General, levels=c("Colon","SI"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-ColonRef-CLR-Muc-ComBat-SeqRunLineSexSite_General-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site_General"), 
                    random_effects = c("MouseID_Line"),
                    min_prevalence=0.15,
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)
df_input_metadata$Site <- factor(df_input_metadata$Site, levels=c("Distal_Colon", "Proximal_Colon", "Cecum", "Ileum","Jejunum","Duodenum"))
fit_data = Maaslin2(input_data=df_input_data, input_metadata=df_input_metadata, output = "Maaslin2_L2/UCLA_O_SPF/L2-DCvsAll-CLR-Muc-ComBat-SeqRunLineSexSite-1-MsID", 
                    fixed_effects = c("Sequencing_Run","Line","Sex", "Site"), 
                    random_effects = c("MouseID_Line"),
                    min_prevalence=0.15,
                    normalization="clr", transform ="none",plot_heatmap = FALSE,plot_scatter = FALSE)

### Maaslin2 ---
