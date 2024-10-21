

devtools::install_github("borenstein-lab/fishtaco-plot")
library(FishTacoPlot)
library(here)
library(tidyverse)
library(dplyr)
library(cowplot)

fp <- "Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/"
taxonomy <- read.delim(here(paste0(fp,"starting_files/UCLA_taxonomy.tsv")))
taxonomy <- taxonomy %>%
  separate(Taxon, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "; ") #%>%
  #mutate(across(everything(), ~ gsub(" ", "", .)))  # Remove spaces
taxonomy <- taxonomy %>% 
  mutate(order=seq(1,323,1))%>%
  mutate(Species=paste0(Species,order))%>%
  select(-c(order)) 
write.table(x=taxonomy,file=here(paste0(fp,"../fish_taco/taxonomy.tsv")),
            col.names = TRUE,row.names=FALSE,sep = "\t")

# Control was Duodenum, Case was Distal Colon
# Depleted in Distal Colon = Enriched in Duodenum
 p <-MultiFunctionTaxaContributionPlots(input_dir=here(paste0(fp,"../fish_taco/muc_duo_dc_0_1")), 
                                        input_prefix="fishtaco_out",
                                        input_taxa_taxonomy=here(paste0(fp,"../fish_taco/taxonomy.tsv")),
                                        input_permutation="single_taxa",
                                        plot_type="bars",
                                        sort_by="list",
                                        add_predicted_da_markers = FALSE,
                                        input_function_filter_list=c("K01667"),
                                        show_only_enriched_taxa = FALSE,
                                        show_only_enriched_functions = TRUE,
                                        
                                    #flip_case_control = TRUE
                                      )
?MultiFunctionTaxaContributionPlots

 big_color_vector <- readRDS(here("global_genera_cols.RDS"))
 names(big_color_vector)<-NULL
 
duo <- p + scale_fill_manual(values=big_color_vector) +
    theme_cowplot(12) +
    theme(legend.position="right")+
    ggtitle("KO depleted in Duo (left) vs. DC (right)")
dev.new(width=10,height=10)
duo
grid::ge

  ?SingleFunctionTaxaContributionPlots
?MultiFunctionTaxaContributionPlots()
