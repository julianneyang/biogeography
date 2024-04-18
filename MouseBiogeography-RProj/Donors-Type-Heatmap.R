library(data.table)
library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)
library(gplots)
library(dplyr)
here::i_am("MouseBiogeography-RProj/Donors-Type-Heatmap.R")


target <- find_features_union_for_type_heatmap(duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Proximal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Distal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv")

#Query the target vector against all_results.tsv for each site
df <- query_type_features_union(
  target,
  duo_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Duodenum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  jej_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Jejunum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  ile_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Ileum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  cec_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Cecum-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  pc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Proximal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv",
  dc_filepath = "Donors-Analysis/differential_genera_type/L6-LumRef-CLR-Distal_Colon-ComBat-SeqRunSexType-1-MsID-DonorID/significant_results.tsv")

