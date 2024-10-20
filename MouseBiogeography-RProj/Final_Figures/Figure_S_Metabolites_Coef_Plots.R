###### The Big and Tragic Rearrangement of Mouse Biogeography ---
### Date: 3.10.2023
### Figure Number: 6
### Figure Contents: GMM Site Heatmaps for all Cohorts 
###### whining ends here ---

library(ggplot2)
library(tidyverse)
#library(rlang)
library(cowplot)
library(viridis)
#library(plyr)
library(gridExtra)
library(paletteer)
library(ComplexUpset)
library(here)

#Replace with filepath to package Microbiome.Biogeography
setwd("/home/julianne/Documents/microbiome.biogeography/")
devtools::document()
library("Microbiome.Biogeography")
setwd("/home/julianne/Documents/biogeography/")

here::i_am("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Coef_Plots.R")

### Upset Plot ---

file_paths <- c("Regional-Mouse-Biogeography-Analysis/melonnpan/Mets-SeqRunLineSexSite_General_1-MouseID/all_results.tsv",
                "CS_SPF/melonnpan/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv",
                "Donors-Analysis/melonnpan/Mets-SeqRunSexSite_General-1-MsID-DonorID/all_results.tsv",
                "Humanized-Biogeography-Analysis/melonnpan/HUM_SD_Gavage/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv",
                "Humanized-Biogeography-Analysis/melonnpan/SPF_Gavage/Mets-SeqRunSexSite_General-1-MsID/all_results.tsv")


cohort_prefixes <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_MD_Gavage",
                     "HUM_SD_Gavage",
                     "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site_General")


id_features <- all_taxa %>% mutate(coef_dir = ifelse(coef > 0, "POS", "NEG"))
id_features <- id_features%>% select(c("feature","Cohort","coef_dir")) %>% unique()
id_f_long <- id_features %>% 
  mutate(value = 1)
id_df_wide <- id_f_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)

id_df_wide <- as.data.frame(id_df_wide)
#id_df_wide <- id_df_wide %>% mutate(SPF_Gavage = 0)

all_taxa <- all_taxa %>% select(c("feature", "Cohort")) %>% unique()
df_long <- all_taxa %>% 
  mutate(value = 1)

df_wide <- df_long %>%
  pivot_wider(names_from = Cohort, values_from = value, values_fill = 0)
df_wide <- as.data.frame(df_wide)
#df_wide <- df_wide %>% mutate(SPF_Gavage = 0)

all_datasets <- cohort_prefixes
all_datasets <- names(df_wide)[-1]
mets_upset <-ComplexUpset::upset(df_wide, all_datasets,width_ratio=0.15, name="",
                                 base_annotations=list(
                                   'Intersection size'=intersection_size(counts=TRUE,mapping=aes(fill='bars_color')) + 
                                     scale_fill_manual(values=c('bars_color'='skyblue'), guide='none')),
                                 themes=list(
                                   default=theme(
                                     axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                   ),
                                   intersections_matrix=theme(
                                     axis.ticks.x=element_blank(),
                                     axis.text.x=element_blank(),
                                   ),
                                   set_sizes=(ylab('Region-specific genera'))
                                 ))

                     #min_degree=3,
                   #themes=upset_default_themes(text=element_text(size=16), colour="black"))
                   #themes=upset_modify_themes(
                     #list(
                       #'intersections_matrix'=theme(text=element_text(size=16,colour = "black"))
                     #)
                   #))
mets_upset

id_df_wide$count_ones <- rowSums(id_df_wide %>% select(all_of(cohort_prefixes)))
df_filtered <- id_df_wide[id_df_wide$count_ones >= 5, ]
length(df_filtered$feature)
length(unique(df_filtered$feature))
df_filtered <- df_filtered[, -which(names(df_filtered) == "count_ones")]
metabs_of_interest <- df_filtered$feature
write_rds(metabs_of_interest,here("Highlighted_metabolites.RDS"))

spf <- c("UCLA_O_SPF", "CS_SPF")
spf <- id_df_wide %>%
  mutate(SPF_only = case_when(
    UCLA_O_SPF == 1 & CS_SPF == 1 &
      HUM_MD_Gavage == 0 & HUM_SD_Gavage == 0 & SPF_Gavage == 0 ~ "SPF",
    TRUE ~ NA_character_  # Assign NA if conditions are not met
  ))

gavage <- c("HUM_MD_Gavage", "HUM_SD_Gavage",
            "SPF_Gavage")
gavage <- id_df_wide %>%
  mutate(gavage_only = case_when(
    UCLA_O_SPF == 0 & CS_SPF == 0 &
      rowSums(select(., HUM_MD_Gavage, HUM_SD_Gavage, SPF_Gavage) == 1) >= 2 ~ "Gavage",
    TRUE ~ NA_character_
  ))

gavage_filtered <- gavage %>% 
  filter(gavage_only == "Gavage")

spf_filtered <- spf %>% 
  filter(SPF_only == "SPF")

length(spf_filtered$feature)
length(unique(spf_filtered$feature))

unique_spf <- spf_filtered$feature
unique_gavage <- gavage_filtered$feature
write_rds(unique_spf,here("SPF_only_metabolites.RDS"))
write_rds(unique_gavage, here("Gavage_only_metabolites.RDS"))

### Make a big barplot for shared_metabs ---
gmm <- readRDS(here("Highlighted_GMM_Fig_6.RDS"))
names(gmm)
gbm <- readRDS(here("Highlighted_Luminal_GBM.RDS"))
names(gbm)
mets <- readRDS(here("Highlighted_metabolites.RDS"))
names(mets) <-  c(
  "Lipids",                 # "LPC.O.16.0"
  "Lipids",                 # "D.erythro.N.stearoylsphingosine"
  "Small Organic Compounds",# "phosphate"
  "Amino Acids",            # "Thr.Ala"
  "Lipids",                 # "Cer.38.3.O4.Cer.24.2.O3.14.1.2OH."
  "Lipids",                 # "phosphoinositol"
  "Lipids",                 # "Cer.36.1.O3.Cer.18.0.O2.18.1"
  "Amino Acids",            # "Ala.Gln"
  "Lipids",                 # "Cer.36.3.O2.Cer.18.1.O2.18.2"
  "Amino Acids",            # "Phe.Thr"
  "Lipids",                 # "Cer.17.1.2O.24.1"
  "Amino Acids",            # "Ala.Ala.Lys"
  "Amino Acids",            # "Ala.Thr"
  "Small Organic Compounds",# "Cholesterol"
  "Amino Acids",            # "Leu.Trp"
  "Amino Acids",            # "Val.Ser"
  "Amino Acids",            # "Ile.Ala"
  "Lipids",                 # "X1..1Z.Hexadecenyl..sn.glycero.3.phosphocholine"
  "Amino Acids",            # "Phe.Val"
  "Amino Acids",            # "Ile.Gln"
  "Small Organic Compounds",# "X3..1.Pyrazolyl..alanine"
  "Small Organic Compounds",# "X3.Methyl_gamma.butyrolactone"
  "Small Organic Compounds",# "piperidone"
  "Amino Acids",            # "Ile.Val.Arg"
  "Small Organic Compounds",# "X2.Hydroxyphenethyl.alcohol"
  "Lipids",                 # "Cer.d44.2"
  "Lipids",                 # "Cer.42.2.O2.Cer.18.1.O2.24.1"
  "Lipids",                 # "Cer.38.2.O2.Cer.18.1.O2.20.1"
  "Lipids",                 # "Cer.38.5.O3.Cer.18.1.O2.20.4.O"
  "Lipids",                 # "SPB.18.1.O2"
  "Small Organic Compounds",# "X5.aminovaleric.acid"
  "Small Organic Compounds",# "glyceric.acid"
  "Lipids",                 # "Cer.d40.2"
  "Small Organic Compounds",# "lactulose"
  "Small Organic Compounds",# "glycerol.3.galactoside"
  "Lipids",                 # "HexCer.42.2.O4.HexCer.18.1.O3.24.1.2OH."
  "Amino Acids",            # "Isoleucine"
  "Amino Acids",            # "Tryptophan"
  "Small Organic Compounds",# "Deoxycholate"
  "Small Organic Compounds",# "guanidinosuccinate"
  "Lipids",                 # "Cer.NDS.d36.0"
  "Lipids",                 # "Cer.36.1.O2.Cer.18.0.O2.18.1"
  "Lipids"                  # "Cer.38.0.O3.Cer.20.0.O2.18.0.O"
)

spf_mets <- readRDS(here("SPF_only_metabolites.RDS"))
names(spf_mets) <- c(
    "Small Organic Compounds",   # "Isorhamnetin.3.O.rutinoside"    
    "Amino Acids",               # "Glu.Ala.Lys"
    "Amino Acids",               # "Lys.Ser"
    "Amino Acids",               # "Gln.Leu.Lys"
    "Amino Acids",               # "Tyr.Asn"
    "Small Organic Compounds",   # "SPERMIDINE"
    "Small Organic Compounds",   # "Catechol"
    "Lipids",                    # "Taurocholic.acid"
    "Lipids",                    # "SM.29.2.2O.8.0"
    "Amino Acids",               # "Met.Met"
    "Small Organic Compounds",   # "X3..3.Hydroxyphenyl.propionic.acid"
    "Amino Acids",               # "Thr.Asp"
    "Small Organic Compounds",   # "X2.Amino.3.methoxybenzoic.acid"
    "Small Organic Compounds",   # "N8.acetylspermidine"
    "Small Organic Compounds",   # "X1.Methylnicotinamide"
    "Amino Acids",               # "beta.glutamic.acid"
    "Small Organic Compounds",   # "adenosine.1"
    "Amino Acids",               # "N2.acetyl.histidine"
    "Lipids",                    # "PC.42.5"
    "Small Organic Compounds",   # "Mandelic.acid"
    "Lipids",                    # "SHexCer.34.1.3O"
    "Lipids",                    # "SM.d42.1"
    "Small Organic Compounds",   # "E.N.Deoxyfructosyllysine"
    "Lipids",                    # "PC.38.4.Isomer.B"
    "Lipids",                    # "NAGlySer.34.1.O2.NAGlySer.17.0.O.FA.17.0"
    "Lipids",                    # "PE.34.1.1"
    "Small Organic Compounds",   # "dihydro.3.coumaric.acid"
    "Lipids",                    # "HexCer.34.1.O2.HexCer.18.1.O2.16.0"
    "Amino Acids",               # "Tyr.Thr"
    "Lipids",                    # "MGDG.17.0_18.2"
    "Small Organic Compounds",   # "maleic.acid"
    "Small Organic Compounds",   # "X1.2.Cyclohexanedione"
    "Lipids",                    # "PG.O.36.2.PG.O.18.1_18.1"
    "Small Organic Compounds",   # "lactic.acid"
    "Lipids",                    # "Cer.18.2.2O.22.0"
    "Lipids",                    # "Cer.42.0.O3.Cer.18.0.O2.24.0.O"
    "Small Organic Compounds",   # "Methyl.beta.galactopyranoside"
    "Lipids",                    # "lauric.acid"
    "Lipids",                    # "FA.41.1"
    "Small Organic Compounds",   # "dibenzylamine"
    "Lipids",                    # "FA.37.1"
    "Lipids",                    # "SM.d33.1"
    "Amino Acids",               # "beta.homoleucine"
    "Lipids",                    # "TG.53.0"
    "Amino Acids",               # "Asn.Tyr"
    "Vitamins and Cofactors",    # "tocopherol.gamma."
    "Lipids",                    # "Nervonic.acid"
    "Lipids",                    # "Cer.d40.1"
    "Lipids",                    # "FA.39.1"
    "Lipids",                    # "FA.33.1"
    "Lipids",                    # "Cer.18.0.2O_20.0"
    "Small Organic Compounds",   # "Caffeic.acid"
    "Small Organic Compounds",   # "X13.KODE"
    "Lipids",                    # "myristic.acid"
    "Lipids",                    # "Cer.d42.1"
    "Amino Acids",               # "N.N.epsiloN.dimethyl.lysine"
    "Lipids",                    # "LPC.18.1.Isomer.B"
    "Lipids",                    # "CAR.16.0"
    "Lipids"                     # "TG.58.6"
  )

spf_mets <- readRDS(here("SPF_only_metabolites.RDS"))
names(spf_mets) <- c(
  "Small Organic Compounds",   # "Isorhamnetin.3.O.rutinoside"    
  "Amino Acids",               # "Glu.Ala.Lys"
  "Amino Acids",               # "Lys.Ser"
  "Amino Acids",               # "Gln.Leu.Lys"
  "Amino Acids",               # "Tyr.Asn"
  "Small Organic Compounds",   # "SPERMIDINE"
  "Small Organic Compounds",   # "Catechol"
  "Lipids",                    # "Taurocholic.acid"
  "Lipids",                    # "SM.29.2.2O.8.0"
  "Amino Acids",               # "Met.Met"
  "Small Organic Compounds",   # "X3..3.Hydroxyphenyl.propionic.acid"
  "Amino Acids",               # "Thr.Asp"
  "Small Organic Compounds",   # "X2.Amino.3.methoxybenzoic.acid"
  "Small Organic Compounds",   # "N8.acetylspermidine"
  "Small Organic Compounds",   # "X1.Methylnicotinamide"
  "Amino Acids",               # "beta.glutamic.acid"
  "Small Organic Compounds",   # "adenosine.1"
  "Amino Acids",               # "N2.acetyl.histidine"
  "Lipids",                    # "PC.42.5"
  "Small Organic Compounds",   # "Mandelic.acid"
  "Lipids",                    # "SHexCer.34.1.3O"
  "Lipids",                    # "SM.d42.1"
  "Small Organic Compounds",   # "E.N.Deoxyfructosyllysine"
  "Lipids",                    # "PC.38.4.Isomer.B"
  "Lipids",                    # "NAGlySer.34.1.O2.NAGlySer.17.0.O.FA.17.0"
  "Lipids",                    # "PE.34.1.1"
  "Small Organic Compounds",   # "dihydro.3.coumaric.acid"
  "Lipids",                    # "HexCer.34.1.O2.HexCer.18.1.O2.16.0"
  "Amino Acids",               # "Tyr.Thr"
  "Lipids",                    # "MGDG.17.0_18.2"
  "Small Organic Compounds",   # "maleic.acid"
  "Small Organic Compounds",   # "X1.2.Cyclohexanedione"
  "Lipids",                    # "PG.O.36.2.PG.O.18.1_18.1"
  "Small Organic Compounds",   # "lactic.acid"
  "Lipids",                    # "Cer.18.2.2O.22.0"
  "Lipids",                    # "Cer.42.0.O3.Cer.18.0.O2.24.0.O"
  "Small Organic Compounds",   # "Methyl.beta.galactopyranoside"
  "Lipids",                    # "lauric.acid"
  "Lipids",                    # "FA.41.1"
  "Small Organic Compounds",   # "dibenzylamine"
  "Lipids",                    # "FA.37.1"
  "Lipids",                    # "SM.d33.1"
  "Amino Acids",               # "beta.homoleucine"
  "Lipids",                    # "TG.53.0"
  "Amino Acids",               # "Asn.Tyr"
  "Small Organic Compounds",    # "tocopherol.gamma."
  "Lipids",                    # "Nervonic.acid"
  "Lipids",                    # "Cer.d40.1"
  "Lipids",                    # "FA.39.1"
  "Lipids",                    # "FA.33.1"
  "Lipids",                    # "Cer.18.0.2O_20.0"
  "Small Organic Compounds",   # "Caffeic.acid"
  "Small Organic Compounds",   # "X13.KODE"
  "Lipids",                    # "myristic.acid"
  "Lipids",                    # "Cer.d42.1"
  "Amino Acids",               # "N.N.epsiloN.dimethyl.lysine"
  "Lipids",                    # "LPC.18.1.Isomer.B"
  "Lipids",                    # "CAR.16.0"
  "Lipids"                     # "TG.58.6"
)

gavage_mets <- readRDS(here("Gavage_only_metabolites.RDS"))
names(gavage_mets) <-c(
  "Lipid",                        # 1 decanedioic.acid
  "Lipid",                        # 2 TG.43.1
  "Amino Acid",                  # 3 Tryptoline
  "Amino Acid",                  # 4 carnitine
  "Lipid",                        # 5 campesterol
  "Small Organic Compound",       # 6 deprenyl
  "Lipid",                        # 7 X.Z..5.8.11.trihydroxyoctadec.9.enoic.acid
  "Lipid",                        # 8 SM.22.1.2O.12.0
  "Lipid",                        # 9 FA.24.0
  "Small Organic Compound",       # 10 X5.Aminosalicylic.acid
  "Small Organic Compound",       # 11 X6.Aminoindazole
  "Small Organic Compound",       # 12 X13.KODE
  "Lipid",                        # 13 FA.25.0
  "Lipid",                        # 14 Cer.NDS.d41.0
  "Small Organic Compound",       # 15 NAGlySer.31.1.O2.NAGlySer.16.0.O.FA.15.0.
  "Small Organic Compound",       # 16 PGE2
  "Small Organic Compound",       # 17 X8.Oxo.2.deoxyadenosine
  "Lipid",                        # 18 FA.27.0
  "Small Organic Compound",       # 19 X7.Methylguanine
  "Lipid",                        # 20 FA.26.0
  "Lipid",                        # 21 MGDG.17.0_18.2
  "Small Organic Compound",       # 22 X2..2.6.dimethylmorpholin.4.yl.ethanol
  "Lipid",                        # 23 FA.22.0
  "Small Organic Compound",       # 24 Isorhamnetin.3.O.rutinoside
  "Lipid",                        # 25 X4.HDoHE
  "Small Organic Compound",       # 26 Biliverdin
  "Small Organic Compound",       # 27 X5.Acetylamino.6.amino.3.methyluracil
  "Lipid",                        # 28 PG.33.0.PG.16.0_17.0
  "Small Organic Compound",       # 29 X6.deoxyglucose
  "Small Organic Compound",       # 30 Caffeic.acid
  "Lipid",                        # 31 SM.d34.1
  "Lipid",                        # 32 lauric.acid
  "Small Organic Compound",       # 33 arabitol
  "Lipid",                        # 34 SM.d36.0.Isomer.B
  "Lipid",                        # 35 CAR.16.0
  "Small Organic Compound",       # 36 lyxitol
  "Lipid",                        # 37 PC.34.0
  "Lipid",                        # 38 FA.24.2
  "Lipid",                        # 39 SM.24.1
  "Lipid",                        # 40 TG.61.3.TG.18.0_25.0_18.3
  "Lipid",                        # 41 CAR.18.0.Isomer.A
  "Lipid",                        # 42 PC.38.4.Isomer.B
  "Small Organic Compound",       # 43 NAGlySer.30.1.O2.NAGlySer.15.0.O.FA.15.0.
  "Lipid",                        # 44 FA.23.0
  "Lipid",                        # 45 FA.28.0
  "Small Organic Compound",       # 46 X2..deoxyguanosine
  "Amino Acid",                  # 47 N.Acetylserine
  "Lipid",                        # 48 Cer.d39.1
  "Lipid",                        # 49 TG.54.5.O2.TG.18.1_18.2_18.2.O2
  "Lipid",                        # 50 Nervonic.acid
  "Lipid",                        # 51 LPE.22.6
  "Small Organic Compound",       # 52 N.N..Diacetylcystine
  "Lipid",                        # 53 LPG.16.0
  "Lipid",                        # 54 PE.34.0
  "Lipid",                        # 55 Cer.43.0.O2.Cer.17.0.O2.26.0
  "Lipid",                        # 56 Cer.d38.1
  "Lipid",                        # 57 PC.36.4.Isomer.B
  "Lipid",                        # 58 Cer.35.0.O3.Cer.18.0.O2.17.0.O
  "Lipid",                        # 59 TG.34.0
  "Amino Acid",                  # 60 tyramine
  "Small Organic Compound",       # 61 X2..Deoxyinosine
  "Small Organic Compound",       # 62 X4.Pentylaniline
  "Small Organic Compound",       # 63 Histamine
  "Lipid",                        # 64 FA.29.1
  "Lipid",                        # 65 FA.27.1
  "Lipid",                        # 66 FA.33.1
  "Small Organic Compound",       # 67 inosine
  "Amino Acid",                  # 68 Met.His
  "Lipid",                        # 69 TG.45.0
  "Lipid",                        # 70 Phosphatidylethanolamine.alkenyl.16.0.20.4
  "Amino Acid",                  # 71 delta.Hydroxylysine
  "Small Organic Compound",       # 72 Tomatidine
  "Small Organic Compound",       # 73 X2.deoxyuridine
  "Small Organic Compound",       # 74 dihydro.3.coumaric.acid
  "Lipid",                        # 75 FA.31.1
  "Lipid",                        # 76 Cer.d34.2.Isomer.A
  "Lipid",                        # 77 PE.36.3
  "Small Organic Compound",       # 78 X2.Amino.3.methoxybenzoic.acid
  "Small Organic Compound",       # 79 p.coumaric.acid
  "Small Organic Compound",       # 80 X2..n.ethyl.n.m.toluidino.ethanol
  "Small Organic Compound",       # 81 X3..3.Hydroxyphenyl.propionic.acid
  "Small Organic Compound",       # 82 Mandelic.acid
  "Lipid",                        # 83 FA.32.0
  "Small Organic Compound",       # 84 X1.Methylnicotinamide
  "Lipid",                        # 85 FA.25.1
  "Amino Acid",                  # 86 N.fructosyl.isoleucine
  "Small Organic Compound",       # 87 cadaverine
  "Amino Acid",                  # 88 Val.Val.Lys
  "Lipid",                        # 89 SM.d40.0.Isomer.B
  "Lipid",                        # 90 Cer.18.0.2O_20.0
  "Lipid",                        # 91 FA.35.1
  "Lipid",                        # 92 SHexCer.34.1.3O
  "Lipid",                        # 93 FA.34.0
  "Amino Acid",                  # 94 Gln.Tyr
  "Lipid",                        # 95 GlcCer.d40.1
  "Lipid",                        # 96 FA.37.1
  "Lipid",                        # 97 Cholic.acid
  "Lipid",                        # 98 PE.32.1.PE.16.0_16.1
  "Amino Acid",                  # 99 Arg.Trp
  "Small Organic Compound",       # 100 FA.20.5..eicosapentaenoic.acid.
  "Small Organic Compound",       # 101 Choline
  "Small Organic Compound",       # 102 guanosine
  "Lipid",                        # 103 Taurocholic.acid
  "Amino Acid",                  # 104 Ile.Glu
  "Lipid",                        # 105 PE.O.16.1_22.6
  "Small Organic Compound",       # 106 X3..3.hydroxyphenyl.propionic.acid
  "Amino Acid",                  # 107 Tyr.Thr
  "Amino Acid",                  # 108 Gln.Pro
  "Lipid",                        # 109 FA.22.4
  "Lipid",                        # 110 LPC.O.22.0
  "Lipid",                        # 111 PI.36.2
  "Small Organic Compound",       # 112 X2.Salicylic.acid
  "Small Organic Compound",       # 113 3..phenylalanine
  "Small Organic Compound",       # 114 4..Cyanocobalamin
  "Lipid",                        # 115 PI.34.2
  "Small Organic Compound",       # 116 2.3..3..beta..carotene
  "Amino Acid",                  # 117 Histamine.5..
  "Small Organic Compound",       # 118 Histidine
  "Lipid",                        # 119 PG.29.0.PG.18.1_12.1
  "Amino Acid",                  # 120 His.Lys
  "Lipid",                        # 121 PG.33.0.PG.18.0_16.0
  "Lipid",                        # 122 Sphingomyelin.d38.2
  "Small Organic Compound",       # 123 N.Acetylspermidine
  "Lipid",                        # 124 PE.36.0
  "Amino Acid",                  # 125 Glu.His
  "Amino Acid",                  # 126 Val.Ala
  "Lipid",                        # 127 Acylcarnitine.d35.1
  "Lipid",                        # 128 PC.37.3
  "Lipid",                        # 129 Cer.25.0
  "Amino Acid"                   # 130 Gly.Gln
)

### Do this for Metabolites concordant across all 5 ---
# Load all the tsv files, add Cohort column, and merge them
all_data <- lapply(1:length(file_paths), function(i) {
  df <- read_tsv(file_paths[i])
  df <- df %>% mutate(Cohort = cohort_prefixes[i])
  return(df)
}) %>% bind_rows()


# Filter the dataframe based on values in "feature" column matching `mets`
filtered_data <- all_data %>% filter(feature %in% mets)%>% 
  filter(metadata == "Site_General") %>%
  filter(qval<0.05)

# Create new column "Site" based on "coef"
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "Colon", "SI"))
filtered_data <- filtered_data %>% 
  mutate(data_type = "16S")
# Create horizontal barplot with faceting by "Cohort"
cols = viridis(2)
metab_plot <- ggplot(filtered_data, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Region-Specific Metabolites") +
  theme_cowplot(16)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))


dev.new(width=10, height=20)
metab_plot

### add shotgun data to big plot ---

shotgun_fp <- c("Shotgun/melonnpan/UCLA_O_SPF/Mets-LineSexSite_General-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/CS_SPF/Mets-SexSite-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/HUM_SD_Gavage/Mets-SexSite-1-MsID/all_results.tsv",
                              "Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/all_results.tsv")

shotgun_prefix <- c("UCLA_O_SPF",
                     "CS_SPF",
                     "HUM_SD_Gavage",
                     "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site")

# Given 0 significant metabolites pick out the features mapping to mets
all_shotgun_data <- lapply(1:length(shotgun_fp), function(i) {
  df <- read_tsv(shotgun_fp[i])
  df <- df %>% mutate(Cohort = shotgun_prefix[i])
  return(df)
}) %>% bind_rows()

# Filter the dataframe based on values in "feature" column matching `mets`
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% mets)%>% 
  filter(metadata == "Site") 
filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Distal_Colon", "Jejunum"))
filtered_shotgun_data <-filtered_shotgun_data %>% 
  mutate(data_type = "Shotgun")
# Bind the dataframe together
full_df <- rbind(filtered_shotgun_data, filtered_data)
full_df$annotation <- names(mets)[match(full_df$feature, mets)]
full_df$feature <- gsub("X","", full_df$feature)
full_df$feature <- gsub("\\.","-", full_df$feature)
full_df$Cohort <- factor(full_df$Cohort, 
                         levels=c("UCLA_O_SPF","CS_SPF",
                                  "SPF_Gavage",
                                  "HUM_SD_Gavage","HUM_MD_Gavage"))

aa <- full_df %>% filter(annotation=="Amino Acids")
lipids <- full_df %>% filter(annotation=="Lipids")
org <- full_df %>% filter(annotation=="Small Organic Compounds")
cols <- viridis::viridis(4)
names(cols) <- c("Colon", "Distal_Colon",
                 "Jejunum", "SI")
aa_plot <- ggplot(aa, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Amino Acid Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

lipids_plot <- ggplot(lipids, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Lipid Metabolites") +
  theme_cowplot(12)+
  #theme(strip.text = element_text(
    #size = 8))+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

soc <- ggplot(org, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Small Organic Compounds") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

dev.new(width=10,height=10)

bottom <- plot_grid(lipids_plot,
                 soc, labels=c("C","D"),
                 rel_widths = c(1,0.85),
                 label_size = 20)
top <- plot_grid(mets_upset,
                 aa_plot, labels=c("A","B"), 
                 rel_widths = c(1,0.85),
                 label_size = 20)

dev.new(width=10,height=10)
top

write_rds(top, here("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Top.RDS"))
write_rds(bottom, here("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Bottom.RDS"))

#########################################

### Do this for SPF metabolites only --- 

#########################################

all_data <- lapply(1:length(file_paths), function(i) {
  df <- read_tsv(file_paths[i])
  df <- df %>% mutate(Cohort = cohort_prefixes[i])
  return(df)
}) %>% bind_rows()


# Filter the dataframe based on values in "feature" column matching `mets`
filtered_data <- all_data %>% filter(feature %in% spf_mets)%>% 
  filter(metadata == "Site_General") %>% 
  filter(qval < 0.05) %>%
  filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")

# Create new column "Site" based on "coef"
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "Colon", "SI"))
filtered_data <- filtered_data %>% 
  mutate(data_type = "16S")

### add shotgun data to big plot ---

shotgun_fp <- c("Shotgun/melonnpan/UCLA_O_SPF/Mets-LineSexSite_General-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/CS_SPF/Mets-SexSite-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/HUM_SD_Gavage/Mets-SexSite-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/all_results.tsv")

shotgun_prefix <- c("UCLA_O_SPF",
                    "CS_SPF",
                    "HUM_SD_Gavage",
                    "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site")

# Given 0 significant metabolites pick out the features mapping to mets
all_shotgun_data <- lapply(1:length(shotgun_fp), function(i) {
  df <- read_tsv(shotgun_fp[i])
  df <- df %>% mutate(Cohort = shotgun_prefix[i])
  return(df)
}) %>% bind_rows()

# Filter the dataframe based on values in "feature" column matching `mets`
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% spf_mets)%>% 
  filter(metadata == "Site") %>% 
  filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")
filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Distal_Colon", "Jejunum"))
filtered_shotgun_data <-filtered_shotgun_data %>% 
  mutate(data_type = "Shotgun")
# Bind the dataframe together
full_df <- rbind(filtered_shotgun_data, filtered_data)
full_df$annotation <- names(spf_mets)[match(full_df$feature, spf_mets)]
full_df$feature <- gsub("X","", full_df$feature)
full_df$feature <- gsub("\\.","-", full_df$feature)
full_df$Cohort <- factor(full_df$Cohort, 
                         levels=c("UCLA_O_SPF","CS_SPF",
                                  "SPF_Gavage",
                                  "HUM_SD_Gavage","HUM_MD_Gavage"))

aa <- full_df %>% filter(annotation=="Amino Acids")
lipids <- full_df %>% filter(annotation=="Lipids")
org <- full_df %>% filter(annotation=="Small Organic Compounds")
vit <- full_df %>% filter(annotation=="Vitamins and Cofactors")

cols <- viridis::viridis(4)
names(cols) <- c("Colon", "Distal_Colon",
                 "Jejunum", "SI")
aa_plot <- ggplot(aa, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Amino Acid Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))
aa_plot

lipids_plot <- ggplot(lipids, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Lipid Metabolites") +
  theme_cowplot(12)+
  #theme(strip.text = element_text(
  #size = 8))+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))
lipids_plot
soc <- ggplot(org, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Small Organic Compounds") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))
soc


dev.new(width=10,height=10)
plot_grid(aa_plot,soc,
         labels=c("A","B"),
         label_size = 20)

dev.new(width=10,height=10)
plot_grid(lipids_plot,
          labels=c("C"),
          label_size = 20)


#########################################

### Do this for Gavage metabolites only --- 

#########################################
gavage_mets

all_data <- lapply(1:length(file_paths), function(i) {
  df <- read_tsv(file_paths[i])
  df <- df %>% mutate(Cohort = cohort_prefixes[i])
  return(df)
}) %>% bind_rows()


# Filter the dataframe based on values in "feature" column matching `mets`
filtered_data <- all_data %>% filter(feature %in% gavage_mets)%>% 
  filter(metadata == "Site_General") #%>%
#filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")

# Create new column "Site" based on "coef"
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "Colon", "SI"))
filtered_data <- filtered_data %>% 
  mutate(data_type = "16S")

### add shotgun data to big plot ---

shotgun_fp <- c("Shotgun/melonnpan/UCLA_O_SPF/Mets-LineSexSite_General-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/CS_SPF/Mets-SexSite-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/HUM_SD_Gavage/Mets-SexSite-1-MsID/all_results.tsv",
                "Shotgun/melonnpan/SPF_Gavage/Mets-SexSite-1-MsID/all_results.tsv")

shotgun_prefix <- c("UCLA_O_SPF",
                    "CS_SPF",
                    "HUM_SD_Gavage",
                    "SPF_Gavage")


all_taxa <- process_results_for_upset_plot(file_paths = file_paths,
                                           cohort_prefixes = cohort_prefixes,
                                           filter_by = "Site")

# Given 0 significant metabolites pick out the features mapping to mets
all_shotgun_data <- lapply(1:length(shotgun_fp), function(i) {
  df <- read_tsv(shotgun_fp[i])
  df <- df %>% mutate(Cohort = shotgun_prefix[i])
  return(df)
}) %>% bind_rows()

# Filter the dataframe based on values in "feature" column matching `mets`
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% gavage_mets)%>% 
  filter(metadata == "Site") #%>% 
#filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")
filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Distal_Colon", "Jejunum"))
filtered_shotgun_data <-filtered_shotgun_data %>% 
  mutate(data_type = "Shotgun")
# Bind the dataframe together
full_df <- rbind(filtered_shotgun_data, filtered_data)
full_df$annotation <- names(gavage_mets)[match(full_df$feature, gavage_mets)]
full_df$feature <- gsub("X","", full_df$feature)
full_df$feature <- gsub("\\.","-", full_df$feature)
full_df$Cohort <- factor(full_df$Cohort, 
                         levels=c("UCLA_O_SPF","CS_SPF",
                                  "SPF_Gavage",
                                  "HUM_SD_Gavage","HUM_MD_Gavage"))

aa <- full_df %>% filter(annotation=="Amino Acid")
lipids <- full_df %>% filter(annotation=="Lipid")
org <- full_df %>% filter(annotation=="Small Organic Compound")
cols <- viridis::viridis(4)
names(cols) <- c("Colon", "Distal_Colon",
                 "Jejunum", "SI")
aa_plot <- ggplot(aa, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Amino Acid Metabolites") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))
aa_plot

lipids_plot <- ggplot(lipids, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Lipid Metabolites") +
  theme_cowplot(12)+
  #theme(strip.text = element_text(
  #size = 8))+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))

lipids_plot
soc <- ggplot(org, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
  geom_bar(stat = "identity",position = "dodge") +
  geom_hline(yintercept = 0, linetype="dashed")+
  coord_flip() +
  facet_wrap(~ Cohort,nrow = 1) +
  labs(x = "", y = "Effect size (SI/Colon)", title = "Small Organic Compounds") +
  theme_cowplot(12)+
  theme(legend.position = "top",legend.justification = "center") +
  #theme(axis.text.y = element_text(face = "italic")) +
  theme(legend.background = element_rect(fill="lightblue", size=0.5, linetype="solid")) +
  scale_fill_manual(values = cols) +
  theme(plot.title = element_text(hjust = 0.5))




bottom <- plot_grid(lipids_plot,
                    soc, labels=c("C","D"),
                    rel_widths = c(1,0.85),
                    label_size = 20)
top <- plot_grid(mets_upset,
                 aa_plot, labels=c("A","B"), 
                 rel_widths = c(1,0.85),
                 label_size = 20)

dev.new(width=10,height=10)
top

