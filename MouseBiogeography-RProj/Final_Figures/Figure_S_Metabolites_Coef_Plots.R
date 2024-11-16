###### The Final Frontier of Mouse Biogeography ---
### Date: 10.18.24
### Figure Number: Supplementary #? 
### Figure Contents: Metabolites Coefficient Pltos  
###### whining ends here ---

library(ggplot2)
library(tidyverse)
#library(rlang)
library(cowplot)
library(viridis)
#library(plyr)
library(gridExtra)
#library(paletteer)
library(ComplexUpset)
library(here)

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
mets_upset <-upset(id_df_wide, all_datasets,width_ratio=0.15, name="",
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

spf_2_out_of_3 <- id_df_wide %>%
  mutate(SPF_only = case_when(
    UCLA_O_SPF == 1 & CS_SPF == 1 &
      rowSums(select(., HUM_MD_Gavage, HUM_SD_Gavage, SPF_Gavage) == 0) >= 2 ~ "SPF",
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

spf_relaxed <- spf_2_out_of_3 %>% 
  filter(SPF_only == "SPF")

length(spf_filtered$feature)
length(unique(spf_filtered$feature))

unique_spf <- spf_filtered$feature
unique_spf_relaxed <- spf_relaxed$feature
unique_gavage <- gavage_filtered$feature
write_rds(unique_spf,here("SPF_only_metabolites.RDS"))
write_rds(unique_gavage, here("Gavage_only_metabolites.RDS"))
write_rds(unique_spf_relaxed, here("SPF_2outof3_metabolites.RDS"))

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

spf_metab <- readRDS(here("SPF_only_metabolites.RDS"))
names(spf_metab) <- c(
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
  "Lipids",                        # 1 decanedioic.acid
  "Lipids",                        # 2 TG.43.1
  "Amino Acids",                  # 3 Tryptoline
  "Amino Acids",                  # 4 carnitine
  "Lipids",                        # 5 campesterol
  "Small Organic Compounds",       # 6 deprenyl
  "Lipids",                        # 7 X.Z..5.8.11.trihydroxyoctadec.9.enoic.acid
  "Lipids",                        # 8 SM.22.1.2O.12.0
  "Lipids",                        # 9 FA.24.0
  "Small Organic Compounds",       # 10 X5.Aminosalicylic.acid
  "Small Organic Compounds",       # 11 X6.Aminoindazole
  "Small Organic Compounds",       # 12 X13.KODE
  "Lipids",                        # 13 FA.25.0
  "Lipids",                        # 14 Cer.NDS.d41.0
  "Small Organic Compounds",       # 15 NAGlySer.31.1.O2.NAGlySer.16.0.O.FA.15.0.
  "Small Organic Compounds",       # 16 PGE2
  "Small Organic Compounds",       # 17 X8.Oxo.2.deoxyadenosine
  "Lipids",                        # 18 FA.27.0
  "Small Organic Compounds",       # 19 X7.Methylguanine
  "Lipids",                        # 20 FA.26.0
  "Lipids",                        # 21 MGDG.17.0_18.2
  "Small Organic Compounds",       # 22 X2..2.6.dimethylmorpholin.4.yl.ethanol
  "Lipids",                        # 23 FA.22.0
  "Small Organic Compounds",       # 24 Isorhamnetin.3.O.rutinoside
  "Lipids",                        # 25 X4.HDoHE
  "Small Organic Compounds",       # 26 Biliverdin
  "Small Organic Compounds",       # 27 X5.Acetylamino.6.amino.3.methyluracil
  "Lipids",                        # 28 PG.33.0.PG.16.0_17.0
  "Small Organic Compounds",       # 29 X6.deoxyglucose
  "Small Organic Compounds",       # 30 Caffeic.acid
  "Lipids",                        # 31 SM.d34.1
  "Lipids",                        # 32 lauric.acid
  "Small Organic Compounds",       # 33 arabitol
  "Lipids",                        # 34 SM.d36.0.Isomer.B
  "Lipids",                        # 35 CAR.16.0
  "Small Organic Compounds",       # 36 lyxitol
  "Lipids",                        # 37 PC.34.0
  "Lipids",                        # 38 FA.24.2
  "Lipids",                        # 39 SM.24.1
  "Lipids",                        # 40 TG.61.3.TG.18.0_25.0_18.3
  "Lipids",                        # 41 CAR.18.0.Isomer.A
  "Lipids",                        # 42 PC.38.4.Isomer.B
  "Small Organic Compounds",       # 43 NAGlySer.30.1.O2.NAGlySer.15.0.O.FA.15.0.
  "Lipids",                        # 44 FA.23.0
  "Lipids",                        # 45 FA.28.0
  "Small Organic Compounds",       # 46 X2..deoxyguanosine
  "Amino Acids",                  # 47 N.Acetylserine
  "Lipids",                        # 48 Cer.d39.1
  "Lipids",                        # 49 TG.54.5.O2.TG.18.1_18.2_18.2.O2
  "Lipids",                        # 50 Nervonic.acid
  "Lipids",                        # 51 LPE.22.6
  "Small Organic Compounds",       # 52 N.N..Diacetylcystine
  "Lipids",                        # 53 LPG.16.0
  "Lipids",                        # 54 PE.34.0
  "Lipids",                        # 55 Cer.43.0.O2.Cer.17.0.O2.26.0
  "Lipids",                        # 56 Cer.d38.1
  "Lipids",                        # 57 PC.36.4.Isomer.B
  "Lipids",                        # 58 Cer.35.0.O3.Cer.18.0.O2.17.0.O
  "Lipids",                        # 59 TG.34.0
  "Amino Acids",                  # 60 tyramine
  "Small Organic Compounds",       # 61 X2..Deoxyinosine
  "Small Organic Compounds",       # 62 X4.Pentylaniline
  "Small Organic Compounds",       # 63 Histamine
  "Lipids",                        # 64 FA.29.1
  "Lipids",                        # 65 FA.27.1
  "Lipids",                        # 66 FA.33.1
  "Small Organic Compounds",       # 67 inosine
  "Amino Acids",                  # 68 Met.His
  "Lipids",                        # 69 TG.45.0
  "Lipids",                        # 70 Phosphatidylethanolamine.alkenyl.16.0.20.4
  "Amino Acids",                  # 71 delta.Hydroxylysine
  "Small Organic Compounds",       # 72 Tomatidine
  "Small Organic Compounds",       # 73 X2.deoxyuridine
  "Small Organic Compounds",       # 74 dihydro.3.coumaric.acid
  "Lipids",                        # 75 FA.31.1
  "Lipids",                        # 76 Cer.d34.2.Isomer.A
  "Lipids",                        # 77 PE.36.3
  "Small Organic Compounds",       # 78 X2.Amino.3.methoxybenzoic.acid
  "Small Organic Compounds",       # 79 p.coumaric.acid
  "Small Organic Compounds",       # 80 X2..n.ethyl.n.m.toluidino.ethanol
  "Small Organic Compounds",       # 81 X3..3.Hydroxyphenyl.propionic.acid
  "Small Organic Compounds",       # 82 Mandelic.acid
  "Lipids",                        # 83 FA.32.0
  "Small Organic Compounds",       # 84 X1.Methylnicotinamide
  "Lipids",                        # 85 FA.25.1
  "Amino Acids",                  # 86 N.fructosyl.isoleucine
  "Small Organic Compounds",       # 87 cadaverine
  "Amino Acids",                  # 88 Val.Val.Lys
  "Lipids",                        # 89 SM.d40.0.Isomer.B
  "Lipids",                        # 90 Cer.18.0.2O_20.0
  "Lipids",                        # 91 FA.35.1
  "Lipids",                        # 92 SHexCer.34.1.3O
  "Lipids",                        # 93 FA.34.0
  "Amino Acids",                  # 94 Gln.Tyr
  "Lipids",                        # 95 GlcCer.d40.1
  "Lipids",                        # 96 FA.37.1
  "Lipids",                        # 97 Cholic.acid
  "Lipids",                        # 98 PE.32.1.PE.16.0_16.1
  "Amino Acids",                  # 99 Arg.Trp
  "Small Organic Compounds",       # 100 FA.20.5..eicosapentaenoic.acid.
  "Small Organic Compounds",       # 101 Choline
  "Small Organic Compounds",       # 102 guanosine
  "Lipids",                        # 103 Taurocholic.acid
  "Amino Acids",                  # 104 Ile.Glu
  "Lipids",                        # 105 PE.O.16.1_22.6
  "Small Organic Compounds",       # 106 X3..3.hydroxyphenyl.propionic.acid
  "Amino Acids",                  # 107 Tyr.Thr
  "Amino Acids",                  # 108 Gln.Pro
  "Lipids",                        # 109 FA.22.4
  "Lipids",                        # 110 LPC.O.22.0
  "Lipids",                        # 111 PI.36.2
  "Small Organic Compounds",       # 112 X2.Salicylic.acid
  "Small Organic Compounds",       # 113 3..phenylalanine
  "Small Organic Compounds",       # 114 4..Cyanocobalamin
  "Lipids",                        # 115 PI.34.2
  "Small Organic Compounds",       # 116 2.3..3..beta..carotene
  "Amino Acids",                  # 117 Histamine.5..
  "Small Organic Compounds",       # 118 Histidine
  "Lipids",                        # 119 PG.29.0.PG.18.1_12.1
  "Amino Acids",                  # 120 His.Lys
  "Lipids",                        # 121 PG.33.0.PG.18.0_16.0
  "Lipids",                        # 122 Sphingomyelin.d38.2
  "Small Organic Compounds",       # 123 N.Acetylspermidine
  "Lipids",                        # 124 PE.36.0
  "Amino Acids",                  # 125 Glu.His
  "Amino Acids",                  # 126 Val.Ala
  "Lipids",                        # 127 Acylcarnitine.d35.1
  "Lipids",                        # 128 PC.37.3
  "Lipids",                        # 129 Cer.25.0
  "Amino Acids"                   # 130 Gly.Gln
)

cols <- viridis::viridis(4)
names(cols) <- c("16S_Colon", "Shotgun_DC",
                 "Shotgun_Jej", "16S_SI")
#######################################################
### Do this for Metabolites concordant across all 5 ---
#######################################################

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
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "16S_Colon", "16S_SI"))
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
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% mets)%>% 
  filter(metadata == "Site") 

filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Shotgun_DC", "Shotgun_Jej"))
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

bottom <- plot_grid(lipids_plot,
                 soc, labels=c("C","D"),
                 rel_widths = c(1,0.85),
                 label_size = 20)
top <- plot_grid(mets_upset,
                 aa_plot, labels=c("A","B"), 
                 rel_widths = c(1,0.85),
                 label_size = 20)

dev.new(width=10,height=10)
bottom

write_rds(top, here("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Top.RDS"))
write_rds(bottom, here("MouseBiogeography-RProj/Final_Figures/Figure_S_Metabolites_Bottom.RDS"))

#########################################

### Do this for SPF metabolites only --- 

#########################################
spf_relaxed_mets <- readRDS(here("SPF_2outof3_metabolites.RDS"))
not_amino_acids <- c("Myristoyl.L.carnitine", 
                     "Tomatidine", 
                     "Catechol", 
                     "Uracil", 
                     "guanosine", 
                     "Choline", 
                     "dibenzylamine", 
                     "Adenine", 
                     "pyrrolidine", 
                     "Mannitol")
classify_metabolite <- function(metabolite) {
  if (grepl("acid$", metabolite)) {
    return("Organic Acids")
  } else if (grepl("FA", metabolite)) {
    return("Fatty Acids")
  } else if (grepl("^PC\\.|^PE\\.|^PS\\.|^LPC\\.|^DGDG\\.|^PG\\.", metabolite)) {
    return("Phospholipids")
  } else if (grepl("^TG\\.", metabolite)) {
    return("Triacylglycerols")
  } else if (grepl("^Cer\\.|^HexCer\\.|^SM\\.|^SL\\.",  metabolite)) {
    return("Sphingolipids")
  } else if (!metabolite %in% not_amino_acids && (grepl("^[A-Z][a-z]*$", metabolite) ||
             grepl("ine$", metabolite))) {  # Added condition for amino acids ending with "ine"
    return("Amino Acids and Derivatives")  # Added condition for amino acids ending with "ine"
  } else if (metabolite %in% c("tryptophan")) {
    return("Amino Acids and Derivatives")
  } else if (grepl("^[A-Z][a-z]*\\.[A-Z][a-z]*$", metabolite) ||
             grepl("^[A-Z][a-z]*\\.[A-Z][a-z]*\\.[A-Z][a-z]*$", metabolite) ||
             metabolite %in% c("tyr.arg")) {
    return("Dipeptides and Tripeptides")
  } else if (metabolite %in% not_amino_acids) {
    return("Small Organic Compounds")
  } else {
    return("Small Organic Compounds")
  }
}

classification <- sapply(spf_relaxed_mets, classify_metabolite)
names(spf_relaxed_mets) <- classification

spf_mets <- readRDS(here("SPF_only_metabolites.RDS"))
classification <- sapply(spf_mets, classify_metabolite)
names(spf_mets) <- classification

spf_mets <- spf_relaxed_mets
#spf_mets <- spf_metab

all_data <- lapply(1:length(file_paths), function(i) {
  df <- read_tsv(file_paths[i])
  df <- df %>% mutate(Cohort = cohort_prefixes[i])
  return(df)
}) %>% bind_rows()


# Filter the dataframe based on values in "feature" column matching `mets`
filtered_data <- all_data %>% filter(feature %in% spf_mets)%>% 
  filter(metadata == "Site_General") %>% 
  filter(qval < 0.05) #%>%
  #filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")

# Create new column "Site" based on "coef"
filtered_data <- filtered_data %>% mutate(Site = ifelse(coef < 0, "16S_Colon", "16S_SI"))
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


# Given 0 significant metabolites pick out the features mapping to mets
all_shotgun_data <- lapply(1:length(shotgun_fp), function(i) {
  df <- read_tsv(shotgun_fp[i])
  df <- df %>% mutate(Cohort = shotgun_prefix[i])
  return(df)
}) %>% bind_rows()

# Filter the dataframe based on values in "feature" column matching `mets`
filtered_shotgun_data <- all_shotgun_data %>% filter(feature %in% spf_mets)%>% 
  filter(metadata == "Site") #%>% 
  #filter(Cohort=="UCLA_O_SPF" | Cohort =="CS_SPF")
filtered_shotgun_data <- filtered_shotgun_data %>% mutate(Site = ifelse(coef < 0, "Shotgun_DC", "Shotgun_Jej"))
filtered_shotgun_data <-filtered_shotgun_data %>% 
  mutate(data_type = "Shotgun")
# Bind the dataframe together
full_df <- rbind(filtered_shotgun_data, filtered_data)
full_df$annotation <- names(spf_mets)[match(full_df$feature, spf_mets)]
full_df <- full_df %>% filter(feature!="Guanosine") #duplicate entry
full_df$feature <- gsub("X","", full_df$feature)
full_df$feature <- gsub("\\.","-", full_df$feature)
full_df$Cohort <- factor(full_df$Cohort, 
                         levels=c("UCLA_O_SPF","CS_SPF",
                                  "SPF_Gavage",
                                  "HUM_SD_Gavage","HUM_MD_Gavage"))

cols <- viridis::viridis(4)
names(cols) <- c("16S_Colon", "Shotgun_DC",
                 "Shotgun_Jej", "16S_SI")

# Plot metabolites
plot_metabolites <- function(data, title) {
  ggplot(data, aes(x = reorder(feature, coef), y = coef, fill = Site)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    coord_flip() +
    facet_wrap(~ Cohort, nrow = 1) +
    labs(x = "", y = "Effect size (SI/Colon)", title = title) +
    theme_cowplot(12) +
    theme(legend.position = "top", legend.justification = "center") +
    #theme(legend.background = element_rect(fill = "lightblue", size = 2, linetype = "solid")) +
    scale_fill_manual(values = cols) +
    theme(plot.title = element_text(hjust = 0.5))+
    scale_y_continuous(labels=scales::label_number(accuracy = 0.1)) 

}

# Define the different categories
categories <- unique(full_df$annotation)
# Loop through the categories and plot
plots <- list()
for (cat in categories) {
  data <- full_df %>% filter(annotation == cat)
  title <- cat
  plots[[cat]] <- plot_metabolites(data, title)
}

dev.new(width=10,height=10)
plot_grid(plots[[3]], plots[[1]],
          labels=c("A","B"),
          label_size = 20)
dev.new(width=10,height=10)
plot_grid(plots[[2]],
          labels=c("C"),
          label_size = 20)
# Access individual plots
phospholipid <- plots[[1]] 
soc <- plots[[2]]
pep <- plots[[3]]

sphingo <- plots[[4]]
aa <- plots[[5]] 
fa <- plots[[6]]
tg <- plots[[7]]
oa<-plots[[8]]

dev.new(width=10,height=10)
plot_grid(phospholipid, sphingo,
         labels=c("A","B"),
         label_size = 20,
         rel_widths = c(0.9,1))

dev.new(width=10,height=10)
tg_whitespace <- plot_grid(NULL,tg,
                           rel_widths = c(0.24,1))
fa_tg <- plot_grid(fa,tg_whitespace,
          labels=c("D","E"),
          label_size = 20,
          rel_heights = c(2,1),
          rel_widths = c(1,0.5),
          nrow=2)


dev.new(width=10,height=10)
plot_grid(aa,fa_tg,
          labels=c("C",""),
          label_size = 20,
          rel_widths = c(0.8,1))

dev.new(width=10,height=10)
plot_grid( pep, soc,
          labels=c("A","B"),
          label_size = 20,
          rel_widths = c(1,1))

dev.new(width=10,height=10)
plot_grid(oa,
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

