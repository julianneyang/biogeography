library(janitor)
library(stringi)
library(stringr)
library(funrar)
library(lessR)
library(ggplot2)
library(tidyr)
library(cowplot)
library(dplyr)
library(plyr)
library(Microbiome.Biogeography)
here::i_am("MouseBiogeography-RProj/Shotgun-GMM-SITE-Barplot.R")


#Query the target vector against all_results.tsv and generate a heatmap 
cols=c("#440154FF","#FDE725FF")

shotgun_lum_ucla <- generate_interregional_GMM_barplot_shotgun("Shotgun/UCLA_O_SPF/GMM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/significant_results.tsv",
                                                         "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                         ystring="metabolic_map",
                                                         titlestring= "UCLA O. SPF Shotgun",
                                                         cols)

shotgun_lum_cs <- generate_interregional_GMM_barplot_shotgun("Shotgun/CS_SPF/GMM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv",
                                                               "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                               ystring="metabolic_map",
                                                               titlestring= "CS SPF Shotgun",
                                                               cols)

shotgun_lum_spf <- generate_interregional_GMM_barplot_shotgun("Shotgun/SPF_Gavage/GMM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv",
                                                             "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                             ystring="metabolic_map",
                                                             titlestring= "SPF Gavage Shotgun",
                                                             cols)

shotgun_lum_hum <- generate_interregional_GMM_barplot_shotgun("Shotgun/HUM_Gavage/GMM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/significant_results.tsv",
                                                              "Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/GOMIXER/Revised_Module_Key.csv",
                                                              ystring="metabolic_map",
                                                              titlestring= "HUM Gavage Shotgun",
                                                              cols)

### Final Figure ---
dev.new()
plot_grid(shotgun_lum_ucla, shotgun_lum_hum,
          shotgun_lum_cs, shotgun_lum_spf,
          nrow=2, ncol=2,
          labels=c("A","B","C","D"),
          label_size =16,
          rel_heights = c(1,1,0.25,0.25))

