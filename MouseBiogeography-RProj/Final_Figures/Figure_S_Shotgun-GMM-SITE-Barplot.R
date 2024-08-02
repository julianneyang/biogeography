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
                                                              titlestring= "HUM SD Gavage Shotgun",
                                                              cols)

### Final Figure ---
dev.new()
plot_grid(shotgun_lum_ucla,  shotgun_lum_cs,
          shotgun_lum_spf,shotgun_lum_hum,
          nrow=2, ncol=2,
          labels=c("A","B","C","D"),
          label_size =16,
          rel_heights = c(1,1,0.25,0.25))

shotgun_lum_ucla <- generate_interregional_GBM_barplot_shotgun("Shotgun/UCLA_O_SPF/GBM-DCvsJej-CLR-UCLA-ComBat-SeqRunLineSexSite-1-MsID/all_results.tsv",
                                                               "GBM_Module_Key.csv",
                                                               colorvector = cols,
                                                               titlestring = "UCLA O. SPF Shotgun")
shotgun_lum_ucla$plot

shotgun_lum_cs <- generate_interregional_GBM_barplot_shotgun("Shotgun/CS_SPF/GBM-DCvsJej-CLR-CS-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                               "GBM_Module_Key.csv",
                                                               colorvector = cols,
                                                               titlestring = "UCLA O. SPF Shotgun")
shotgun_lum_cs$plot

shotgun_lum_spf <- generate_interregional_GBM_barplot_shotgun("Shotgun/SPF_Gavage/GBM-DCvsJej-CLR-SPF-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                              "GBM_Module_Key.csv",
                                                              colorvector = cols,
                                                              titlestring= "SPF Gavage Shotgun")
shotgun_lum_spf$plot

shotgun_lum_hum <- generate_interregional_GBM_barplot_shotgun("Shotgun/HUM_Gavage/GBM-DCvsJej-CLR-HUM-ComBat-SeqRunSexSite-1-MsID/all_results.tsv",
                                                             "GBM_Module_Key.csv",
                                                             colorvector = cols,
                                                             titlestring= "SPF Gavage Shotgun")
shotgun_lum_hum$plot
### Final Figure ---
dev.new()
plot_grid(shotgun_lum_ucla$plot,  shotgun_lum_cs$plot,
           shotgun_lum_spf$plot,shotgun_lum_hum$plot,
          nrow=2, ncol=2,
          labels=c("A","B","C","D"),
          label_size =16,
          rel_heights = c(1,1,0.25,0.25))

