# Mouse Biogeography 

Scripts used in  "Biogeographical distribution of gut microbiome composition and function is only partially recapitulated by fecal transplantation into germ-free mice".

This repository concerns reproducing the figures found in the manuscript.

Raw data can be found in the NCBI database: http://www.ncbi.nlm.nih.gov/bioproject/944800.

Repository layout:
- `CS_SPF` : Cedars-Sinai SPF dataset files 
- `Donors_Analysis` : HUM MD Gavage dataset files
- `Humanized-Biogeography-Analysis` : HUM SD Gavage and SPF Gavage dataset files
- `Regional-Mouse-Biogeography-Analysis` : UCLA O. SPF Dataset files
- `Shotgun` : shotgun metagenomics for 4 datsets
- `UCLA_V_SPF_Analysis` : UCLA V. SPF dataset files
- `melonnpan_model` : Trained model used for predicting metabolites from KO microbiome data 
-  Miscellaneous files    
	- `global_phyla_cols.RDS` corresponds to the named color vector used for all phyla stacked column charts in this manuscript.
	- `global_genera_cols.RDS` corresponds to the named color vector used for all genera stacked column charts in this manuscript.


To clone this repository: 
```bash
git clone https://github.com/julianneyang/biogeography
```

To reproduce the figures in this repository: 

- Navigate to `MouseBiogeography-RProj/Final_Figures` and run each script to produce the corresponding figure.
  - Note that filepaths are built to the top level of the repo.
  - Each of the following scripts corresponds to the main figures:
    - Figure 2: Figure_2_Taxa_Barplots_Aggregated.R
    - Figure 3: Figure_3_Luminal_Alpha_Beta_Diversity.R
    - Figure 4: Figure_4_Genus_Site_Heatmaps_Clustering.R
      - Note that the final versions of these heatmaps can be found in `Figure_4_Heatmaps` as cosmetic changes were made in Inkscape or Adobe Illustrator
    - Figure 6: Figure_6_GMM_Coef_Plots.R
    - Figure 7: Figure_7_GBM_Coef_Plots.R
