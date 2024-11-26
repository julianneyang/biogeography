# Mouse Biogeography 

Scripts used in  "Biogeographical distribution of gut microbiome composition and function is only partially recapitulated by fecal transplantation into germ-free mice".

Raw data can be found in the NCBI database: http://www.ncbi.nlm.nih.gov/bioproject/944800.

Repository layout:
- 
- `global_phyla_cols.RDS` corresponds to the named color vector used for all phyla stacked column charts in this manuscript.
- `global_genera_cols.RDS` corresponds to the named color vector used for all genera stacked column charts in this manuscript.


To clone this repository: 
```bash
git clone https://github.com/julianneyang/biogeography
```

To reproduce the figures in this repository: 

- Navigate to `MouseBiogeography-RProj/Final_Figures`
  - Each of the following scripts corresponds to the main figures:
    - Figure 2: Figure_2_Taxa_Barplots_Aggregated.R
    - Figure 3: Figure_3_Luminal_Alpha_Beta_Diversity.R
    - Figure 4: Figure_4_Genus_Site_Heatmaps_Clustering.R
      - Note that the final versions of these heatmaps can be found in `Figure_4_Heatmaps` as cosmetic changes were made in Inkscape or Adobe Illustrator
    - Figure 6: Figure_6_GMM_Coef_Plots.R
    - Figure 7: Figure_7_GBM_Coef_Plots.R
