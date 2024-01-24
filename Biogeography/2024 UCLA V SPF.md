```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/UCLA_V_SPF_Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/1-qiime-tools-import.sh UCLA_V_SPF_ComBat_Adjusted_ASV.tsv 
Enter filepath to tsv containing ASV count data UCLA_V_SPF_ComBat_Adjusted_ASV.tsv
Imported UCLA_V_SPF_ComBat_Adjusted_ASV.biom as BIOMV210Format to UCLA_V_SPF_ComBat_Adjusted_ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/UCLA_V_SPF_Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/import-taxonomy.sh ucla_v_taxonomy.tsv 
Imported ucla_v_taxonomy.tsv as TSVTaxonomyDirectoryFormat to taxonomy.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/UCLA_V_SPF_Analysis/starting_files$ nano site_metadata.tsv
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/UCLA_V_SPF_Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Site_UCLA_V_SPF_ComBat_Adjusted_ASV.qza taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_UCLA_V_SPF_ComBat_Adjusted_ASV.qza_dir.qzv

```