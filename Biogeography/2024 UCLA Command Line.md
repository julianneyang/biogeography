```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis$ bash ../../../fast-16s-analysis/shell_scripts/1-qiime-tools-import.sh UCLA-ComBat-Adjusted-ASV.tsv 
Enter filepath to tsv containing ASV count data UCLA-ComBat-Adjusted-ASV.tsv
Imported UCLA-ComBat-Adjusted-ASV.biom as BIOMV210Format to UCLA-ComBat-Adjusted-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh UCLA-ComBat-Adjusted-ASV.qza Regional-Combat-Metadata.tsv Type Luminal
Enter filepath to the .qza table UCLA-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format Regional-Combat-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Luminal
Saved FeatureTable[Frequency] to: Luminal_UCLA-ComBat-Adjusted-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh UCLA-ComBat-Adjusted-ASV.qza Regional-Combat-Metadata.tsv Type Mucosal
Enter filepath to the .qza table UCLA-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format Regional-Combat-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal

```

```
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis$ bash ../../../fast-16s-analysis/shell_scripts/import-taxonomy.sh UCLA_taxonomy.tsv 
Imported UCLA_taxonomy.tsv as TSVTaxonomyDirectoryFormat to taxonomy.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis$ mv Luminal_UCLA-ComBat-Adjusted-ASV.qza Mucosal_UCLA-ComBat-Adjusted-ASV.qza taxa_barplots/

qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/groupsamples.sh $file ../Regional-Combat-Metadata.tsv Site; done
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Luminal_UCLA-ComBat-Adjusted-ASV.qza
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Mucosal_UCLA-ComBat-Adjusted-ASV.qza


```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots$ bash ../../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Site_Luminal_UCLA-ComBat-Adjusted-ASV.qza ../taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Luminal_UCLA-ComBat-Adjusted-ASV.qza_dir.qzv
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Regional-Mouse-Biogeography-Analysis/2021-8-Microbiome-Batch-Correction-Analysis/taxa_barplots$ bash ../../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Site_Mucosal_UCLA-ComBat-Adjusted-ASV.qza ../taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Mucosal_UCLA-ComBat-Adjusted-ASV.qza_dir.qzv

```