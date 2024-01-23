```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/1-qiime-tools-import.sh Colonized-ComBat-Adjusted-ASV.tsv 
Enter filepath to tsv containing ASV count data Colonized-ComBat-Adjusted-ASV.tsv
Imported Colonized-ComBat-Adjusted-ASV.biom as BIOMV210Format to Colonized-ComBat-Adjusted-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Colonized-ComBat-Adjusted-ASV.qza Humanized-Metadata.tsv Microbiota Cedars_SPF
Enter filepath to the .qza table Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format Humanized-Metadata.tsv
Enter the column name by which to subset the data Microbiota
Enter the value in the column by which to subset the data Cedars_SPF
Saved FeatureTable[Frequency] to: Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Colonized-ComBat-Adjusted-ASV.qza Humanized-Metadata.tsv Microbiota Humanized
Enter filepath to the .qza table Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format Humanized-Metadata.tsv
Enter the column name by which to subset the data Microbiota
Enter the value in the column by which to subset the data Humanized
Saved FeatureTable[Frequency] to: Humanized_Colonized-ComBat-Adjusted-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files$ bash ../../../fast-16s-analysis/shell_scripts/import-taxonomy.sh Humanized_taxonomy.tsv 
Imported Humanized_taxonomy.tsv as TSVTaxonomyDirectoryFormat to taxonomy.qza

```
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files/Site_Subsets$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh $file ../Humanized-Metadata.tsv Type Mucosal; done
Enter filepath to the .qza table Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format ../Humanized-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal
Saved FeatureTable[Frequency] to: Mucosal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the .qza table Humanized_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format ../Humanized-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal
Saved FeatureTable[Frequency] to: Mucosal_Humanized_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the .qza table Luminal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format ../Humanized-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal
Saved FeatureTable[Frequency] to: Mucosal_Luminal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the .qza table Luminal_Humanized_Colonized-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format ../Humanized-Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal
Saved FeatureTable[Frequency] to: Mucosal_Luminal_Humanized_Colonized-ComBat-Adjusted-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files/Site_Subsets$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/groupsamples.sh $file ../Humanized-Metadata.tsv Site;done
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Luminal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Luminal_Humanized_Colonized-ComBat-Adjusted-ASV.qza
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Mucosal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Mucosal_Humanized_Colonized-ComBat-Adjusted-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Humanized-Biogeography-Analysis/starting_files/Site_Subsets/taxa_barplots$ for file in *; do bash ../../../../../fast-16s-analysis/shell_scripts/taxabarplot.sh $file ../../taxonomy.qza site_metadata.tsv ;done
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Luminal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza_dir.qzv
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Luminal_Humanized_Colonized-ComBat-Adjusted-ASV.qza_dir.qzv
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Mucosal_Cedars_SPF_Colonized-ComBat-Adjusted-ASV.qza_dir.qzv
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Mucosal_Humanized_Colonized-ComBat-Adjusted-ASV.qza_dir.qzv

```