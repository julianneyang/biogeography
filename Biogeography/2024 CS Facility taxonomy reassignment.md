```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis$ bash ../../fast-16s-analysis/shell_scripts/1-qiime-tools-import.sh CS-Facility-ComBat-Adjusted-ASV.tsv 
Enter filepath to tsv containing ASV count data CS-Facility-ComBat-Adjusted-ASV.tsv
Imported CS-Facility-ComBat-Adjusted-ASV.biom as BIOMV210Format to CS-Facility-ComBat-Adjusted-ASV.qza


```

Separate datasets into subsets by Type : Luminal
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh CS-Facility-ComBat-Adjusted-ASV.qza CS_Facility_Metadata.tsv Type Luminal
Enter filepath to the .qza table CS-Facility-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format CS_Facility_Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Luminal

```
Separate datasets into subsets by Type:Mucosal
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh CS-Facility-ComBat-Adjusted-ASV.qza CS_Facility_Metadata.tsv Type Mucosal
Enter filepath to the .qza table CS-Facility-ComBat-Adjusted-ASV.qza
Enter filepath to the metadata file in tsv format CS_Facility_Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal

```

Group samples by Site for barplot visualization 
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/groupsamples.sh Luminal_CS-Facility-ComBat-Adjusted-ASV.qza ../CS_Facility_Metadata.tsv Site
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/groupsamples.sh Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza ../CS_Facility_Metadata.tsv Site
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
```
Generate taxa barplot 
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Site_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza ../taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza_dir.qzv
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/taxabarplot.sh groupby_Site_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza ../taxonomy.qza site_metadata.tsv 
Call bash taxabarplot.sh table taxonomy metadata
Saved Visualization to: barplot_groupby_Site_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza_dir.qzv

```