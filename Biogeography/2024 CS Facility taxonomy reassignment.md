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

Collapse to higher order phylogenies
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/collapse-taxa.sh Luminal_CS-Facility-ComBat-Adjusted-ASV.qza ../taxonomy.qza 
Saved FeatureTable[Frequency] to: L1_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L2_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L3_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L4_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L5_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L6_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L7_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ bash ../../../fast-16s-analysis/shell_scripts/collapse-taxa.sh Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza ../taxonomy.qza 
Saved FeatureTable[Frequency] to: L1_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L2_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L3_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L4_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L5_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L6_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Saved FeatureTable[Frequency] to: L7_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza

```

Move subsets to Site Subsets and Type Subsets and export all qza files 
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ mv L*_Luminal* Site_subsets/
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ mv L*_Mucosal* Site_subsets/
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots$ cd Site_subsets/
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots/Site_subsets$ ls
L1_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L1_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L2_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L2_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L3_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L3_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L4_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L4_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L5_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L5_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L6_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L6_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
L7_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
L7_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/CS-Facility-Analysis/Taxa-Barplots/Site_subsets$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/4-qiime-tools-export.sh $file;done
Takes qza input file as input and cranks out tsv and summary.txt file L1_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L1_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L1_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L1_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L1_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L1_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L2_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L2_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L2_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L2_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L2_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L2_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L3_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L3_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L3_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L3_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L3_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L3_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L4_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L4_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L4_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L4_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L4_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L4_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L5_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L5_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L5_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L5_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L5_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L5_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L6_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L6_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L6_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L6_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L6_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L6_Mucosal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L7_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L7_Luminal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L7_Luminal_CS-Facility-ComBat-Adjusted-ASV
Takes qza input file as input and cranks out tsv and summary.txt file L7_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza
Exported L7_Mucosal_CS-Facility-ComBat-Adjusted-ASV.qza as BIOMV210DirFmt to directory export_L7_Mucosal_CS-Facility-ComBat-Adjusted-ASV

```