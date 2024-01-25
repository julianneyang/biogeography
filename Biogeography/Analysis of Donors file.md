### Preprocessing details 




### making of taxa barplots 
- for Mouse dataset 
	- Separate by Donor ID (14 total)
	- Split each Donor ID dataset into Mucosal and Luminal subsets 
	- Group each Mucosal and Luminal subset data by Site 
	- use `taxa_barplot.sh` to generate .qzv file of Mucosal and Luminal files grouped by Site 
```bash
for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/taxabarplot.sh $file ../../taxonomy.qza ../Site_Metadata.tsv; done ```
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Donors-Mice-1xPrev0.15-ComBat-ASV.qza Donors_Metadata.tsv Type Luminal
Enter filepath to the .qza table Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Enter filepath to the metadata file in tsv format Donors_Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Luminal
Saved FeatureTable[Frequency] to: Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh Donors-Mice-1xPrev0.15-ComBat-ASV.qza Donors_Metadata.tsv Type Mucosal
Enter filepath to the .qza table Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Enter filepath to the metadata file in tsv format Donors_Metadata.tsv
Enter the column name by which to subset the data Type
Enter the value in the column by which to subset the data Mucosal
Saved FeatureTable[Frequency] to: Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis/site_subsets$ for file in *; do bash ../../../fast-16s-analysis/shell_scripts/collapse-taxa.sh $file ../taxonomy.qza ; done
Saved FeatureTable[Frequency] to: L1_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L2_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L3_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L4_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L5_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L6_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L7_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L1_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L2_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L3_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L4_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L5_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L6_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
Saved FeatureTable[Frequency] to: L7_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis/site_subsets$ for file in *; do bash ../../../fast-16s-analysis/shell_scripts/4-qiime-tools-export.sh $file; done

```