### Preprocessing details 




### making of taxa barplots 
- for Mouse dataset 
	- Separate by Donor ID (14 total)
	- Split each Donor ID dataset into Mucosal and Luminal subsets 
	- Group each Mucosal and Luminal subset data by Site 
	- use `taxa_barplot.sh` to generate .qzv file of Mucosal and Luminal files grouped by Site 


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
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis$ bash ../../fast-16s-analysis/shell_scripts/groupsamples.sh Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza Donors_Metadata.tsv Site
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Luminal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Donors-Analysis$ bash ../../fast-16s-analysis/shell_scripts/groupsamples.sh Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza Donors_Metadata.tsv Site
Call bash groupsamples.sh asv-table.qza metadata-file metadata-column
Saved FeatureTable[Frequency] to: groupby_Site_Mucosal_Donors-Mice-1xPrev0.15-ComBat-ASV.qza

```