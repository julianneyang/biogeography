### Preprocessing details 




### making of taxa barplots 
- for Mouse dataset 
	- Separate by Donor ID (14 total)
	- Split each Donor ID dataset into Mucosal and Luminal subsets 
	- Group each Mucosal and Luminal subset data by Site 
	- use `taxa_barplot.sh` to generate .qzv file of Mucosal and Luminal files grouped by Site 
```bash
for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/taxabarplot.sh $file ../../taxonomy.qza ../Site_Metadata.tsv; done ```
