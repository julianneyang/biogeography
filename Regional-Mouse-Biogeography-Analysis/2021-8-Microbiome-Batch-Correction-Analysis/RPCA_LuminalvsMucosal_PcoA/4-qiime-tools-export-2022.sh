echo "Enter the path to the qza input file." $1 
qiime tools export --input-path $1 --output-path  "${1%.qza}_export"
cd $output
biom summarize-table -i feature-table.biom -o biom-summary.txt 
biom convert -i feature-table.biom -o feature-table.tsv --to-tsv --header-key taxonomy
 

