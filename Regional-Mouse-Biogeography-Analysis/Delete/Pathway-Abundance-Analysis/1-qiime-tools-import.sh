echo "Enter filepath to tsv containing ASV count data" $1
echo "Enter name of output .biom file" 
read output
biom convert -i $1 -o $output --table-type="OTU table" --to-hdf5 
echo "Enter name of output .qza file"
read output2
qiime tools import --input-path $output --type 'FeatureTable[Frequency]' --input-format BIOMV210Format --output-path $output2
echo "Enter name of biom summary .txt file" 
read output3
biom summarize-table -i $output -o $output3 

