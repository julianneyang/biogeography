echo "Enter filepath to ASV count table in .qza format" $1 
echo "Enter filepath to sample metadata file in .tsv format" $2 
echo "Enter the name of the column in the metadata that indicates the individual subject/site (e.g. subject ID) that was sampled repeatedly" $3 
echo "Enter the name of the column in the metadata that indicates the numeric repeated measure (e.g., Time in months/days) or non-numeric category" $4 
echo "Enter the name of the output directory" 
read output 
qiime gemelli ctf  --i-table  $1 --m-sample-metadata-file $2 --p-state-column $4 --p-individual-id-column $3 --output-dir $output
