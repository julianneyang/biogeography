echo "Enter filepath to the .qza table" $1 
echo "Enter filepath to the metadata file in tsv format" $2 
echo "Enter the column name by which to subset the data" $3 
echo "Enter the value in the column by which to subset the data" $4 
echo "Enter the output name for the subsetted .qza file" $output 
read output 
qiime feature-table filter-samples --i-table $1 --m-metadata-file $2 --p-where "[$3]='$4'" --o-filtered-table $output 
