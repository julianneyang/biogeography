echo "Enter filepath for qza file" $1
echo "Enter minimum number of samples a taxon must be observed in to be retained" $2
echo "Enter output.qza file name"
read output 
qiime feature-table filter-features --i-table $1 --p-min-samples $2 --o-filtered-table $output 
