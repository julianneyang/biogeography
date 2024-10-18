echo "Enter input filepath for qza table" $1
echo "Enter the minimum number of reads a sample must have in order to be retained" $2
echo "Enter minimum number of samples a taxon must be observed in to be retained" $3
echo "Enter the output name for your new qza file with low sequencing depth samples discarded" 
read output
qiime feature-table filter-samples --i-table $1 --p-min-frequency $2  --o-filtered-table $output
echo "Enter the output name for qza file with low depth samples discarded and prevalence filtered taxa"
read output2 
qiime feature-table filter-features --i-table $output --p-min-samples $3 --o-filtered-table $output2 
 

