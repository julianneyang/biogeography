echo "Enter filepath of filtered ASV table.qza file" $1 
echo "Enter filepath of tree.qza" $2
echo "Enter metadata file" $3
echo "Enter depth at which to rarefy" $4 
echo "Enter name of output directory" 
read output
echo "Enter output name of Robust Aitchison Biplot"
read output2
echo "Enter output name of Robust Aitchison Distance Matrix"
read output3 
qiime diversity core-metrics-phylogenetic --i-table $1 --i-phylogeny $2 --m-metadata-file $3 --p-sampling-depth $4 --output-dir $output
cd $output 
qiime tools export --input-path bray_curtis_pcoa_results.qza --output-path bray_curtis_pcoa_results
cd ..
qiime deicode rpca --i-table $1 --o-biplot $output2 --o-distance-matrix $output3

