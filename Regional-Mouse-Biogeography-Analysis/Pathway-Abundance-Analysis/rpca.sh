echo "Enter filepath of filtered ASV table.qza file" $1 
echo "Enter output name of Robust Aitchison Biplot"
read output2
echo "Enter output name of Robust Aitchison Distance Matrix"
read output3 

qiime deicode rpca --i-table $1 --o-biplot $output2 --o-distance-matrix $output3
