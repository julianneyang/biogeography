echo "enter file path for distance matrix .qza" $1
echo "enter output name for PcoA with .qza file extension" 
read output
qiime diversity pcoa --i-distance-matrix $1 --o-pcoa $output

echo "Enter your desired output directory name for your principal coordinates .txt file" 
read output2 
qiime tools export --input-path $output --output-path  $output2

