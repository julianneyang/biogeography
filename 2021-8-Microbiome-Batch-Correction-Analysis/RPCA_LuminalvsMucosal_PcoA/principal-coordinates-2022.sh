echo "enter file path for distance matrix .qza" $1
qiime diversity pcoa --i-distance-matrix $1 --o-pcoa "${1%.qza}_pcoa.qza" 

