echo "Enter filepath to filtered ASV table" $1
qiime diversity beta --i-table $1 --p-metric braycurtis --o-distance-matrix "${1%.qza}_DM.qza"
qiime diversity pcoa --i-distance-matrix "${1%.qza}_DM.qza" --o-pcoa "${1%.qza}_pcoa.qza"
qiime tools export --input-path "${1%.qza}_pcoa.qza" --output-path  "${1%.qza}_pcoa/"
