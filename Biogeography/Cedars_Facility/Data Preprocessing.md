```shell
conda activate qiime2-2020.6

```
Removed all index files 
```shell
rm -rf *_I*
```
Removed all Tg mice 
```shell
rm -rf *B_*
```
Manually moved files 1BDTT13 - 1BDTT18 to the folder since they belong to WT mouse 3A

Replace all underscores with hyphens for use in QIIME2
```shell
for file in *; do mv "$file" `echo $file | tr '_' '-'` ; done
```

Export and then make a manifest file 
```shell
ls -d "$PWD"/* > filepaths.tsv

```
Run part 1 of deblur workflow
```shell
bash deblur.sh CS_Facility_manifest.csv
```

Run part 2 of deblur workflow (neg and pos filter)
```shell
bash deblur-denoise.sh CS_Facility_manifest.csv_demux_joined_filtered.qza 250

```