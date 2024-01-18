
```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh BioGeo_Shotgun_ASV.qza BioGeo_Shotgun_Metadata.tsv Dataset UCLA_O_SPF
Enter filepath to the .qza table BioGeo_Shotgun_ASV.qza
Enter filepath to the metadata file in tsv format BioGeo_Shotgun_Metadata.tsv
Enter the column name by which to subset the data Dataset
Enter the value in the column by which to subset the data UCLA_O_SPF
Saved FeatureTable[Frequency] to: UCLA_O_SPF_BioGeo_Shotgun_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh BioGeo_Shotgun_ASV.qza BioGeo_Shotgun_Metadata.tsv Dataset CS_SPF
Enter filepath to the .qza table BioGeo_Shotgun_ASV.qza
Enter filepath to the metadata file in tsv format BioGeo_Shotgun_Metadata.tsv
Enter the column name by which to subset the data Dataset
Enter the value in the column by which to subset the data CS_SPF
Saved FeatureTable[Frequency] to: CS_SPF_BioGeo_Shotgun_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh BioGeo_Shotgun_ASV.qza BioGeo_Shotgun_Metadata.tsv Dataset SPF_Gavage
Enter filepath to the .qza table BioGeo_Shotgun_ASV.qza
Enter filepath to the metadata file in tsv format BioGeo_Shotgun_Metadata.tsv
Enter the column name by which to subset the data Dataset
Enter the value in the column by which to subset the data SPF_Gavage
Saved FeatureTable[Frequency] to: SPF_Gavage_BioGeo_Shotgun_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ bash ../../fast-16s-analysis/shell_scripts/filter-ASV-by-metadata.sh BioGeo_Shotgun_ASV.qza BioGeo_Shotgun_Metadata.tsv Dataset HUM_Gavage
Enter filepath to the .qza table BioGeo_Shotgun_ASV.qza
Enter filepath to the metadata file in tsv format BioGeo_Shotgun_Metadata.tsv
Enter the column name by which to subset the data Dataset
Enter the value in the column by which to subset the data HUM_Gavage
Saved FeatureTable[Frequency] to: HUM_Gavage_BioGeo_Shotgun_ASV.qza
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ mkdir Site_RPCA
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ mv *_BioGeo_* Site_RPCA/
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun$ cd Site_RPCA/

```

```bash
(qiime2-2022.2) julianne@laptop:~/Documents/biogeography/Shotgun/Site_RPCA/rpca$ for file in *; do bash ../../../../fast-16s-analysis/shell_scripts/rpca.sh $file; done
Enter filepath of filtered ASV table.qza file s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza
This will perform Robust Aitchison PCA on a count matrix
QIIME is caching your current deployment for improved performance. This may take a few moments and should only happen once per deployment.
Saved PCoAResults % Properties('biplot') to: biplot_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza
Saved DistanceMatrix to: dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza
Saved PCoAResults to: pcoa_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza
Exported pcoa_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza as OrdinationDirectoryFormat to directory pcoa_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza.txt
Exported dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza as DistanceMatrixDirectoryFormat to directory dm_rpca_s2_min400000_SPF_Gavage_BioGeo_Shotgun_ASV.qza.txt
Enter filepath of filtered ASV table.qza file s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza
This will perform Robust Aitchison PCA on a count matrix
Saved PCoAResults % Properties('biplot') to: biplot_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza
Saved DistanceMatrix to: dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza
Saved PCoAResults to: pcoa_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza
Exported pcoa_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza as OrdinationDirectoryFormat to directory pcoa_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza.txt
Exported dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza as DistanceMatrixDirectoryFormat to directory dm_rpca_s2_min50000_HUM_Gavage_BioGeo_Shotgun_ASV.qza.txt
Enter filepath of filtered ASV table.qza file s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza
This will perform Robust Aitchison PCA on a count matrix
Saved PCoAResults % Properties('biplot') to: biplot_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza
Saved DistanceMatrix to: dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza
Saved PCoAResults to: pcoa_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza
Exported pcoa_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza as OrdinationDirectoryFormat to directory pcoa_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza.txt
Exported dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza as DistanceMatrixDirectoryFormat to directory dm_rpca_s3_min500000_CS_SPF_BioGeo_Shotgun_ASV.qza.txt
Enter filepath of filtered ASV table.qza file s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza
This will perform Robust Aitchison PCA on a count matrix
Saved PCoAResults % Properties('biplot') to: biplot_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza
Saved DistanceMatrix to: dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza
Saved PCoAResults to: pcoa_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza
Exported pcoa_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza as OrdinationDirectoryFormat to directory pcoa_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza.txt
Exported dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza as DistanceMatrixDirectoryFormat to directory dm_rpca_s9_min1000000_UCLA_O_SPF_BioGeo_Shotgun_ASV.qza.txt

```