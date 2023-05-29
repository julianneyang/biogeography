# Luminal
## Site General
observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             180.79225  12.86174 42 14.056589  0.0000
Sequencing_Run2014_Sept -14.66527  19.07268 42 -0.768915  0.4462
Sequencing_Run2015_Sept   6.29435  15.76558  6  0.399247  0.7035
SexMale                  12.02778  13.17383  6  0.913005  0.3964
Site_GeneralSI          -23.25163  12.66587 42 -1.835771  0.0735
```
pielou_e
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6203418 0.03400573 42 18.242272  0.0000
Sequencing_Run2014_Sept -0.0220456 0.03903687 42 -0.564738  0.5753
Sequencing_Run2015_Sept  0.0223811 0.04567965  6  0.489958  0.6416
SexMale                 -0.0084451 0.03921162  6 -0.215373  0.8366
Site_GeneralSI          -0.0852408 0.02593724 42 -3.286425  0.0021
```
## Site
observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             173.64949  35.62651 38  4.874165  0.0000
Sequencing_Run2014_Sept  -7.52250  38.45659 38 -0.195610  0.8460
Sequencing_Run2015_Sept   8.32273  17.30800  6  0.480861  0.6476
SexMale                  12.02778  13.54449  6  0.888020  0.4087
SiteProximal_Colon        1.92694  36.24873 38  0.053159  0.9579
SiteCecum                12.59361  36.24873 38  0.347422  0.7302
SiteIleum                -8.08219  37.18737 38 -0.217337  0.8291
SiteJejunum             -28.62861  36.24873 38 -0.789783  0.4346
SiteDuodenum            -11.85084  36.24873 38 -0.326931  0.7455
```
pielou_e
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6692383 0.07580387 38  8.828551  0.0000
Sequencing_Run2014_Sept -0.0709421 0.07864043 38 -0.902107  0.3727
Sequencing_Run2015_Sept  0.0116421 0.04735144  6  0.245867  0.8140
SexMale                 -0.0084451 0.03915113  6 -0.215706  0.8364
SiteProximal_Colon      -0.0512558 0.07411548 38 -0.691567  0.4934
SiteCecum               -0.0502437 0.07411548 38 -0.677911  0.5019
SiteIleum               -0.1529612 0.07613602 38 -2.009052  0.0517
SiteJejunum             -0.1409148 0.07411548 38 -1.901286  0.0649
SiteDuodenum            -0.1048293 0.07411548 38 -1.414405  0.1654
```

# Mucosal
## Site General
observed otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             215.47406  18.96005 38 11.364633  0.0000
Sequencing_Run2014_Sept  -3.14475  15.98815 38 -0.196693  0.8451
Sequencing_Run2015_Sept  -7.39611  22.22479  6 -0.332787  0.7506
SexMale                  14.83302  15.48443  6  0.957932  0.3751
Site_GeneralSI          -43.93322  13.62934 38 -3.223429  0.0026
```

pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6732106 0.04657764 38 14.453514  0.0000
Sequencing_Run2014_Sept -0.0520112 0.03389552 38 -1.534458  0.1332
Sequencing_Run2015_Sept  0.0135494 0.05920621  6  0.228851  0.8266
SexMale                 -0.0356734 0.04556982  6 -0.782830  0.4635
Site_GeneralSI           0.0230907 0.02778495 38  0.831050  0.4111
```

## Site
observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             222.44980  21.32636 34 10.430743  0.0000
Sequencing_Run2014_Sept  -0.42266  16.93447 34 -0.024959  0.9802
Sequencing_Run2015_Sept  -4.52418  22.75636  6 -0.198809  0.8490
SexMale                  14.63454  15.41313  6  0.949485  0.3790
SiteProximal_Colon      -15.53871  20.15810 34 -0.770842  0.4461
SiteCecum               -11.53871  20.15810 34 -0.572411  0.5708
SiteIleum               -49.88906  20.84925 34 -2.392847  0.0224
SiteJejunum             -56.81689  21.39133 34 -2.656072  0.0119
SiteDuodenum            -48.88765  21.95704 34 -2.226514  0.0327
```

pielou e 
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6686756 0.04759529 34 14.049196  0.0000
Sequencing_Run2014_Sept -0.0356114 0.03213729 34 -1.108103  0.2756
Sequencing_Run2015_Sept  0.0269919 0.05948476  6  0.453761  0.6660
SexMale                 -0.0381581 0.04648897  6 -0.820798  0.4431
SiteProximal_Colon      -0.0258341 0.03402253 34 -0.759324  0.4529
SiteCecum                0.0011330 0.03402253 34  0.033303  0.9736
SiteIleum               -0.0441229 0.03564322 34 -1.237906  0.2242
SiteJejunum              0.0601304 0.03683963 34  1.632221  0.1119
SiteDuodenum             0.0585842 0.03799733 34  1.541797  0.1324
```