# Luminal
## Site General
otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             102.29459  19.52396 42  5.239439  0.0000
Sequencing_Run2014_Sept   2.72922  19.86359 42  0.137398  0.8914
Sequencing_Run2015_Sept -49.55142  26.91539  6 -1.841007  0.1152
SexMale                  20.45833  23.26321  6  0.879429  0.4130
Site_GeneralSI           79.16774  13.19968 42  5.997703  0.0000
```
pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.4322125 0.04481554 42  9.644255  0.0000
Sequencing_Run2014_Sept -0.0121867 0.05281213 42 -0.230755  0.8186
Sequencing_Run2015_Sept -0.0366942 0.05979621  6 -0.613654  0.5620
SexMale                  0.1294999 0.05123528  6  2.527553  0.0448
Site_GeneralSI          -0.0724278 0.03508869 42 -2.064135  0.0452
```

## Site
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                            Value Std.Error DF    t-value p-value
(Intercept)              88.65782  40.15714 38  2.2077722  0.0334
Sequencing_Run2014_Sept  16.36599  40.68774 38  0.4022340  0.6898
Sequencing_Run2015_Sept -46.61615  27.70791  6 -1.6824129  0.1435
SexMale                  20.45833  23.24533  6  0.8801050  0.4127
SiteProximal_Colon       18.17355  38.34524 38  0.4739453  0.6383
SiteCecum                10.17355  38.34524 38  0.2653145  0.7922
SiteIleum                86.17355  38.34524 38  2.2473074  0.0305
SiteJejunum              93.95133  38.34524 38  2.4501429  0.0190
SiteDuodenum             97.21513  39.40409 38  2.4671332  0.0182
```
pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.5323831 0.09734246 38  5.469176  0.0000
Sequencing_Run2014_Sept -0.1123573 0.10054929 38 -1.117435  0.2708
Sequencing_Run2015_Sept -0.0478205 0.06204273  6 -0.770768  0.4701
SexMale                  0.1294999 0.05146303  6  2.516367  0.0455
SiteProximal_Colon      -0.1079174 0.09476298 38 -1.138814  0.2619
SiteCecum               -0.1072664 0.09476298 38 -1.131944  0.2648
SiteIleum               -0.1875702 0.09476298 38 -1.979362  0.0551
SiteJejunum             -0.2175469 0.09476298 38 -2.295695  0.0273
SiteDuodenum            -0.0978356 0.09735379 38 -1.004949  0.3213
```

# Mucosal

## Site General 
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             185.62288  19.79252 38  9.378436  0.0000
Sequencing_Run2014_Sept -46.09069  17.27179 38 -2.668554  0.0111
Sequencing_Run2015_Sept -31.51980  22.87386  6 -1.377983  0.2174
SexMale                 -26.11111  14.53721  6 -1.796157  0.1226
Site_GeneralSI           27.76405  15.09807 38  1.838913  0.0738
```

```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6007546 0.04294577 38 13.988681  0.0000
Sequencing_Run2014_Sept -0.1355998 0.03160947 38 -4.289847  0.0001
Sequencing_Run2015_Sept -0.1148688 0.05501898  6 -2.087803  0.0818
SexMale                  0.0781333 0.04105818  6  1.902990  0.1057
Site_GeneralSI          -0.0009532 0.02562360 38 -0.037200  0.9705
```

## Site
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             204.04601  21.99133 34  9.278473  0.0000
Sequencing_Run2014_Sept -34.97552  16.89742 34 -2.069874  0.0461
Sequencing_Run2015_Sept -14.36938  22.37020  6 -0.642345  0.5444
SexMale                 -26.11111  13.71763  6 -1.903471  0.1057
SiteProximal_Colon      -31.95181  23.25478 34 -1.373989  0.1784
SiteCecum               -49.90889  22.82795 34 -2.186306  0.0358
SiteIleum                30.56349  23.77048 34  1.285775  0.2072
SiteJejunum             -16.34245  22.88173 34 -0.714214  0.4800
SiteDuodenum             -4.06763  22.99060 34 -0.176926  0.8606
```

```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6161133 0.04505804 34 13.673771  0.0000
Sequencing_Run2014_Sept -0.1523451 0.03155280 34 -4.828259  0.0000
Sequencing_Run2015_Sept -0.1236119 0.05370479  6 -2.301692  0.0610
SexMale                  0.0781333 0.03910952  6  1.997808  0.0927
SiteProximal_Colon      -0.0059676 0.03917494 34 -0.152333  0.8798
SiteCecum                0.0027277 0.03847152 34  0.070902  0.9439
SiteIleum               -0.0539714 0.04025179 34 -1.340844  0.1889
SiteJejunum             -0.0217369 0.03856932 34 -0.563581  0.5767
SiteDuodenum             0.0509176 0.03863495 34  1.317917  0.1963
```