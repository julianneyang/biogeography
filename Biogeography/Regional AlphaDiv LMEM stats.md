# Site Differences
### Luminal Dataset: Site_General
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site_General 
                                     Value Std.Error  DF    t-value p-value
(Intercept)                       5.046404 0.1869831 230  26.988550  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.618353 0.2463746  43  -2.509807  0.0159
Sequencing_RunNovaSeq_Mar_Twenty  0.589189 0.5633134  43   1.045934  0.3014
SexMale                          -0.055429 0.2238138  43  -0.247659  0.8056
Site_GeneralSI                   -1.439385 0.1433985 230 -10.037660  0.0000

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                                      Value Std.Error  DF    t-value p-value
(Intercept)                       282.99532  11.35413 230  24.924444  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -32.79445  15.03846  43  -2.180706  0.0347
Sequencing_RunNovaSeq_Mar_Twenty    0.86154  34.34707  43   0.025083  0.9801
SexMale                             9.65992  13.66145  43   0.707094  0.4833
Site_GeneralSI                   -115.20840   8.43697 230 -13.655182  0.0000

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site_General 
                                     Value Std.Error  DF    t-value p-value
(Intercept)                       333.1979  13.54867 230  24.592660  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -35.5298  18.12603  43  -1.960151  0.0565
Sequencing_RunNovaSeq_Mar_Twenty   16.6037  41.31333  43   0.401896  0.6898
SexMale                            16.5133  16.46650  43   1.002841  0.3215
Site_GeneralSI                   -117.2762   9.40234 230 -12.473088  0.0000

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.6175785 0.02089234 230 29.560051  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0618424 0.02759381  43 -2.241170  0.0302
Sequencing_RunNovaSeq_Mar_Twenty  0.0800921 0.06305969  43  1.270101  0.2109
SexMale                          -0.0135784 0.02506707  43 -0.541681  0.5908
Site_GeneralSI                   -0.1272948 0.01579750 230 -8.057906  0.0000
```
### Mucosal Dataset: Site_General
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site_General 
                                     Value Std.Error  DF    t-value p-value
(Intercept)                       5.053586 0.3221615 221  15.686497  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.404508 0.4574832  44  -0.884204  0.3814
Sequencing_RunNovaSeq_Mar_Twenty  1.731545 1.0465153  44   1.654582  0.1051
SexMale                          -0.776099 0.4208949  44  -1.843925  0.0719
Site_GeneralSI                   -1.616926 0.1154658 221 -14.003504  0.0000

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                                      Value Std.Error  DF    t-value p-value
(Intercept)                       302.08166  17.76538 221  17.003952  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   -7.20383  25.03831  44  -0.287712  0.7749
Sequencing_RunNovaSeq_Mar_Twenty   73.83750  57.13044  44   1.292437  0.2030
SexMale                           -43.77164  23.02019  44  -1.901445  0.0638
Site_GeneralSI                   -112.96171   7.93804 221 -14.230430  0.0000

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site_General 
                                     Value Std.Error  DF    t-value p-value
(Intercept)                       347.0675  20.63381 221  16.820329  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   10.1339  29.08151  44   0.348466  0.7292
Sequencing_RunNovaSeq_Mar_Twenty   88.4543  66.35625  44   1.333021  0.1894
SexMale                           -44.0905  26.73754  44  -1.649012  0.1063
Site_GeneralSI                   -116.5352   9.21628 221 -12.644491  0.0000

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                                      Value  Std.Error  DF    t-value p-value
(Intercept)                       0.6082636 0.03878529 221  15.682845  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0420446 0.05517480  44  -0.762025  0.4501
Sequencing_RunNovaSeq_Mar_Twenty  0.1999398 0.12629137  44   1.583163  0.1205
SexMale                          -0.0837715 0.05077025  44  -1.650012  0.1061
Site_GeneralSI                   -0.1616101 0.01294773 221 -12.481726  0.0000
```
### Luminal Dataset:Site
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       4.733369 0.2319355 226 20.408128  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.625492 0.2458553  43 -2.544149  0.0146
Sequencing_RunNovaSeq_Mar_Twenty  0.629951 0.5618573  43  1.121194  0.2684
SexMale                          -0.064143 0.2233475  43 -0.287190  0.7753
SiteProximal_Colon                0.312607 0.2400919 226  1.302029  0.1942
SiteCecum                         0.641680 0.2400919 226  2.672642  0.0081
SiteIleum                        -1.334510 0.2415425 226 -5.524948  0.0000
SiteJejunum                      -1.324213 0.2400919 226 -5.515442  0.0000
SiteDuodenum                     -0.678935 0.2446535 226 -2.775088  0.0060

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                       266.12996  13.99030 226 19.022456  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -33.06082  15.01257  43 -2.202209  0.0331
Sequencing_RunNovaSeq_Mar_Twenty    2.85477  34.28117  43  0.083275  0.9340
SexMale                             9.26819  13.63824  43  0.679574  0.5004
SiteProximal_Colon                 22.17021  14.25864 226  1.554862  0.1214
SiteCecum                          29.04255  14.25864 226  2.036839  0.0428
SiteIleum                        -101.64604  14.34512 226 -7.085758  0.0000
SiteJejunum                      -115.10638  14.25864 226 -8.072746  0.0000
SiteDuodenum                      -76.20112  14.53045 226 -5.244235  0.0000

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       321.3243  16.33198 226 19.674553  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -35.9243  18.09548  43 -1.985264  0.0535
Sequencing_RunNovaSeq_Mar_Twenty   19.2442  41.23785  43  0.466664  0.6431
SexMale                            15.9754  16.43906  43  0.971797  0.3366
SiteProximal_Colon                 15.3605  15.90946 226  0.965492  0.3353
SiteCecum                          21.1433  15.90946 226  1.328975  0.1852
SiteIleum                        -112.7987  16.00689 226 -7.046888  0.0000
SiteJejunum                      -124.2390  15.90946 226 -7.809125  0.0000
SiteDuodenum                      -76.4352  16.21531 226 -4.713766  0.0000

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.5852670 0.02576003 226 22.719966  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0627419 0.02753045  43 -2.279001  0.0277
Sequencing_RunNovaSeq_Mar_Twenty  0.0849603 0.06288218  43  1.351103  0.1837
SexMale                          -0.0146369 0.02501013  43 -0.585237  0.5614
SiteProximal_Colon                0.0310319 0.02639238 226  1.175791  0.2409
SiteCecum                         0.0677844 0.02639238 226  2.568333  0.0109
SiteIleum                        -0.1238816 0.02655224 226 -4.665581  0.0000
SiteJejunum                      -0.1141412 0.02639238 226 -4.324779  0.0000
SiteDuodenum                     -0.0420237 0.02689494 226 -1.562514  0.1196
```

###  Mucosal : Site
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       4.686239 0.3421392 217 13.696880  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.426646 0.4572576  44 -0.933054  0.3559
Sequencing_RunNovaSeq_Mar_Twenty  1.734101 1.0462853  44  1.657388  0.1046
SexMale                          -0.780793 0.4207378  44 -1.855772  0.0702
SiteProximal_Colon                0.877405 0.1913606 217  4.585085  0.0000
SiteCecum                         0.195997 0.1913606 217  1.024227  0.3069
SiteIleum                        -1.372335 0.1932855 217 -7.100043  0.0000
SiteJejunum                      -1.258579 0.1958796 217 -6.425270  0.0000
SiteDuodenum                     -1.076349 0.1986645 217 -5.417924  0.0000

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                       305.76575  19.67708 217 15.539185  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   -7.51911  25.09580  44 -0.299616  0.7659
Sequencing_RunNovaSeq_Mar_Twenty   73.89394  57.25353  44  1.290644  0.2036
SexMale                           -44.68579  23.07402  44 -1.936628  0.0592
SiteProximal_Colon                  4.80321  13.75329 217  0.349241  0.7272
SiteCecum                         -13.67595  13.75329 217 -0.994377  0.3211
SiteIleum                        -122.92857  13.89155 217 -8.849161  0.0000
SiteJejunum                      -123.07679  14.07507 217 -8.744311  0.0000
SiteDuodenum                     -100.96530  14.27429 217 -7.073228  0.0000

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       357.4016  22.84377 217 15.645479  0.0000
Sequencing_RunNovaSeq_Jan_Twenty    9.9214  29.14616  44  0.340400  0.7352
Sequencing_RunNovaSeq_Mar_Twenty   88.5049  66.49577  44  1.330985  0.1900
SexMale                           -45.4326  26.79825  44 -1.695356  0.0971
SiteProximal_Colon                 -9.0706  15.94326 217 -0.568928  0.5700
SiteCecum                         -18.2202  15.94326 217 -1.142818  0.2544
SiteIleum                        -134.8745  16.10353 217 -8.375463  0.0000
SiteJejunum                      -136.6876  16.31630 217 -8.377365  0.0000
SiteDuodenum                     -105.0083  16.54725 217 -6.345967  0.0000

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.5568829 0.04078847 217 13.652947  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0448280 0.05509639  44 -0.813628  0.4202
Sequencing_RunNovaSeq_Mar_Twenty  0.2002557 0.12615140  44  1.587423  0.1196
SexMale                          -0.0838919 0.05070435  44 -1.654531  0.1051
SiteProximal_Colon                0.1130778 0.02121323 217  5.330533  0.0000
SiteCecum                         0.0359795 0.02121323 217  1.696088  0.0913
SiteIleum                        -0.1235240 0.02142664 217 -5.764971  0.0000
SiteJejunum                      -0.1076627 0.02171519 217 -4.957943  0.0000
SiteDuodenum                     -0.0955890 0.02202422 217 -4.340175  0.0000
```

### Colon Data
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site + Type 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       5.042106 0.2764287 227 18.240164  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.329861 0.3530214  45 -0.934395  0.3551
Sequencing_RunNovaSeq_Mar_Twenty  1.102511 0.8044915  45  1.370445  0.1773
SexMale                          -0.821586 0.3118602 227 -2.634469  0.0090
SiteProximal_Colon                0.594711 0.1809752 227  3.286148  0.0012
SiteCecum                         0.413226 0.1809752 227  2.283327  0.0233
TypeMucosal                      -0.311706 0.1474337 227 -2.114213  0.0356

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site + Type 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                      289.48216  14.92943 227 19.390039  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   1.78874  18.70470  45  0.095631  0.9242
Sequencing_RunNovaSeq_Mar_Twenty  31.72996  42.57561  45  0.745261  0.4600
SexMale                          -33.17094  16.62874 227 -1.994796  0.0473
SiteProximal_Colon                13.71421  10.38648 227  1.320391  0.1880
SiteCecum                          7.77737  10.38648 227  0.748798  0.4548
TypeMucosal                       -1.24341   8.45819 227 -0.147006  0.8833

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site + Type 
                                    Value Std.Error  DF   t-value p-value
(Intercept)                      342.1036  16.33081 227 20.948359  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   8.6147  20.28014  45  0.424784  0.6730
Sequencing_RunNovaSeq_Mar_Twenty  36.7953  46.13519  45  0.797553  0.4293
SexMale                          -24.3647  18.07276 227 -1.348146  0.1790
SiteProximal_Colon                 3.5008  11.65792 227  0.300296  0.7642
SiteCecum                          1.7388  11.65792 227  0.149152  0.8816
TypeMucosal                       -4.7411   9.49197 227 -0.499490  0.6179

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site + Type 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.6108578 0.03114950 227 19.610517  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0373539 0.03960932  45 -0.943058  0.3507
Sequencing_RunNovaSeq_Mar_Twenty  0.1299892 0.09024155  45  1.440458  0.1567
SexMale                          -0.0902669 0.03504625 227 -2.575650  0.0106
SiteProximal_Colon                0.0720662 0.02068747 227  3.483565  0.0006
SiteCecum                         0.0512940 0.02068747 227  2.479473  0.0139
TypeMucosal                      -0.0409308 0.01685178 227 -2.428873  0.0159
```

Pairwise comparisons using K-W test: **Shannon**
```
> kruskal.test(shannon ~ Type, DC)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 0.013239, df = 1, p-value = 0.9084

> kruskal.test(shannon ~ Type, PC)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 3.8034, df = 1, p-value = 0.05115

> kruskal.test(shannon ~ Type, cec)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 10.245, df = 1, p-value = 0.001371


```
**OTUs**
```R
> kruskal.test(observed_otus ~ Type, DC)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 2.6084, df = 1, p-value = 0.1063

> kruskal.test(observed_otus ~ Type, PC)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 0.11983, df = 1, p-value = 0.7292

> kruskal.test(observed_otus ~ Type, cec)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 4.2985, df = 1, p-value = 0.03814

> 
```
**pielou**
```R
> kruskal.test(pielou_e ~ Type, DC)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 0.23512, df = 1, p-value = 0.6278

> kruskal.test(pielou_e ~ Type, PC)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 4.8546, df = 1, p-value = 0.02757

> kruskal.test(pielou_e ~ Type, cec)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 10.293, df = 1, p-value = 0.001336
```
**chao1**
```R
> kruskal.test(chao1 ~ Type, DC)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 1.3429, df = 1, p-value = 0.2465

> kruskal.test(chao1 ~ Type, PC)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 0.03463, df = 1, p-value = 0.8524

> kruskal.test(chao1 ~ Type, cec)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 1.7359, df = 1, p-value = 0.1877
```
### SI Data
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site + Type 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       3.348770 0.2245483 217 14.913364  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.726573 0.2649479  45 -2.742323  0.0087
Sequencing_RunNovaSeq_Mar_Twenty  1.220062 0.6068008  45  2.010646  0.0504
SexMale                           0.165556 0.2387547 217  0.693417  0.4888
SiteJejunum                       0.051798 0.1820898 217  0.284466  0.7763
SiteDuodenum                      0.428941 0.1853970 217  2.313634  0.0216
TypeMucosal                      -0.508098 0.1508265 217 -3.368761  0.0009

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site + Type 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                      167.90205  13.35238 217 12.574691  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -42.31294  14.24583  45 -2.970198  0.0048
Sequencing_RunNovaSeq_Mar_Twenty  46.88988  32.61987  45  1.437463  0.1575
SexMale                            4.38770  12.95067 217  0.338801  0.7351
SiteJejunum                       -7.21730  12.76958 217 -0.565195  0.5725
SiteDuodenum                      21.66622  12.98576 217  1.668459  0.0967
TypeMucosal                       -1.06499  10.55395 217 -0.100910  0.9197

Fixed effects: chao1 ~ Sequencing_Run + Sex + Site + Type 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                      213.05369  16.31541 217 13.058434  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -33.40954  17.61331  45 -1.896835  0.0643
Sequencing_RunNovaSeq_Mar_Twenty  71.34927  40.32925  45  1.769169  0.0836
SexMale                            5.45210  16.00210 217  0.340711  0.7337
SiteJejunum                       -6.72800  15.36491 217 -0.437881  0.6619
SiteDuodenum                      31.33129  15.62736 217  2.004901  0.0462
TypeMucosal                       -6.48712  12.70203 217 -0.510715  0.6101

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site + Type 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.4529870 0.02767588 217 16.367575  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0707686 0.03391422  45 -2.086693  0.0426
Sequencing_RunNovaSeq_Mar_Twenty  0.1456198 0.07769643  45  1.874215  0.0674
SexMale                           0.0203883 0.03035307 217  0.671706  0.5025
SiteJejunum                       0.0118491 0.02053381 217  0.577056  0.5645
SiteDuodenum                      0.0498862 0.02091545 217  2.385135  0.0179
TypeMucosal                      -0.0763934 0.01702405 217 -4.487383  0.0000
```

**Shannon**
```R
> kruskal.test(shannon ~ Type, duo)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 5.5627, df = 1, p-value = 0.01835

> kruskal.test(shannon ~ Type, jej)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 0.30906, df = 1, p-value = 0.5783

> kruskal.test(shannon ~ Type, ile)

	Kruskal-Wallis rank sum test

data:  shannon by Type
Kruskal-Wallis chi-squared = 0.71126, df = 1, p-value = 0.399
```
**OTUs**
```R
> kruskal.test(observed_otus ~ Type, duo)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 0.046652, df = 1, p-value = 0.829

> kruskal.test(observed_otus ~ Type, jej)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 0.7633, df = 1, p-value = 0.3823

> kruskal.test(observed_otus ~ Type, ile)

	Kruskal-Wallis rank sum test

data:  observed_otus by Type
Kruskal-Wallis chi-squared = 0.079041, df = 1, p-value = 0.7786
```
**Pielou**
```R
> kruskal.test(pielou_e ~ Type, duo)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 9.3533, df = 1, p-value = 0.002226

> kruskal.test(pielou_e ~ Type, jej)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 0.95423, df = 1, p-value = 0.3286

> kruskal.test(pielou_e ~ Type, ile)

	Kruskal-Wallis rank sum test

data:  pielou_e by Type
Kruskal-Wallis chi-squared = 1.3177, df = 1, p-value = 0.251
```
**chao1**
```R
> kruskal.test(chao1 ~ Type, duo)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 0.036125, df = 1, p-value = 0.8493

> kruskal.test(chao1 ~ Type, jej)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 0.52231, df = 1, p-value = 0.4699

> kruskal.test(chao1 ~ Type, ile)

	Kruskal-Wallis rank sum test

data:  chao1 by Type
Kruskal-Wallis chi-squared = 0.10757, df = 1, p-value = 0.7429

```