
## Luminal by Site_General
Nbr of ASVs - does not include DonorID
```R
fixed effects:  observed_features ~ Sequencing_Run + Sex + Site_General 
                           Value Std.Error  DF   t-value p-value
(Intercept)             19.97421  81.58174 214  0.244837  0.8068
Sequencing_RunDec_2017  64.84290  82.09524  38  0.789850  0.4345
Sequencing_RunDec_2018 122.38047  83.20355  38  1.470856  0.1496
Sequencing_RunJan_2017  48.80431  81.81340 214  0.596532  0.5515
Sequencing_RunMay_2018  68.92480  83.83265  38  0.822171  0.4161
Sequencing_RunOct_2017  38.99061  82.58142  38  0.472148  0.6395
SexMale                 -5.18032  16.25892  38 -0.318614  0.7518
Site_GeneralSI          35.51004   9.87850 214  3.594679  0.0004
```
Pielou's evenness
```R
Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Site_General 
                             Value  Std.Error  DF   t-value p-value
(Intercept)             0.29000788 0.12264044 214  2.364700  0.0189
Sequencing_RunDec_2017  0.18028919 0.12391243  38  1.454973  0.1539
Sequencing_RunDec_2018  0.17117184 0.12593911  38  1.359163  0.1821
Sequencing_RunJan_2017  0.16246348 0.12188677 214  1.332905  0.1840
Sequencing_RunMay_2018  0.15908935 0.12747628  38  1.247992  0.2197
Sequencing_RunOct_2017  0.21038462 0.12551517  38  1.676169  0.1019
SexMale                 0.04487501 0.03145530  38  1.426628  0.1619
Site_GeneralSI         -0.07266255 0.01433674 214 -5.068277  0.0000

```

## Mucosal by Site_General
Nbr of ASVs
```R
Fixed effects:  observed_features ~ Sequencing_Run + Sex + Site_General 
                           Value Std.Error  DF    t-value p-value
(Intercept)             67.09738  36.71532 213  1.8275034  0.0690
Sequencing_RunDec_2017  43.03596  40.34632 213  1.0666638  0.2873
Sequencing_RunDec_2018  99.91670  41.80525  40  2.3900515  0.0216
Sequencing_RunJan_2017  31.80046  34.79514 213  0.9139339  0.3618
Sequencing_RunMay_2018 -10.70471  41.50750  40 -0.2578982  0.7978
Sequencing_RunOct_2017  17.83962  42.72362 213  0.4175587  0.6767
SexMale                -35.33555  23.59240  40 -1.4977511  0.1420
Site_GeneralSI          31.02049  12.21250 213  2.5400610  0.0118
```
Pielou's evenness
```R
Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Site_General 
                            Value  Std.Error  DF   t-value p-value
(Intercept)             0.5633732 0.05689646 213  9.901728  0.0000
Sequencing_RunDec_2017 -0.0285178 0.06649865 213 -0.428848  0.6685
Sequencing_RunDec_2018 -0.1912870 0.07003431  40 -2.731333  0.0093
Sequencing_RunJan_2017 -0.0040065 0.04133349 213 -0.096930  0.9229
Sequencing_RunMay_2018 -0.1893276 0.06995049  40 -2.706594  0.0099
Sequencing_RunOct_2017 -0.0288898 0.07310443 213 -0.395185  0.6931
SexMale                -0.0069246 0.04932991  40 -0.140373  0.8891
Site_GeneralSI         -0.0118579 0.01398869 213 -0.847678  0.3976
```

## Luminal by Site
Nbr ASVs
```R
Fixed effects:  observed_features ~ Sequencing_Run + Sex + Site 
                           Value Std.Error  DF   t-value p-value
(Intercept)            -22.94697  80.36925 210 -0.285519  0.7755
Sequencing_RunDec_2017 108.31270  80.34597  38  1.348079  0.1856
Sequencing_RunDec_2018 166.26003  81.44813  38  2.041299  0.0482
Sequencing_RunJan_2017  92.97903  79.99070 210  1.162373  0.2464
Sequencing_RunMay_2018 112.84432  82.09312  38  1.374589  0.1773
Sequencing_RunOct_2017  82.58380  80.88459  38  1.021008  0.3137
SexMale                 -5.83302  16.30169  38 -0.357817  0.7225
SiteProximal_Colon       0.53821  16.55044 210  0.032519  0.9741
SiteCecum               -2.02273  16.34541 210 -0.123749  0.9016
SiteIleum                4.37645  16.44709 210  0.266093  0.7904
SiteJejunum             25.90909  16.34541 210  1.585099  0.1144
SiteDuodenum            76.00387  16.55141 210  4.591988  0.0000
```
Pielou's evenness
```R
Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Site 
                             Value  Std.Error  DF    t-value p-value
(Intercept)             0.26904635 0.12523646 210  2.1483069  0.0328
Sequencing_RunDec_2017  0.19486345 0.12568118  38  1.5504584  0.1293
Sequencing_RunDec_2018  0.18593190 0.12771009  38  1.4558904  0.1536
Sequencing_RunJan_2017  0.17727407 0.12374653 210  1.4325580  0.1535
Sequencing_RunMay_2018  0.17388269 0.12922534  38  1.3455773  0.1864
Sequencing_RunOct_2017  0.22493283 0.12727430  38  1.7673075  0.0852
SexMale                 0.04465617 0.03144255  38  1.4202465  0.1637
SiteProximal_Colon      0.00928312 0.02504847 210  0.3706064  0.7113
SiteCecum               0.01002465 0.02472382 210  0.4054653  0.6855
SiteIleum              -0.07652186 0.02488583 210 -3.0749171  0.0024
SiteJejunum            -0.06876464 0.02472382 210 -2.7813114  0.0059
SiteDuodenum           -0.05316006 0.02505109 210 -2.1220659  0.0350
```

## mucosal by Site
Number of ASVs
```R
Fixed effects:  observed_features ~ Sequencing_Run + Sex + Site 
                           Value Std.Error  DF    t-value p-value
(Intercept)             75.78724  38.20299 209  1.9838039  0.0486
Sequencing_RunDec_2017  52.53766  40.35561 209  1.3018675  0.1944
Sequencing_RunDec_2018 110.28068  41.83169  40  2.6362949  0.0119
Sequencing_RunJan_2017  43.14471  34.59298 209  1.2472098  0.2137
Sequencing_RunMay_2018  -0.39098  41.54970  40 -0.0094099  0.9925
Sequencing_RunOct_2017  28.02673  42.85635 209  0.6539692  0.5139
SexMale                -36.23761  23.80071  40 -1.5225434  0.1357
SiteProximal_Colon     -27.22273  20.48077 209 -1.3291847  0.1852
SiteCecum              -28.53349  20.34861 209 -1.4022327  0.1623
SiteIleum              -13.43018  20.25416 209 -0.6630825  0.5080
SiteJejunum             12.76365  20.62756 209  0.6187670  0.5367
SiteDuodenum            42.45850  20.82792 209  2.0385378  0.0428
```
Evenness
```R
Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Site 
                            Value  Std.Error  DF   t-value p-value
(Intercept)             0.5547375 0.05838006 209  9.502173  0.0000
Sequencing_RunDec_2017 -0.0211139 0.06659039 209 -0.317071  0.7515
Sequencing_RunDec_2018 -0.1835697 0.07014870  40 -2.616865  0.0125
Sequencing_RunJan_2017  0.0044852 0.04156049 209  0.107921  0.9142
Sequencing_RunMay_2018 -0.1810231 0.07007327  40 -2.583340  0.0135
Sequencing_RunOct_2017 -0.0224226 0.07332759 209 -0.305787  0.7601
SexMale                -0.0083412 0.04937325  40 -0.168942  0.8667
SiteProximal_Colon      0.0107827 0.02371063 209  0.454764  0.6498
SiteCecum              -0.0054710 0.02354963 209 -0.232319  0.8165
SiteIleum              -0.0310361 0.02345393 209 -1.323279  0.1872
SiteJejunum            -0.0145871 0.02393221 209 -0.609516  0.5428
SiteDuodenum            0.0185102 0.02417888 209  0.765553  0.**4448**
```



## Type Differences 

### Duodenum
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
       AIC      BIC    logLik
  1039.101 1064.882 -508.5503

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev: 0.002878991

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept) Residual
StdDev: 0.009752592 156.3266

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF    t-value p-value
(Intercept)            103.08547  75.14874 38  1.3717524  0.1782
Sequencing_RunDec_2017  10.86654  84.11617 38  0.1291849  0.8979
Sequencing_RunDec_2018 235.56510  93.02838 10  2.5321853  0.0298
Sequencing_RunJan_2017   0.46324  83.89737 38  0.0055215  0.9956
Sequencing_RunMay_2018 -14.73722  94.97233 38 -0.1551738  0.8775
Sequencing_RunOct_2017  -2.00705  86.24524 10 -0.0232714  0.9819
SexMale                  0.66422  50.30906 10  0.0132027  0.9897
TypeMucosal            -20.35683  34.45401 38 -0.5908407  0.5581

```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
        AIC      BIC  logLik
  -1.949209 23.83265 11.9746

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:  0.06300892

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept)  Residual
StdDev: 3.281941e-06 0.1761701

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.4415511 0.09747211 38  4.530025  0.0001
Sequencing_RunDec_2017 -0.0642318 0.11091079 38 -0.579131  0.5659
Sequencing_RunDec_2018 -0.0653836 0.12419329 10 -0.526467  0.6100
Sequencing_RunJan_2017 -0.0215640 0.10128009 38 -0.212914  0.8325
Sequencing_RunMay_2018 -0.1320459 0.12537527 38 -1.053206  0.2989
Sequencing_RunOct_2017 -0.1040519 0.11733920 10 -0.886762  0.3960
SexMale                 0.0895394 0.07427804 10  1.205462  0.2558
TypeMucosal             0.0517235 0.03897690 38  1.327029  0.1924
```

### Jejunum
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
      AIC      BIC  logLik
  1001.98 1028.044 -489.99

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:    31.06065

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept) Residual
StdDev:     30.5155 98.48452

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF    t-value p-value
(Intercept)             96.21736  61.82539 38  1.5562760  0.1279
Sequencing_RunDec_2017   0.48517  64.85095 38  0.0074813  0.9941
Sequencing_RunDec_2018 118.77089  69.41243 11  1.7110897  0.1151
Sequencing_RunJan_2017  -7.03873  59.53459 38 -0.1182292  0.9065
Sequencing_RunMay_2018  -2.28987  69.06090 38 -0.0331572  0.9737
Sequencing_RunOct_2017  22.75314  70.34684 38  0.3234423  0.7481
SexMale                -18.86463  39.45646 11 -0.4781126  0.6419
TypeMucosal             -0.78210  22.01111 38 -0.0355320  0.9718
```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
       AIC      BIC  logLik
  -24.0612 2.002728 23.0306

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:   0.1068019

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept)  Residual
StdDev: 2.282611e-06 0.1471606

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.5533505 0.11061860 38  5.002327  0.0000
Sequencing_RunDec_2017 -0.2185937 0.11926219 38 -1.832883  0.0747
Sequencing_RunDec_2018 -0.0847863 0.13636373 11 -0.621766  0.5468
Sequencing_RunJan_2017 -0.1278635 0.09114471 38 -1.402862  0.1688
Sequencing_RunMay_2018 -0.1024419 0.12976447 38 -0.789445  0.4348
Sequencing_RunOct_2017 -0.0888794 0.13372916 38 -0.664623  0.5103
SexMale                -0.0139790 0.08686056 11 -0.160936  0.8751
TypeMucosal             0.0144256 0.03303039 38  0.436737  0.6648
```

### Ileum
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
       AIC      BIC    logLik
  875.2089 901.4112 -426.6044

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev: 0.008377141

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept) Residual
StdDev: 0.003616319 44.12138

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF   t-value p-value
(Intercept)            106.55173  33.34953 40  3.195000  0.0027
Sequencing_RunDec_2017  -4.32552  33.48848 40 -0.129165  0.8979
Sequencing_RunDec_2018  20.09720  33.51515 10  0.599645  0.5621
Sequencing_RunJan_2017  -1.69891  32.94487 40 -0.051568  0.9591
Sequencing_RunMay_2018 -25.04529  33.89090 40 -0.738997  0.4642
Sequencing_RunOct_2017 -22.62611  34.83195 10 -0.649579  0.5306
SexMale                -15.70954  12.86563 10 -1.221047  0.2501
TypeMucosal             -2.69696   9.57643 40 -0.281625  0.7797
```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
        AIC       BIC   logLik
  -28.19956 -1.997268 25.09978

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:  0.04977453

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept)  Residual
StdDev: 2.206543e-06 0.1518896

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.5353435 0.11986863 40  4.466085  0.0001
Sequencing_RunDec_2017 -0.1223021 0.12256776 40 -0.997833  0.3244
Sequencing_RunDec_2018 -0.0928970 0.12582149 10 -0.738324  0.4773
Sequencing_RunJan_2017 -0.0846047 0.11424587 40 -0.740550  0.4633
Sequencing_RunMay_2018 -0.1198129 0.12600447 40 -0.950862  0.3474
Sequencing_RunOct_2017 -0.0846471 0.12963090 10 -0.652986  0.5285
SexMale                -0.0446751 0.05799683 10 -0.770302  0.4589
TypeMucosal             0.0215814 0.03302273 40  0.653532  0.5172
```

### Cecum
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
      AIC      BIC    logLik
  727.543 751.4875 -353.7715

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:    33.47139

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept) Residual
StdDev: 0.0005670486 14.02367

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF   t-value p-value
(Intercept)            108.33780  21.36951 42  5.069737  0.0000
Sequencing_RunDec_2018  14.64682  25.90341  9  0.565440  0.5856
Sequencing_RunJan_2017  -5.92800  25.54052  9 -0.232102  0.8217
Sequencing_RunMay_2018  26.08236  11.54933 42  2.258344  0.0292
Sequencing_RunOct_2017 -32.41547  32.13064  9 -1.008865  0.3394
SexMale                -38.37848  24.15844  9 -1.588616  0.1466
TypeMucosal            -11.34465   3.03559 42 -3.737209  0.0006
```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
        AIC       BIC   logLik
  -109.6442 -85.69972 64.82211

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:   0.1481758

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept)   Residual
StdDev: 2.73481e-07 0.08228156

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.4866274 0.09614676 42  5.061298  0.0000
Sequencing_RunDec_2018  0.0017908 0.11746613  9  0.015246  0.9882
Sequencing_RunJan_2017  0.0593720 0.11486907  9  0.516867  0.6177
Sequencing_RunMay_2018  0.1574549 0.06492035 42  2.425355  0.0197
Sequencing_RunOct_2017  0.0761746 0.14406248  9  0.528761  0.6098
SexMale                -0.0657351 0.10848454  9 -0.605940  0.5595
TypeMucosal            -0.0423317 0.01778895 42 -2.379664  0.0219
```

### Proximal Colon
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
       AIC      BIC    logLik
  711.2895 734.8566 -345.6447

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:    28.22122

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept) Residual
StdDev: 0.000668914 15.23927

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF   t-value p-value
(Intercept)            113.88954 18.277692 39  6.231068  0.0000
Sequencing_RunDec_2018   9.31841 22.336994  9  0.417174  0.6863
Sequencing_RunJan_2017  -9.55964 21.839133  9 -0.437730  0.6719
Sequencing_RunMay_2018  13.79347 12.098255 39  1.140121  0.2612
Sequencing_RunOct_2017 -39.83870 27.395319  9 -1.454216  0.1799
SexMale                -40.22899 20.635880  9 -1.949468  0.0830
TypeMucosal             -8.26834  3.357355 39 -2.462755  0.0183
```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
        AIC       BIC   logLik
  -138.1848 -114.6177 79.09242

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:    0.124191

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept)   Residual
StdDev: 8.841917e-07 0.06562536

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.5427145 0.08031452 39  6.757365  0.0000
Sequencing_RunDec_2018 -0.0794971 0.09808185  9 -0.810518  0.4385
Sequencing_RunJan_2017 -0.0178559 0.09596808  9 -0.186061  0.8565
Sequencing_RunMay_2018  0.0194972 0.05229597 39  0.372824  0.7113
Sequencing_RunOct_2017 -0.0093226 0.12041615  9 -0.077420  0.9400
SexMale                -0.0559870 0.09069138  9 -0.617335  0.5523
TypeMucosal            -0.0101906 0.01445948 39 -0.704769  0.4851
```

### Distal Colon
```R
> output=lme(fixed= observed_features ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
       AIC      BIC    logLik
  979.2825 1005.485 -478.6412

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:    23.71723

 Formula: ~1 | MouseID %in% Donor_ID
        (Intercept) Residual
StdDev:     0.01062 83.17638

Fixed effects:  observed_features ~ Sequencing_Run + Sex + Type 
                           Value Std.Error DF    t-value p-value
(Intercept)            104.76012  87.70769 40  1.1944234  0.2393
Sequencing_RunDec_2017   8.93593  89.96430 40  0.0993275  0.9214
Sequencing_RunDec_2018  23.16049  92.89352 10  0.2493229  0.8082
Sequencing_RunJan_2017  -6.95638  87.58158 40 -0.0794274  0.9371
Sequencing_RunMay_2018 -25.84971  92.48642 40 -0.2794973  0.7813
Sequencing_RunOct_2017 -21.07569  91.14527 10 -0.2312319  0.8218
SexMale                -24.80457  30.20182 10 -0.8212941  0.4306
TypeMucosal             15.13115  17.92187 40  0.8442840  0.4035
```

```R
> output=lme(fixed= pielou_evenness ~ Sequencing_Run + Sex + Type, random = ~1|(Donor_ID/MouseID), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
        AIC       BIC   logLik
  -102.8533 -76.65098 62.42664

Random effects:
 Formula: ~1 | Donor_ID
        (Intercept)
StdDev:   0.1101631

 Formula: ~1 | MouseID %in% Donor_ID
         (Intercept)   Residual
StdDev: 1.356651e-07 0.08672832

Fixed effects:  pielou_evenness ~ Sequencing_Run + Sex + Type 
                            Value  Std.Error DF   t-value p-value
(Intercept)             0.6043871 0.11542420 40  5.236225  0.0000
Sequencing_RunDec_2017 -0.0389993 0.12621632 40 -0.308988  0.7589
Sequencing_RunDec_2018 -0.2049415 0.14158092 10 -1.447522  0.1784
Sequencing_RunJan_2017 -0.0574601 0.09376089 40 -0.612837  0.5435
Sequencing_RunMay_2018 -0.1040569 0.13112384 40 -0.793577  0.4321
Sequencing_RunOct_2017 -0.0892985 0.14049261 10 -0.635610  0.5393
SexMale                -0.0298101 0.08314891 10 -0.358515  0.7274
TypeMucosal            -0.0230002 0.01883646 40 -1.221046  0.2292
```