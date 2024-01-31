
## Luminal by Site_General
Nbr of ASVs
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