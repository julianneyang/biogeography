# Full Dataset
### Differences in Microbiota* Site_General
Shannon
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept)  Residual
StdDev:     0.45863 0.7784938

Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Microbiota * Site_General 
                              Value  Std.Error  DF   t-value p-value
(Intercept)                4.023492 0.14001777 181 28.735578  0.0000
Sequencing_Run1            0.276826 0.12898906 181  2.146121  0.0332
Sequencing_Run2           -0.246880 0.13641552 181 -1.809764  0.0720
Sex1                      -0.154955 0.13802893  15 -1.122628  0.2792
Type1                     -0.404185 0.06315725 181 -6.399668  0.0000
Microbiota1                0.441448 0.12126926  15  3.640228  0.0024
Site_General1              0.098634 0.06123548 181  1.610731  0.1090
Microbiota1:Site_General1  0.150227 0.05466987 181  2.747894  0.0066
```
Observed_OTUs
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:     17.2816 43.67883

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Microbiota * Site_General 
                              Value Std.Error  DF   t-value p-value
(Intercept)               164.91611  5.930818 181 27.806638  0.0000
Sequencing_Run1             9.55628  6.080884 181  1.571528  0.1178
Sequencing_Run2            -1.34594  6.563570 181 -0.205062  0.8378
Sex1                       -2.60247  5.785256  15 -0.449844  0.6593
Type1                      -9.40300  3.536412 181 -2.658909  0.0085
Microbiota1                20.33345  5.105612  15  3.982568  0.0012
Site_General1              -6.21665  3.429433 181 -1.812735  0.0715
Microbiota1:Site_General1  23.54026  3.066018 181  7.677795  0.0000
```
Chao1
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:     23.2327  52.2999

Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Microbiota * Site_General 
                              Value Std.Error  DF   t-value p-value
(Intercept)               191.14728  7.653910 181 24.973808  0.0000
Sequencing_Run1            10.37915  7.601402 181  1.365425  0.1738
Sequencing_Run2            11.36404  8.159565 181  1.392726  0.1654
Sex1                       -2.64178  7.491412  15 -0.352641  0.7293
Type1                       0.02632  4.237174 181  0.006211  0.9951
Microbiota1                21.79293  6.602008  15  3.300955  0.0048
Site_General1              -8.59690  4.108760 181 -2.092335  0.0378
Microbiota1:Site_General1  24.97399  3.671677 181  6.801795  0.0000
```
Pielou E
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.05606691 0.09374986

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Microbiota * Site_General 
                               Value   Std.Error  DF  t-value p-value
(Intercept)                0.5502095 0.017061876 181 32.24789  0.0000
Sequencing_Run1            0.0289526 0.015658990 181  1.84894  0.0661
Sequencing_Run2           -0.0366302 0.016546974 181 -2.21371  0.0281
Sex1                      -0.0208950 0.016825066  15 -1.24190  0.2333
Type1                     -0.0494430 0.007606165 181 -6.50038  0.0000
Microbiota1                0.0440918 0.014780011  15  2.98321  0.0093
Site_General1              0.0190836 0.007374678 181  2.58772  0.0104
Microbiota1:Site_General1  0.0030605 0.006583692 181  0.46485  0.6426
```
### Differences in Microbiota* Site
Shannon
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept)  Residual
StdDev:   0.4532577 0.7593602

Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Microbiota * Site 
                                           Value Std.Error  DF   t-value p-value
(Intercept)                             4.451544 0.2971240 174 14.982107  0.0000
Sequencing_Run2014_Sept                -0.664722 0.1627714 174 -4.083777  0.0001
Sequencing_Run2015_Sept                -0.322666 0.3206464  14 -1.006297  0.3313
SexMale                                 0.306311 0.2721456  14  1.125540  0.2793
TypeMucosal                             0.879291 0.1255400 174  7.004069  0.0000
MicrobiotaHumanized                    -1.365668 0.3312456  14 -4.122824  0.0010
SiteDistal_Colon                        0.207460 0.2620829 174  0.791583  0.4297
SiteDuodenum                           -0.302045 0.2649148 174 -1.140160  0.2558
SiteIleum                              -0.798201 0.2643747 174 -3.019204  0.0029
SiteJejunum                            -0.485672 0.2610772 174 -1.860263  0.0645
SiteProximal_Colon                     -0.138769 0.2531201 174 -0.548233  0.5842
MicrobiotaHumanized:SiteDistal_Colon    0.382741 0.3667992 174  1.043462  0.2982
MicrobiotaHumanized:SiteDuodenum        0.987541 0.3707910 174  2.663336  0.0085
MicrobiotaHumanized:SiteIleum           0.935633 0.3667914 174  2.550858  0.0116
MicrobiotaHumanized:SiteJejunum         0.495522 0.3609394 174  1.372868  0.1716
MicrobiotaHumanized:SiteProximal_Colon  0.207753 0.3608919 174  0.575665  0.5656
```
OTUs
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    16.76455 43.84149

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Microbiota * Site 
                                            Value Std.Error  DF   t-value p-value
(Intercept)                             203.67597 14.296879 174 14.246184  0.0000
Sequencing_Run2014_Sept                 -15.09532  9.315727 174 -1.620413  0.1070
Sequencing_Run2015_Sept                 -18.04647 13.812183  14 -1.306562  0.2124
SexMale                                   5.12894 11.362604  14  0.451388  0.6586
TypeMucosal                              21.01537  7.228535 174  2.907279  0.0041
MicrobiotaHumanized                    -103.50000 16.613837  14 -6.229747  0.0000
SiteDistal_Colon                          0.03465 15.123708 174  0.002291  0.9982
SiteDuodenum                            -34.07625 15.282735 174 -2.229722  0.0270
SiteIleum                               -33.75587 15.256798 174 -2.212513  0.0282
SiteJejunum                             -47.42861 15.065616 174 -3.148136  0.0019
SiteProximal_Colon                       -7.33333 14.613830 174 -0.501808  0.6164
MicrobiotaHumanized:SiteDistal_Colon     30.60716 21.169072 174  1.445843  0.1500
MicrobiotaHumanized:SiteDuodenum        103.20504 21.385255 174  4.825991  0.0000
MicrobiotaHumanized:SiteIleum           116.53259 21.167892 174  5.505158  0.0000
MicrobiotaHumanized:SiteJejunum         110.51853 20.836171 174  5.304167  0.0000
MicrobiotaHumanized:SiteProximal_Colon   18.84398 20.832493 174  0.904547  0.3670
 Correlation: 
```

Chao1
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    22.86343 52.67042

Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Microbiota * Site 
                                            Value Std.Error  DF   t-value p-value
(Intercept)                             239.74829 17.919395 174 13.379262  0.0000
Sequencing_Run2014_Sept                  -2.15525 11.223395 174 -0.192032  0.8479
Sequencing_Run2015_Sept                 -32.19827 17.874546  14 -1.801347  0.0932
SexMale                                   5.20463 14.849220  14  0.350499  0.7312
TypeMucosal                               1.68957  8.691717 174  0.194389  0.8461
MicrobiotaHumanized                    -106.30716 20.601093  14 -5.160268  0.0001
SiteDistal_Colon                          1.31280 18.172277 174  0.072242  0.9425
SiteDuodenum                            -30.92886 18.365034 174 -1.684117  0.0940
SiteIleum                               -28.40824 18.331858 174 -1.549665  0.1230
SiteJejunum                             -47.26677 18.102479 174 -2.611066  0.0098
SiteProximal_Colon                       -6.88493 17.556806 174 -0.392152  0.6954
MicrobiotaHumanized:SiteDistal_Colon     22.63251 25.435185 174  0.889811  0.3748
MicrobiotaHumanized:SiteDuodenum        101.45537 25.700151 174  3.947656  0.0001
MicrobiotaHumanized:SiteIleum           121.98625 25.433982 174  4.796192  0.0000
MicrobiotaHumanized:SiteJejunum         114.55869 25.033197 174  4.576271  0.0000
MicrobiotaHumanized:SiteProximal_Colon   17.31636 25.029097 174  0.691849  0.4900
 
```
Pielou
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.05618543 0.09034422

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Microbiota * Site 
                                            Value  Std.Error  DF   t-value p-value
(Intercept)                             0.5752477 0.03612426 174 15.924139  0.0000
Sequencing_Run2014_Sept                -0.0820637 0.01937773 174 -4.234951  0.0000
Sequencing_Run2015_Sept                -0.0231199 0.03935459  14 -0.587477  0.5662
SexMale                                 0.0413497 0.03348014  14  1.235051  0.2371
TypeMucosal                             0.1070965 0.01493894 174  7.168953  0.0000
MicrobiotaHumanized                    -0.1064144 0.04010498  14 -2.653397  0.0189
SiteDistal_Colon                        0.0287265 0.03118222 174  0.921245  0.3582
SiteDuodenum                           -0.0183093 0.03151982 174 -0.580881  0.5621
SiteIleum                              -0.0865581 0.03145478 174 -2.751828  0.0066
SiteJejunum                            -0.0327431 0.03106259 174 -1.054099  0.2933
SiteProximal_Colon                     -0.0139896 0.03011474 174 -0.464544  0.6428
MicrobiotaHumanized:SiteDistal_Colon    0.0287926 0.04364084 174  0.659763  0.5103
MicrobiotaHumanized:SiteDuodenum        0.0613370 0.04411806 174  1.390293  0.1662
MicrobiotaHumanized:SiteIleum           0.0404388 0.04364006 174  0.926645  0.3554
MicrobiotaHumanized:SiteJejunum        -0.0190114 0.04294285 174 -0.442715  0.6585
MicrobiotaHumanized:SiteProximal_Colon  0.0119585 0.04293736 174  0.278511  0.7810
```

### Microbiota* Type
Shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Site + Microbiota * Type 
                                    Value Std.Error  DF   t-value p-value
(Intercept)                      4.636547 0.2960312 178 15.662359  0.0000
Sequencing_Run2014_Sept         -0.663055 0.1656579 178 -4.002556  0.0001
Sequencing_Run2015_Sept         -0.318677 0.3230189  14 -0.986559  0.3406
SexMale                          0.303173 0.2739986  14  1.106475  0.2872
SiteProximal_Colon              -0.431432 0.1954125 178 -2.207800  0.0285
SiteCecum                       -0.400455 0.1948691 178 -2.054996  0.0413
SiteIleum                       -0.722783 0.2142748 178 -3.373158  0.0009
SiteJejunum                     -0.637948 0.2115498 178 -3.015594  0.0029
SiteDuodenum                    -0.211325 0.2129594 178 -0.992325  0.3224
MicrobiotaHumanized             -0.937586 0.2620216  14 -3.578278  0.0030
TypeMucosal                      0.822716 0.1661999 178  4.950160  0.0000
MicrobiotaHumanized:TypeMucosal  0.121025 0.2180486 178  0.555036  0.5796
```
OTUs
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site + Microbiota * Type 
                                     Value  Std.Error  DF   t-value p-value
(Intercept)                      0.6159595 0.03560794 178 17.298379  0.0000
Sequencing_Run2014_Sept         -0.0826440 0.01931094 178 -4.279648  0.0000
Sequencing_Run2015_Sept         -0.0241875 0.03935446  14 -0.614606  0.5487
SexMale                          0.0411537 0.03349570  14  1.228627  0.2395
SiteProximal_Colon              -0.0511177 0.02276048 178 -2.245899  0.0259
SiteCecum                       -0.0433397 0.02269747 178 -1.909450  0.0578
SiteIleum                       -0.1094104 0.02496088 178 -4.383273  0.0000
SiteJejunum                     -0.0865233 0.02464413 178 -3.510909  0.0006
SiteDuodenum                    -0.0317339 0.02480860 178 -1.279150  0.2025
MicrobiotaHumanized             -0.1001921 0.03178393  14 -3.152287  0.0071
TypeMucosal                      0.0934991 0.01935970 178  4.829571  0.0000
MicrobiotaHumanized:TypeMucosal  0.0276359 0.02539809 178  1.088109  0.2780
```

Chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Site + Microbiota * Type 
                                    Value Std.Error  DF   t-value p-value
(Intercept)                     214.28628  18.68000 178 11.471429  0.0000
Sequencing_Run2014_Sept          -1.16226  12.34880 178 -0.094120  0.9251
Sequencing_Run2015_Sept         -29.86969  18.54464  14 -1.610691  0.1296
SexMale                           4.80059  15.29038  14  0.313961  0.7582
SiteProximal_Colon              -10.54835  14.67333 178 -0.718880  0.4732
SiteCecum                       -12.74106  14.63102 178 -0.870826  0.3850
SiteIleum                        20.94776  16.07050 178  1.303491  0.1941
SiteJejunum                      -1.36936  15.86263 178 -0.086326  0.9313
SiteDuodenum                      7.87266  15.96680 178  0.493064  0.6226
MicrobiotaHumanized             -31.73684  15.58086  14 -2.036912  0.0610
TypeMucosal                      15.62500  12.47016 178  1.252991  0.2119
MicrobiotaHumanized:TypeMucosal -26.73624  16.36780 178 -1.633466  0.1041
```

Pielou's E
```R
Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.05624965 0.09011229

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site + Microbiota * Type 
                                     Value  Std.Error  DF   t-value p-value
(Intercept)                      0.5726199 0.03398333 178 16.850020  0.0000
Sequencing_Run2014_Sept         -0.0826440 0.01931094 178 -4.279648  0.0000
Sequencing_Run2015_Sept         -0.0241875 0.03935446  14 -0.614606  0.5487
SexMale                          0.0411537 0.03349570  14  1.228627  0.2395
SiteDistal_Colon                 0.0433397 0.02269747 178  1.909450  0.0578
SiteDuodenum                     0.0116058 0.02239052 178  0.518334  0.6049
SiteIleum                       -0.0660707 0.02232427 178 -2.959591  0.0035
SiteJejunum                     -0.0431836 0.02194308 178 -1.967983  0.0506
SiteProximal_Colon              -0.0077781 0.02141098 178 -0.363274  0.7168
MicrobiotaHumanized             -0.1001921 0.03178393  14 -3.152287  0.0071
TypeMucosal                      0.0934991 0.01935970 178  4.829571  0.0000
MicrobiotaHumanized:TypeMucosal  0.0276359 0.02539809 178  1.088109  0.2780
```

# Cedars SPF
### Microbiota * Site_General
shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)              4.614837 0.2816762 90 16.383485  0.0000
Sequencing_Run2014_Sept -0.427850 0.2119257 90 -2.018866  0.0465
Sequencing_Run2015_Sept  0.070869 0.3898668  6  0.181776  0.8617
SexMale                 -0.089417 0.3291647  6 -0.271649  0.7950
TypeMucosal              0.712720 0.1661470 90  4.289695  0.0000
Site_GeneralSI          -0.462537 0.1649144 90 -2.804708  0.0062
```

otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             183.10551 11.002154 90 16.642696  0.0000
Sequencing_Run2014_Sept  -2.39045 11.427512 90 -0.209184  0.8348
Sequencing_Run2015_Sept   2.07443 13.896044  6  0.149282  0.8862
SexMale                  13.18151 11.001999  6  1.198102  0.2761
TypeMucosal              25.53744  9.099814 90  2.806370  0.0061
Site_GeneralSI          -31.09404  9.023168 90 -3.446022  0.0009
```

chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             213.46834  14.33145 90 14.895100  0.0000
Sequencing_Run2014_Sept  13.48416  13.88005 90  0.971478  0.3339
Sequencing_Run2015_Sept  -9.11476  18.60971  6 -0.489785  0.6417
SexMale                  17.27971  15.05323  6  1.147907  0.2947
TypeMucosal               9.24692  10.99249 90  0.841203  0.4025
Site_GeneralSI          -27.59195  10.90366 90 -2.530522  0.0131
```

pielou_e
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6126518 0.03351622 90 18.279265  0.0000
Sequencing_Run2014_Sept -0.0542110 0.02442760 90 -2.219253  0.0290
Sequencing_Run2015_Sept  0.0110417 0.04664319  6  0.236726  0.8207
SexMale                 -0.0201927 0.03951058  6 -0.511071  0.6276
TypeMucosal              0.0806289 0.01913687 90  4.213276  0.0001
Site_GeneralSI          -0.0404525 0.01899584 90 -2.129546  0.0359
```

### For Type and Site Differences
shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)              4.786489 0.3334592 86 14.354050  0.0000
Sequencing_Run2014_Sept -0.485445 0.2212462 86 -2.194140  0.0309
Sequencing_Run2015_Sept  0.042660 0.3923589  6  0.108728  0.9170
SexMale                 -0.096417 0.3312277  6 -0.291090  0.7808
TypeMucosal              0.747025 0.1675288 86  4.459085  0.0000
SiteProximal_Colon      -0.289999 0.2567387 86 -1.129550  0.2618
SiteCecum               -0.151230 0.2567387 86 -0.589044  0.5574
SiteIleum               -0.904394 0.2823100 86 -3.203548  0.0019
SiteJejunum             -0.588247 0.2829162 86 -2.079228  0.0406
SiteDuodenum            -0.407781 0.2848747 86 -1.431440  0.1559
```

otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             183.73414 14.814561 86 12.402266  0.0000
Sequencing_Run2014_Sept  -2.88960 12.137317 86 -0.238075  0.8124
Sequencing_Run2015_Sept   2.14388 13.978620  6  0.153369  0.8831
SexMale                  13.04947 11.011940  6  1.185029  0.2808
TypeMucosal              25.74368  9.357797 86  2.751041  0.0072
SiteProximal_Colon       -4.26942 14.396454 86 -0.296560  0.7675
SiteCecum                 3.06392 14.396454 86  0.212824  0.8320
SiteIleum               -26.98808 15.787762 86 -1.709430  0.0910
SiteJejunum             -40.37392 15.796627 86 -2.555857  0.0123
SiteDuodenum            -27.12736 15.885801 86 -1.707648  0.0913
```

chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             214.36202  18.68983 86 11.469449  0.0000
Sequencing_Run2014_Sept  12.77435  14.74873 86  0.866132  0.3888
Sequencing_Run2015_Sept  -9.00016  18.70283  6 -0.481219  0.6474
SexMale                  17.13239  15.06622  6  1.137140  0.2988
TypeMucosal               9.50781  11.30146 86  0.841291  0.4025
SiteProximal_Colon       -4.32072  17.36405 86 -0.248831  0.8041
SiteCecum                 2.56421  17.36405 86  0.147674  0.8829
SiteIleum               -21.44736  19.05936 86 -1.125293  0.2636
SiteJejunum             -39.74045  19.07996 86 -2.082837  0.0402
SiteDuodenum            -23.40661  19.19575 86 -1.219364  0.2260
```

pielous
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.6349896 0.03910095 86 16.239746  0.0000
Sequencing_Run2014_Sept -0.0614603 0.02516755 86 -2.442047  0.0167
Sequencing_Run2015_Sept  0.0071826 0.04686890  6  0.153248  0.8832
SexMale                 -0.0210449 0.03973618  6 -0.529615  0.6154
TypeMucosal              0.0850635 0.01903755 86  4.468198  0.0000
SiteProximal_Colon      -0.0358649 0.02916860 86 -1.229573  0.2222
SiteCecum               -0.0218753 0.02916860 86 -0.749961  0.4553
SiteIleum               -0.1036941 0.03207884 86 -3.232477  0.0017
SiteJejunum             -0.0494506 0.03215078 86 -1.538083  0.1277
SiteDuodenum            -0.0354730 0.03237578 86 -1.095664  0.2763
```

# Humanized Dataset
### Site General
shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Site_General 
                             Value Std.Error DF   t-value p-value
(Intercept)              3.0629033 0.3506109 90  8.735905  0.0000
Sequencing_Run2014_Sept -0.6170179 0.2359160 90 -2.615413  0.0104
Sequencing_Run2015_Sept -0.6816864 0.4979781  6 -1.368908  0.2201
SexMale                  0.7092858 0.4221116  6  1.680327  0.1439
TypeMucosal              0.9075198 0.1925058 90  4.714248  0.0000
Site_GeneralSI           0.0700434 0.1824789 90  0.383844  0.7020
```
otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             128.15588  16.76163 90  7.645789  0.0000
Sequencing_Run2014_Sept -16.51825  13.18544 90 -1.252764  0.2135
Sequencing_Run2015_Sept -37.56606  23.27522  6 -1.613994  0.1577
SexMale                  -2.82639  19.38423  6 -0.145809  0.8888
TypeMucosal              10.86104  10.78609 90  1.006948  0.3167
Site_GeneralSI           57.70351  10.22904 90  5.641144  0.0000
```
chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Site_General 
                            Value Std.Error DF   t-value p-value
(Intercept)             165.59083  21.03733 90  7.871285  0.0000
Sequencing_Run2014_Sept  -8.13526  15.63645 90 -0.520276  0.6041
Sequencing_Run2015_Sept -54.89282  29.48314  6 -1.861838  0.1119
SexMale                  -6.82489  24.73417  6 -0.275929  0.7919
TypeMucosal             -10.51114  12.77765 90 -0.822619  0.4129
Site_GeneralSI           64.20588  12.11540 90  5.299525  0.0000
```

pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Site_General 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.4438339 0.04061389 90 10.928131  0.0000
Sequencing_Run2014_Sept -0.0788475 0.02920180 90 -2.700091  0.0083
Sequencing_Run2015_Sept -0.0533728 0.05719404  6 -0.933189  0.3867
SexMale                  0.1038166 0.04816246  6  2.155550  0.0745
TypeMucosal              0.1186957 0.02385020 90  4.976717  0.0000
Site_GeneralSI          -0.0368116 0.02261183 90 -1.627979  0.1070
```
### For Type and Site Differences
shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)              3.551963 0.3898301 86  9.111566  0.0000
Sequencing_Run2014_Sept -0.845622 0.2399069 86 -3.524791  0.0007
Sequencing_Run2015_Sept -0.680968 0.4776085  6 -1.425787  0.2038
SexMale                  0.709286 0.4043808  6  1.754005  0.1300
TypeMucosal              1.022198 0.1887016 86  5.417008  0.0000
SiteProximal_Colon      -0.582669 0.2873716 86 -2.027579  0.0457
SiteCecum               -0.656218 0.2858449 86 -2.295713  0.0241
SiteIleum               -0.567427 0.3144356 86 -1.804589  0.0746
SiteJejunum             -0.686568 0.3059775 86 -2.243851  0.0274
SiteDuodenum            -0.003871 0.3076412 86 -0.012584  0.9900
```

otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             148.70583  19.98613 86  7.440451  0.0000
Sequencing_Run2014_Sept -24.17238  13.90939 86 -1.737846  0.0858
Sequencing_Run2015_Sept -38.02680  22.64797  6 -1.679038  0.1441
SexMale                  -2.82639  18.77762  6 -0.150519  0.8853
TypeMucosal              15.15610  10.97499 86  1.380967  0.1709
SiteProximal_Colon      -21.49219  16.73659 86 -1.284144  0.2025
SiteCecum               -32.32663  16.64674 86 -1.941920  0.0554
SiteIleum                47.10417  18.29958 86  2.574057  0.0118
SiteJejunum              28.74617  17.80900 86  1.614138  0.1102
SiteDuodenum             34.02947  17.90481 86  1.900577  0.0607
```

chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Type + Site 
                            Value Std.Error DF   t-value p-value
(Intercept)             181.48342  24.90554 86  7.286871  0.0000
Sequencing_Run2014_Sept -13.39316  16.61901 86 -0.805894  0.4225
Sequencing_Run2015_Sept -55.02675  29.09966  6 -1.890976  0.1075
SexMale                  -6.82489  24.33601  6 -0.280444  0.7886
TypeMucosal              -7.30287  13.09566 86 -0.557655  0.5785
SiteProximal_Colon      -16.44746  19.95916 86 -0.824056  0.4122
SiteCecum               -26.08310  19.85247 86 -1.313847  0.1924
SiteIleum                63.39385  21.82970 86  2.904018  0.0047
SiteJejunum              38.71151  21.24364 86  1.822263  0.0719
SiteDuodenum             41.20726  21.35842 86  1.929322  0.0570
```

pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Type + Site 
                             Value  Std.Error DF   t-value p-value
(Intercept)              0.4974812 0.04598452 86 10.818448  0.0000
Sequencing_Run2014_Sept -0.1053090 0.02942092 86 -3.579391  0.0006
Sequencing_Run2015_Sept -0.0527701 0.05515855  6 -0.956698  0.3757
SexMale                  0.1038166 0.04645222  6  2.234911  0.0668
TypeMucosal              0.1316127 0.02315941 86  5.682905  0.0000
SiteProximal_Colon      -0.0678709 0.03528133 86 -1.923705  0.0577
SiteCecum               -0.0669219 0.03509339 86 -1.906965  0.0599
SiteIleum               -0.1187959 0.03859708 86 -3.077846  0.0028
SiteJejunum             -0.1238420 0.03755971 86 -3.297202  0.0014
SiteDuodenum            -0.0274854 0.03776336 86 -0.727832  0.4687
```