# Luminal Data
## Site General 
 observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Line + Site_General 
                                      Value Std.Error  DF    t-value p-value
(Intercept)                       284.01978  14.22243 230  19.969846  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   14.06281  25.82820  40   0.544475  0.5891
Sequencing_RunNovaSeq_Mar_Twenty   68.39122  45.34246  40   1.508326  0.1393
SexMale                            -1.18161  14.23149  40  -0.083028  0.9342
LineJJWT                          -57.58432  32.42733  40  -1.775796  0.0834
LineLyzCre                        -15.18300  20.17116  40  -0.752708  0.4560
LineVil1Cre                        21.89573  18.92859  40   1.156754  0.2542
Site_GeneralSI                   -115.30860   8.43624 230 -13.668242  0.0000
```
 pielou e 
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Line + Site_General 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.6224884 0.02587909 230 24.053719  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1052799 0.04686105  40 -2.246640  0.0303
Sequencing_RunNovaSeq_Mar_Twenty -0.0084258 0.08231583  40 -0.102359  0.9190
SexMale                          -0.0072378 0.02582171  40 -0.280300  0.7807
LineJJWT                          0.0773085 0.05883958  40  1.313887  0.1964
LineLyzCre                       -0.0330798 0.03660355  40 -0.903731  0.3715
LineVil1Cre                       0.0002547 0.03434228  40  0.007415  0.9941
Site_GeneralSI                   -0.1272117 0.01579512 230 -8.053860  0.0000
```

## Site 

Observed OTUs
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Line + Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                       267.15647  16.39782 226 16.292195  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   13.70202  25.78820  40  0.531329  0.5981
Sequencing_RunNovaSeq_Mar_Twenty   70.33344  45.26814  40  1.553707  0.1281
SexMale                            -1.51910  14.21006  40 -0.106903  0.9154
LineJJWT                          -57.59840  32.37590  40 -1.779052  0.0828
LineLyzCre                        -15.00318  20.13821  40 -0.745011  0.4606
LineVil1Cre                        21.69788  18.89930  40  1.148079  0.2578
SiteProximal_Colon                 22.17021  14.25769 226  1.554966  0.1214
SiteCecum                          29.04255  14.25769 226  2.036975  0.0428
SiteIleum                        -101.83607  14.34420 226 -7.099458  0.0000
SiteJejunum                      -115.10638  14.25769 226 -8.073286  0.0000
SiteDuodenum                      -76.31049  14.52881 226 -5.252357  0.0000
```

Evenness
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Line + Site 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.5901917 0.02995343 226 19.703640  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1060389 0.04681353  40 -2.265134  0.0290
Sequencing_RunNovaSeq_Mar_Twenty -0.0031023 0.08220543  40 -0.037739  0.9701
SexMale                          -0.0082339 0.02579629  40 -0.319189  0.7512
LineJJWT                          0.0767318 0.05877552  40  1.305507  0.1992
LineLyzCre                       -0.0326631 0.03656077  40 -0.893391  0.3770
LineVil1Cre                      -0.0001375 0.03430771  40 -0.004007  0.9968
SiteProximal_Colon                0.0310319 0.02638968 226  1.175911  0.2409
SiteCecum                         0.0677844 0.02638968 226  2.568595  0.0109
SiteIleum                        -0.1234502 0.02654939 226 -4.649831  0.0000
SiteJejunum                      -0.1141412 0.02638968 226 -4.325221  0.0000
SiteDuodenum                     -0.0422443 0.02689029 226 -1.570988  0.1176
```

# Mucosal Data
## Site General
observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Line + Site_General 
                                     Value Std.Error  DF    t-value p-value
(Intercept)                       333.3015  20.33221 221  16.392785  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -21.7557  36.15409  41  -0.601749  0.5507
Sequencing_RunNovaSeq_Mar_Twenty   73.3117  65.09449  41   1.126235  0.2666
SexMale                           -26.2149  20.55399  41  -1.275416  0.2093
LineJJWT                          -48.2759  46.38597  41  -1.040744  0.3041
LineLyzCre                         -1.9110  29.59440  41  -0.064573  0.9488
LineVil1Cre                      -104.6789  27.43021  41  -3.816189  0.0004
Site_GeneralSI                   -112.9115   7.94402 221 -14.213402  0.0000
```
pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Line + Site_General 
                                      Value  Std.Error  DF    t-value p-value
(Intercept)                       0.6572648 0.04607934 221  14.263764  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1948669 0.08286512  41  -2.351616  0.0236
Sequencing_RunNovaSeq_Mar_Twenty -0.0105279 0.14900960  41  -0.070653  0.9440
SexMale                          -0.0282055 0.04694672  41  -0.600799  0.5513
LineJJWT                          0.1058377 0.10600875  41   0.998386  0.3239
LineLyzCre                        0.0182160 0.06775546  41   0.268850  0.7894
LineVil1Cre                      -0.2131558 0.06268307  41  -3.400532  0.0015
Site_GeneralSI                   -0.1614842 0.01295037 221 -12.469470  0.0000
```

## Site
observed_otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Line + Site 
                                     Value Std.Error  DF   t-value p-value
(Intercept)                       336.5698  21.99857 217 15.299622  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -21.3235  36.25954  41 -0.588079  0.5597
Sequencing_RunNovaSeq_Mar_Twenty   74.3593  65.29218  41  1.138870  0.2614
SexMale                           -27.1579  20.61902  41 -1.317127  0.1951
LineJJWT                          -49.5515  46.53391  41 -1.064847  0.2932
LineLyzCre                         -2.7284  29.68290  41 -0.091917  0.9272
LineVil1Cre                      -104.9537  27.51181  41 -3.814859  0.0005
SiteProximal_Colon                  5.8692  13.76585 217  0.426358  0.6703
SiteCecum                         -12.6100  13.76585 217 -0.916034  0.3607
SiteIleum                        -122.2970  13.90213 217 -8.796997  0.0000
SiteJejunum                      -122.1208  14.08299 217 -8.671511  0.0000
SiteDuodenum                     -100.1598  14.28480 217 -7.011639  0.0000
```
pielou_e
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Line + Site 
                                      Value  Std.Error  DF   t-value p-value
(Intercept)                       0.6068318 0.04759583 217 12.749684  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1947809 0.08256543  41 -2.359109  0.0232
Sequencing_RunNovaSeq_Mar_Twenty -0.0049098 0.14846086  41 -0.033072  0.9738
SexMale                          -0.0283932 0.04676635  41 -0.607129  0.5471
LineJJWT                          0.0993174 0.10560888  41  0.940427  0.3525
LineLyzCre                        0.0171966 0.06750548  41  0.254744  0.8002
LineVil1Cre                      -0.2157678 0.06243970  41 -3.455618  0.0013
SiteProximal_Colon                0.1135870 0.02121844 217  5.353222  0.0000
SiteCecum                         0.0364887 0.02121844 217  1.719669  0.0869
SiteIleum                        -0.1232674 0.02143038 217 -5.751994  0.0000
SiteJejunum                      -0.1068347 0.02171726 217 -4.919346  0.0000
SiteDuodenum                     -0.0952896 0.02202821 217 -4.325799  0.0000
```