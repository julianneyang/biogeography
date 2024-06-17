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

# Duodenum 
observed otus
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
       AIC      BIC    logLik
  982.0756 1005.643 -481.0378

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev: 0.007930672 100.8855

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Line + Type 
                                     Value Std.Error DF   t-value p-value
(Intercept)                      185.40294  26.73395 41  6.935112  0.0000
Sequencing_RunNovaSeq_Jan_Twenty   2.38604  41.26109 41  0.057828  0.9542
Sequencing_RunNovaSeq_Mar_Twenty 173.05318  78.30894 41  2.209878  0.0328
SexMale                            4.46594  23.76864 37  0.187892  0.8520
LineJJWT                         -79.97345  53.58877 41 -1.492355  0.1433
LineLyzCre                        10.96958  34.93692 41  0.313982  0.7551
LineVil1Cre                       24.88972  32.19140 41  0.773179  0.4439
TypeMucosal                      -16.42292  21.90719 37 -0.749659  0.4582
```
pielou_e
```R
 output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
        AIC       BIC   logLik
  -25.41986 -1.852769 22.70993

Random effects:
 Formula: ~1 | MouseID_Line
         (Intercept)  Residual
StdDev: 4.348744e-07 0.1581332

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Line + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.5476955 0.04190418 41 13.070188  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1485828 0.06467477 41 -2.297384  0.0268
Sequencing_RunNovaSeq_Mar_Twenty -0.0102072 0.12274549 41 -0.083158  0.9341
SexMale                           0.0316624 0.03725619 37  0.849857  0.4009
LineJJWT                          0.1048005 0.08399781 41  1.247657  0.2192
LineLyzCre                       -0.0360430 0.05476193 41 -0.658177  0.5141
LineVil1Cre                      -0.0527192 0.05045846 41 -1.044804  0.3022
TypeMucosal                      -0.1214349 0.03433846 37 -3.536409  0.0011
```

# Jejunum

observed_otus
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
      AIC      BIC    logLik
  1016.73 1040.918 -498.3649

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev: 0.006699499 86.15639

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Line + Type 
                                     Value Std.Error DF   t-value p-value
(Intercept)                      157.95410  21.00438 42  7.520056  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -20.84285  34.92066 42 -0.596863  0.5538
Sequencing_RunNovaSeq_Mar_Twenty  51.65550  62.99908 42  0.819941  0.4169
SexMale                            1.87446  20.08818 41  0.093312  0.9261
LineJJWT                         -37.61745  45.53443 42 -0.826132  0.4134
LineLyzCre                        -3.79041  28.08732 42 -0.134951  0.8933
LineVil1Cre                        2.95977  26.42676 42  0.111999  0.9114
TypeMucosal                       10.26677  18.17068 41  0.565019  0.5751
```

pielou_e
```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
        AIC       BIC   logLik
  -26.77588 -2.587471 23.38794

Random effects:
 Formula: ~1 | MouseID_Line
         (Intercept)  Residual
StdDev: 1.165165e-05 0.1604112

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Line + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.4661713 0.03910722 42 11.920339  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1357049 0.06501740 42 -2.087209  0.0430
Sequencing_RunNovaSeq_Mar_Twenty  0.0018625 0.11729550 42  0.015879  0.9874
SexMale                           0.0284935 0.03740138 41  0.761830  0.4505
LineJJWT                          0.1045368 0.08477876 42  1.233054  0.2244
LineLyzCre                       -0.0146913 0.05229467 42 -0.280934  0.7801
LineVil1Cre                      -0.0523499 0.04920295 42 -1.063958  0.2934
TypeMucosal                      -0.0483948 0.03383128 41 -1.430476  0.1602
```

# Ileum
observed_otus
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
       AIC      BIC    logLik
  1023.269 1047.577 -501.6345

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev: 0.008210716 83.46051

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Line + Type 
                                     Value Std.Error DF   t-value p-value
(Intercept)                      176.22448  20.17250 42  8.735875  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -21.64919  33.77768 42 -0.640932  0.5250
Sequencing_RunNovaSeq_Mar_Twenty  80.36954  60.27338 42  1.333417  0.1896
SexMale                          -15.32318  19.10242 42 -0.802159  0.4270
LineJJWT                         -11.61242  43.41577 42 -0.267470  0.7904
LineLyzCre                        -3.31987  26.94852 42 -0.123193  0.9025
LineVil1Cre                       -4.54378  25.32036 42 -0.179452  0.8584
TypeMucosal                       -2.31684  17.41552 42 -0.133033  0.8948
```
pielou_e
```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line+ Type, random = ~1|(MouseID_Line), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
        AIC      BIC  logLik
  -21.63841 2.669762 20.8192

Random effects:
 Formula: ~1 | MouseID_Line
         (Intercept)  Residual
StdDev: 5.798785e-06 0.1660754

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Line + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.4749699 0.04014063 42 11.832647  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1404072 0.06721313 42 -2.088984  0.0428
Sequencing_RunNovaSeq_Mar_Twenty  0.0941410 0.11993610 42  0.784926  0.4369
SexMale                          -0.0115054 0.03801130 42 -0.302684  0.7636
LineJJWT                          0.0909946 0.08639168 42  1.053280  0.2982
LineLyzCre                       -0.0281421 0.05362401 42 -0.524805  0.6025
LineVil1Cre                      -0.0221278 0.05038419 42 -0.439182  0.6628
TypeMucosal                      -0.0552541 0.03465459 42 -1.594423  0.1183
```

# Cecum
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line+ Type, random = ~1|(MouseID_Line), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
       AIC      BIC    logLik
  1019.326 1043.985 -499.6629

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev:      14.696 65.11949

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Line + Type 
                                    Value Std.Error DF   t-value p-value
(Intercept)                      331.3476  16.34348 45 20.273985  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -22.6990  27.60826 42 -0.822182  0.4156
Sequencing_RunNovaSeq_Mar_Twenty   4.1158  48.85670 42  0.084243  0.9333
SexMale                          -11.7547  15.39190 45 -0.763693  0.4490
LineJJWT                         -17.9460  34.93406 42 -0.513712  0.6101
LineLyzCre                       -17.3019  22.07101 42 -0.783919  0.4375
LineVil1Cre                      -73.8003  20.55943 42 -3.589610  0.0009
TypeMucosal                      -22.5253  13.37470 45 -1.684174  0.0991
```
pielou
```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
        AIC       BIC   logLik
  -40.47711 -15.81803 30.23856

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept)  Residual
StdDev: 0.009674158 0.1504846

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Line + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.7123840 0.03640998 45 19.565625  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.1980929 0.06109651 42 -3.242295  0.0023
Sequencing_RunNovaSeq_Mar_Twenty -0.1015920 0.10804377 42 -0.940285  0.3524
SexMale                          -0.0538525 0.03406594 45 -1.580832  0.1209
LineJJWT                          0.1594535 0.07727658 42  2.063414  0.0453
LineLyzCre                        0.0276833 0.04878967 42  0.567400  0.5735
LineVil1Cre                      -0.1325836 0.04545248 42 -2.916973  0.0057
TypeMucosal                      -0.0866616 0.03089984 45 -2.804597  0.0074
```

# Proximal Colon
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Line + Type, random = ~1|(MouseID_Line), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
       AIC      BIC    logLik
  1019.239 1043.898 -499.6193

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev:    30.42038 60.51185

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Line + Type 
                                    Value Std.Error DF   t-value p-value
(Intercept)                      327.9306  17.34299 45 18.908534  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  39.3838  29.87421 42  1.318320  0.1945
Sequencing_RunNovaSeq_Mar_Twenty  63.8357  52.99118 42  1.204648  0.2351
SexMale                          -23.6539  16.63289 45 -1.422116  0.1619
LineJJWT                         -83.6081  37.85643 42 -2.208558  0.0327
LineLyzCre                       -20.6186  23.97284 42 -0.860081  0.3946
LineVil1Cre                      -80.2965  22.32214 42 -3.597166  0.0008
TypeMucosal                        1.4914  12.43924 45  0.119895  0.9051
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
        AIC       BIC   logLik
  -36.82194 -19.32328 25.41097

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept)  Residual
StdDev:  0.07143083 0.1538029

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.6420493 0.03379787 45 18.996741  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0106798 0.04252482 45 -0.251142  0.8028
Sequencing_RunNovaSeq_Mar_Twenty  0.1202180 0.09692236 45  1.240353  0.2213
SexMale                          -0.0896369 0.03878155 45 -2.311328  0.0255
TypeMucosal                       0.0243982 0.03160028 45  0.772088  0.4441
```

# Distal Colon
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
       AIC      BIC    logLik
  1054.402 1071.417 -520.2009

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept) Residual
StdDev:    11.79856 107.7962

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                                     Value Std.Error DF   t-value p-value
(Intercept)                      286.15565  21.20509 45 13.494669  0.0000
Sequencing_RunNovaSeq_Jan_Twenty  -5.48743  26.93195 45 -0.203752  0.8395
Sequencing_RunNovaSeq_Mar_Twenty  52.06941  57.61489 45  0.903749  0.3709
SexMale                          -44.45020  23.95572 39 -1.855515  0.0711
TypeMucosal                       20.45028  22.92213 39  0.892163  0.3778
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID_Line), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
        AIC      BIC   logLik
  -8.915565 8.100153 11.45778

Random effects:
 Formula: ~1 | MouseID_Line
        (Intercept)  Residual
StdDev:  0.06870605 0.1818868

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                                      Value  Std.Error DF   t-value p-value
(Intercept)                       0.6228402 0.03886834 45 16.024357  0.0000
Sequencing_RunNovaSeq_Jan_Twenty -0.0640143 0.05043291 45 -1.269296  0.2109
Sequencing_RunNovaSeq_Mar_Twenty  0.1497504 0.10886945 45  1.375504  0.1758
SexMale                          -0.0855745 0.04493802 39 -1.904279  0.0643
TypeMucosal                      -0.0624075 0.03884792 39 -1.606456  0.1162
```