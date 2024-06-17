# Luminal Data
## Site General 
 observed_otus
```R
 > summary(output)
Linear mixed-effects model fit by REML
 Data: luminaldata 
      AIC      BIC    logLik
  608.917 621.0691 -298.4585

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    9.070334 44.28682

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                       Value Std.Error DF   t-value p-value
(Intercept)        247.65303  10.34800 49 23.932455  0.0000
Sequencing_RunTwo   17.34848  14.93911  7  1.161279  0.2836
SexM                 6.43939  13.97426  7  0.460804  0.6589
Site_GeneralSI    -110.36667  11.43481 49 -9.651816  0.0000
```
pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: luminaldata 
        AIC       BIC   logLik
  -75.63462 -63.48251 43.81731

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.07030576 0.09043891

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.6408944 0.03581736 49 17.893399  0.0000
Sequencing_RunTwo  0.0574716 0.05864888  7  0.979928  0.3598
SexM               0.0043635 0.05486100  7  0.079537  0.9388
Site_GeneralSI    -0.0744898 0.02335123 49 -3.189972  0.0025
 Correlation: 
                  (Intr) Sqn_RT SexM  
Sequencing_RunTwo -0.273              
SexM              -0.438 -0.356       
Site_GeneralSI    -0.326  0.000  0.000

Standardized Within-Group Residuals:
        Min          Q1         Med          Q3         Max 
-3.07325189 -0.52052184  0.04660606  0.49565412  1.92253384 

Number of Observations: 60
Number of Groups: 10 
```

## Site
observed_otus
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: luminaldata 
       AIC      BIC    logLik
  584.3963 603.9087 -282.1981

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    8.144069 45.35408

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                        Value Std.Error DF   t-value p-value
(Intercept)         240.41970  15.67828 45 15.334571  0.0000
Sequencing_RunTwo    17.34848  14.93911  7  1.161279  0.2836
SexM                  6.43939  13.97426  7  0.460804  0.6589
SiteProximal_Colon   12.40000  20.28296 45  0.611351  0.5440
SiteCecum             9.30000  20.28296 45  0.458513  0.6488
SiteIleum          -109.70000  20.28296 45 -5.408480  0.0000
SiteJejunum         -89.70000  20.28296 45 -4.422431  0.0001
SiteDuodenum       -110.00000  20.28296 45 -5.423271  0.0000
```

pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: luminaldata 
        AIC       BIC   logLik
  -49.92787 -30.41543 34.96394

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.06967599 0.09331802

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                        Value  Std.Error DF   t-value p-value
(Intercept)         0.6341451 0.04326956 45 14.655687  0.0000
Sequencing_RunTwo   0.0574716 0.05864897  7  0.979926  0.3598
SexM                0.0043635 0.05486108  7  0.079537  0.9388
SiteProximal_Colon -0.0066621 0.04173309 45 -0.159636  0.8739
SiteCecum           0.0269099 0.04173309 45  0.644809  0.5223
SiteIleum          -0.0679656 0.04173309 45 -1.628578  0.1104
SiteJejunum        -0.0790070 0.04173309 45 -1.893150  0.0648
SiteDuodenum       -0.0562490 0.04173309 45 -1.347827  0.1845
```
# Mucosal Data
## Site General
observed_otus
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: mucosaldata 
       AIC     BIC    logLik
  480.9229 491.628 -234.4614

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:     6.08513 43.77452

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                       Value Std.Error DF   t-value p-value
(Intercept)        270.28220  10.99638 39 24.579202  0.0000
Sequencing_RunTwo   24.48485  15.58839  5  1.570710  0.1770
SexM                 8.65152  15.58839  5  0.554997  0.6028
Site_GeneralSI    -112.29167  12.63662 39 -8.886213  0.0000
```

pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: mucosaldata 
        AIC       BIC   logLik
  -70.92476 -60.21962 41.46238

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.04012527 0.07895392

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.6158774 0.02705388 39 22.764844  0.0000
Sequencing_RunTwo  0.0313157 0.04249859  5  0.736865  0.4943
SexM              -0.0174612 0.04249859  5 -0.410867  0.6982
Site_GeneralSI    -0.0280415 0.02279203 39 -1.230322  0.2259
```

## Site
observed_otus
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: mucosaldata 
       AIC      BIC    logLik
  447.6553 464.5441 -213.8276

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    9.264601 40.29125

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                        Value Std.Error DF   t-value p-value
(Intercept)         296.19886  15.81461 35 18.729443  0.0000
Sequencing_RunTwo    24.48485  15.58839  5  1.570710  0.1770
SexM                  8.65152  15.58839  5  0.554997  0.6028
SiteProximal_Colon  -31.62500  20.14563 35 -1.569820  0.1255
SiteCecum           -46.12500  20.14563 35 -2.289579  0.0282
SiteIleum          -144.37500  20.14563 35 -7.166568  0.0000
SiteJejunum        -158.25000  20.14563 35 -7.855303  0.0000
SiteDuodenum       -112.00000  20.14563 35 -5.559519  0.0000
```

pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: mucosaldata 
        AIC       BIC logLik
  -52.10401 -35.21522 36.052

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.04158231 0.07429291

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                        Value  Std.Error DF   t-value p-value
(Intercept)         0.6288463 0.03430722 35 18.329851  0.0000
Sequencing_RunTwo   0.0313157 0.04249862  5  0.736864  0.4943
SexM               -0.0174612 0.04249862  5 -0.410866  0.6982
SiteProximal_Colon  0.0207487 0.03714645 35  0.558565  0.5800
SiteCecum          -0.0596555 0.03714645 35 -1.605953  0.1173
SiteIleum          -0.0807179 0.03714645 35 -2.172963  0.0366
SiteJejunum        -0.0351590 0.03714645 35 -0.946498  0.3504
SiteDuodenum       -0.0071545 0.03714645 35 -0.192602  0.8484
```


# Duodenum
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex +  Type, random = ~1|(MouseID), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
       AIC      BIC    logLik
  159.5745 163.4088 -73.78723

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    4.711643 36.34728

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF  t-value p-value
(Intercept)       115.99410  13.92187  7 8.331791  0.0001
Sequencing_RunTwo  58.70899  20.26616  7 2.896897  0.0231
SexM               11.48302  19.53146  7 0.587924  0.5751
TypeMucosal        54.28314  17.34492  7 3.129628  0.0166
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=duodata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: duodata 
        AIC        BIC   logLik
  -4.243047 -0.4087032 8.121524

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.07141878 0.08689709

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.5943123 0.04493139  7 13.227107  0.0000
Sequencing_RunTwo  0.0542858 0.07160072  7  0.758174  0.4731
SexM              -0.0342875 0.06815915  7 -0.503050  0.6304
TypeMucosal        0.0265723 0.04229566  7  0.628251  0.5498
```

# Jejunum
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
       AIC      BIC    logLik
  165.5383 169.3726 -76.76915

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    9.224437 44.49015

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF   t-value p-value
(Intercept)       148.64963  17.31844  7  8.583317  0.0001
Sequencing_RunTwo   4.33514  25.39282  7  0.170723  0.8693
SexM               21.37456  24.45238  7  0.874130  0.4110
TypeMucosal        -8.13563  21.25891  7 -0.382693  0.7133
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex +  Type, random = ~1|(MouseID), data=jejdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: jejdata 
        AIC        BIC   logLik
  -4.820588 -0.9862437 8.410294

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:   0.1300459 0.05932362

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                      Value  Std.Error DF  t-value p-value
(Intercept)       0.5475983 0.05988951  7 9.143476  0.0000
Sequencing_RunTwo 0.0680806 0.10172336  7 0.669272  0.5248
SexM              0.0152563 0.09561845  7 0.159554  0.8777
TypeMucosal       0.0285913 0.02944416  7 0.971035  0.3639
```

# Ileum
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
       AIC      BIC    logLik
  169.0051 172.8395 -78.50257

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:     36.0917 41.68776

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF   t-value p-value
(Intercept)       133.63991  22.09361  7  6.048804  0.0005
Sequencing_RunTwo  49.87545  35.36995  7  1.410108  0.2014
SexM              -25.25637  33.64368  7 -0.750702  0.4773
TypeMucosal        21.64444  20.31560  7  1.065410  0.3221
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex +  Type, random = ~1|(MouseID), data=iledata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: iledata 
      AIC      BIC   logLik
  3.31871 7.153054 4.340645

Random effects:
 Formula: ~1 | MouseID
        (Intercept)  Residual
StdDev:   0.1099109 0.1064563

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.5301733 0.06194597  7  8.558640  0.0001
Sequencing_RunTwo  0.1018489 0.10067271  7  1.011683  0.3454
SexM               0.0610962 0.09550418  7  0.639723  0.5427
TypeMucosal       -0.0357908 0.05209601  7 -0.687017  0.5142
```

# Cecum
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
       AIC      BIC    logLik
  168.5388 172.3731 -78.26938

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    15.37453   48.429

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF   t-value p-value
(Intercept)       249.31857  19.50157  7 12.784535  0.0000
Sequencing_RunTwo  -5.99841  28.99565  7 -0.206873  0.8420
SexM               24.95239  27.87583  7  0.895126  0.4005
TypeMucosal         6.08983  23.20310  7  0.262458  0.8005
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
        AIC       BIC   logLik
  -19.21235 -15.37801 15.60618

Random effects:
 Formula: ~1 | MouseID
         (Intercept)   Residual
StdDev: 7.493305e-06 0.06173567

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.6758180 0.02339548  7 28.886690  0.0000
Sequencing_RunTwo  0.0322398 0.03388842  7  0.951351  0.3731
SexM              -0.0136203 0.03267772  7 -0.416806  0.6893
TypeMucosal       -0.1084141 0.02943392  7 -3.683305  0.0078> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=cecdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: cecdata 
        AIC       BIC   logLik
  -19.21235 -15.37801 15.60618

Random effects:
 Formula: ~1 | MouseID
         (Intercept)   Residual
StdDev: 7.493305e-06 0.06173567

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF   t-value p-value
(Intercept)        0.6758180 0.02339548  7 28.886690  0.0000
Sequencing_RunTwo  0.0322398 0.03388842  7  0.951351  0.3731
SexM              -0.0136203 0.03267772  7 -0.416806  0.6893
TypeMucosal       -0.1084141 0.02943392  7 -3.683305  0.0078
```

# Proximal Colon
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
       AIC      BIC    logLik
  147.1369 150.9712 -67.56845

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev: 0.002693956 23.48005

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF   t-value p-value
(Intercept)       258.97093  8.898052  7 29.104226  0.0000
Sequencing_RunTwo -13.83789 12.888851  7 -1.073632  0.3186
SexM               14.45110 12.428385  7  1.162750  0.2830
TypeMucosal        17.79912 11.194662  7  1.589965  0.1559
```


```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=pcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: pcdata 
        AIC       BIC logLik
  -31.56199 -27.72765 21.781

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.01608641 0.03699802

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF  t-value p-value
(Intercept)        0.6476630 0.01562521  7 41.44987  0.0000
Sequencing_RunTwo -0.0080091 0.02363761  7 -0.33883  0.7447
SexM               0.0030240 0.02267517  7  0.13336  0.8977
TypeMucosal        0.0103200 0.01778900  7  0.58013  0.5800
```

# Distal Colon
```R
> output=lme(fixed= observed_otus ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
       AIC      BIC    logLik
  163.9903 167.8247 -75.99517

Random effects:
 Formula: ~1 | MouseID
        (Intercept) Residual
StdDev:    12.31089 41.35465

Fixed effects:  observed_otus ~ Sequencing_Run + Sex + Type 
                      Value Std.Error DF   t-value p-value
(Intercept)       239.27405  16.53765  7 14.468442  0.0000
Sequencing_RunTwo  30.36794  24.52153  7  1.238420  0.2555
SexM               -0.46109  23.58230  7 -0.019552  0.9849
TypeMucosal        58.66413  19.80303  7  2.962382  0.0210
```

```R
> output=lme(fixed= pielou_e ~ Sequencing_Run + Sex + Type, random = ~1|(MouseID), data=dcdata)
> summary(output)
Linear mixed-effects model fit by REML
  Data: dcdata 
        AIC       BIC   logLik
  -28.89097 -25.05663 20.44549

Random effects:
 Formula: ~1 | MouseID
        (Intercept)   Residual
StdDev:  0.03590228 0.03318748

Fixed effects:  pielou_e ~ Sequencing_Run + Sex + Type 
                       Value  Std.Error DF  t-value p-value
(Intercept)        0.6647442 0.01984651  7 33.49425  0.0000
Sequencing_RunTwo -0.0109278 0.03237591  7 -0.33753  0.7456
SexM              -0.0208346 0.03069217  7 -0.67882  0.5191
TypeMucosal       -0.0177424 0.01625856  7 -1.09127  0.3113
```