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