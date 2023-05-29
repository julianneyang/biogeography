## Site _ General
Observed OTUs
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: data 
      AIC      BIC    logLik
  1442.71 1460.141 -715.3549

Random effects:
 Formula: ~1 | MouseID_Original
        (Intercept) Residual
StdDev:    12.42228 44.40405

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site_General 
                                  Value Std.Error  DF    t-value p-value
(Intercept)                   228.09481  8.562297 112  26.639442  0.0000
Sequencing_Run2014-September   30.08141  9.564540  23   3.145097  0.0045
SexMale                        -3.59949  9.285628  23  -0.387641  0.7018
Site_GeneralSI               -122.62367  7.616030 112 -16.100733  0.0000
```

Pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: data 
        AIC    BIC   logLik
  -86.22165 -68.79 49.11082

Random effects:
 Formula: ~1 | MouseID_Original
        (Intercept)  Residual
StdDev:  0.03985962 0.1548184

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site_General 
                                  Value  Std.Error  DF   t-value p-value
(Intercept)                   0.6472712 0.02927459 112 22.110345  0.0000
Sequencing_Run2014-September  0.0091154 0.03256649  23  0.279902  0.7821
SexMale                      -0.0380550 0.03161690  23 -1.203628  0.2410
Site_GeneralSI               -0.1979685 0.02653280 112 -7.461276  0.0000
```
## Site
Observed OTUs
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: data 
       AIC      BIC    logLik
  1400.134 1428.886 -690.0668

Random effects:
 Formula: ~1 | MouseID_Original
        (Intercept) Residual
StdDev:    14.15083  40.7384

Fixed effects: observed_otus ~ Sequencing_Run + Sex + Site 
                                  Value Std.Error  DF    t-value p-value
(Intercept)                   221.02930 11.189679 108  19.752961  0.0000
Sequencing_Run2014-September   28.96770  9.469631  23   3.059010  0.0056
SexMale                        -2.68464  9.195968  23  -0.291937  0.7730
SiteProximal_Colon             16.81271 12.444476 108   1.351018  0.1795
SiteCecum                       4.75286 12.106749 108   0.392580  0.6954
SiteIleum                    -136.31636 12.112238 108 -11.254432  0.0000
SiteJejunum                  -124.69295 12.248340 108 -10.180396  0.0000
SiteDuodenum                  -80.61989 12.531475 108  -6.433392  0.0000
```
pielou_e
```R
> summary(output)
Linear mixed-effects model fit by REML
 Data: data 
        AIC       BIC   logLik
  -117.4772 -88.72525 68.73861

Random effects:
 Formula: ~1 | MouseID_Original
        (Intercept)  Residual
StdDev:  0.05663522 0.1213741

Fixed effects: pielou_e ~ Sequencing_Run + Sex + Site 
                                  Value  Std.Error  DF   t-value p-value
(Intercept)                   0.6285978 0.03583963 108 17.539183  0.0000
Sequencing_Run2014-September  0.0006761 0.03240389  23  0.020864  0.9835
SexMale                      -0.0376454 0.03147443  23 -1.196063  0.2439
SiteProximal_Colon            0.0384658 0.03709988 108  1.036818  0.3021
SiteCecum                     0.0219375 0.03612777 108  0.607219  0.5450
SiteIleum                    -0.3124990 0.03613959 108 -8.647000  0.0000
SiteJejunum                  -0.1806422 0.03655852 108 -4.941178  0.0000
SiteDuodenum                 -0.0070190 0.03743476 108 -0.187498  0.8516
```