## Mucosal 
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   5     88.21  17.642  31.918 0.15363    0.0451 *  
Donor_ID        10    329.44  32.944  59.602 0.57376 9.999e-05 ***
Site_General     1     17.24  17.236  31.184 0.03002 9.999e-05 ***
Residuals      252    139.29   0.553         0.24259 9.999e-05 ***
Total          268    574.18                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
## Luminal
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   5     45.09   9.018  16.652 0.08232    0.2935    
Donor_ID         9    259.96  28.884  53.335 0.47463 9.999e-05 ***
Site_General     1    109.43 109.434 202.072 0.19980 9.999e-05 ***
Residuals      246    133.22   0.542         0.24324 9.999e-05 ***
Total          261    547.71                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
## Mucosal SI
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   5    97.484 19.4969 23.3977 0.32693 9.999e-05 ***
Donor_ID        10    99.654  9.9654 11.9592 0.33420    0.0012 ** 
Site             2     1.053  0.5265  0.6318 0.00353    0.6649    
Residuals      120    99.994  0.8333         0.33534 9.999e-05 ***
Total          137   298.185                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Mucosal Colon
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   4    39.707  9.9268  48.894 0.14229  0.094291 .  
Donor_ID        10   214.778 21.4778 105.789 0.76963 9.999e-05 ***
Site             2     1.436  0.7181   3.537 0.00515  0.006799 ** 
Residuals      114    23.145  0.2030         0.08294 9.999e-05 ***
Total          130   279.066                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## Luminal Colon
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   4    28.713  7.1783  80.888 0.10222    0.3610    
Donor_ID         9   241.950 26.8833 302.934 0.86132 9.999e-05 ***
Site             2     0.125  0.0626   0.705 0.00045    0.5481    
Residuals      114    10.117  0.0887         0.03601 9.999e-05 ***
Total          129   280.905                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

No repeated measures 
```R
adonis2(formula = data.dist ~ Sequencing_Run + Sex + Site, data = metadata, permutations = 10000)
                Df SumOfSqs      R2       F    Pr(>F)    
Sequencing_Run   4   28.713 0.10222  3.8071    0.0003 ***
Sex              1   21.916 0.07802 11.6235 9.999e-05 ***
Site             2    0.246 0.00088  0.0653    0.9976    
Residual       122  230.029 0.81889  
```
## Luminal SI
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   5    49.535  9.9070 10.1667 0.17847  0.008799 ** 
Donor_ID         9   113.768 12.6409 12.9723 0.40989 9.999e-05 ***
Site             2     2.192  1.0958  1.1245 0.00790  0.267773    
Residuals      115   112.062  0.9745         0.40374 9.999e-05 ***
Total          131   277.557                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
```

## Luminal

no Repeated measures 
```R
adonis2(formula = data.dist ~ Sequencing_Run + Sex + Site_General, data = metadata, permutations = 10000)
                Df SumOfSqs      R2      F    Pr(>F)    
Sequencing_Run   5    34.16 0.06237  4.695 9.999e-05 ***
Sex              1    33.38 0.06094 22.936 9.999e-05 ***
Site_General     1   110.55 0.20185 75.972 9.999e-05 ***
Residual       254   369.62 0.67484                     
Total          261   547.71 1.00000 
```

```R
adonis2(formula = data.dist ~ Sequencing_Run + Sex + Donor_ID + Site_General, data = metadata, permutations = 10000)
                Df SumOfSqs      R2       F    Pr(>F)    
Sequencing_Run   5    34.16 0.06237  12.646 9.999e-05 ***
Sex              1    33.38 0.06094  61.779 9.999e-05 ***
Donor_ID         8   237.78 0.43413  55.014 9.999e-05 ***
Site_General     1   109.49 0.19990 202.652 9.999e-05 ***
Residual       246   132.91 0.24266                      
Total          261   547.71 1.00000 
```

## Mucosal
no repeated measures
```R
adonis2(formula = data.dist ~ Sequencing_Run + Sex + Site_General, data = metadata, permutations = 10000)
                Df SumOfSqs      R2       F    Pr(>F)    
Sequencing_Run   5    85.43 0.14879  9.8976 9.999e-05 ***
Sex              1    22.45 0.03911 13.0068 9.999e-05 ***
Site_General     1    15.71 0.02736  9.0991     2e-04 ***
Residual       261   450.58 0.78474                      
Total          268   574.18 1.00000 
```

```R
adonis2(formula = data.dist ~ Sequencing_Run + Sex + Donor_ID + Site_General, data = metadata, permutations = 10000)
                Df SumOfSqs      R2      F    Pr(>F)    
Sequencing_Run   5    85.43 0.14879 32.311 9.999e-05 ***
Sex              1    22.45 0.03911 42.461 9.999e-05 ***
Donor_ID        10   319.10 0.55575 60.341 9.999e-05 ***
Site_General     1    14.46 0.02518 27.338 9.999e-05 ***
Residual       251   132.73 0.23117                     
Total          268   574.18 1.00000            
```