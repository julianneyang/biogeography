## Mucosal SI
```R
adonis2(formula = D ~ ., data = mtdat[, metadata_order, drop = F], permutations = 0)
                Df SumOfSqs      R2       F Pr(>F) aov.tab
Sequencing_Run   5   97.484 0.32693 13.1486      1        
Sex              1    7.830 0.02626  5.2806      1        
Site             2    1.589 0.00533  0.5356      1        
Residual       129  191.282 0.64149                       
Total          137  298.185 1.00000
```

## Mucosal Colon
```R
adonis2(formula = D ~ ., data = mtdat[, metadata_order, drop = F], permutations = 0)
                Df SumOfSqs      R2      F Pr(>F) aov.tab
Sequencing_Run   4   39.707 0.14229 5.5069      1        
Sex              1   16.518 0.05919 9.1633      1        
Site             2    1.119 0.00401 0.3104      1        
Residual       123  221.722 0.79451                      
Total          130  279.066 1.00000 
```

## Luminal Colon
```R
adonis2(formula = D ~ ., data = mtdat[, metadata_order, drop = F], permutations = 0)
                Df SumOfSqs      R2       F Pr(>F) aov.tab
Sequencing_Run   4   28.713 0.10222  3.8071      1        
Sex              1   21.916 0.07802 11.6235      1        
Site             2    0.246 0.00088  0.0653      1        
Residual       122  230.029 0.81889                       
Total          129  280.905 1.00000  
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
adonis2(formula = D ~ ., data = mtdat[, metadata_order, drop = F], permutations = 0)
                Df SumOfSqs      R2      F Pr(>F) aov.tab
Sequencing_Run   5   49.535 0.17847 5.7814      1        
Sex              1   15.059 0.05425 8.7877      1        
Site             2    2.192 0.00790 0.6395      1        
Residual       123  210.772 0.75938                      
Total          131  277.557 1.00000  
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