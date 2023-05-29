# WT Validation Cohort
## Before CombatSeq
## PERMANOVA, Repeated Measures
### Mucosal
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run   1     0.943   0.943   0.681 0.00338     0.7112    
Sex              1     4.211   4.211   3.042 0.01508     0.2494    
Site_General     1    92.763  92.763  67.012 0.33218 0.00009999 ***
Residuals      131   181.339   1.384         0.64937 0.00009999 ***
Total          134   279.256                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    35.075  35.075  35.472 0.25510 0.00009999 ***
Sex             1     0.387   0.387   0.391 0.00281     0.8987    
Site            2    39.741  19.871  20.095 0.28903 0.00009999 ***
Residuals      63    62.296   0.989         0.45306 0.00009999 ***
Total          67   137.499                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sequencing_Run  1     1.545  1.5445  0.8068 0.01091 0.6800    
Sex             1     3.093  3.0931  1.6157 0.02185 0.4578    
Site            2    18.201  9.1004  4.7536 0.12860 0.0002 ***
Residuals      62   118.694  1.9144         0.83864 0.0375 *  
Total          66   141.533                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
-----
## PERMANOVA
### Mucosal Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1     0.943   0.943   0.681 0.00338    0.5159    
Sex              1     4.211   4.211   3.042 0.01508    0.0459 *  
Site_General     1    92.763  92.763  67.012 0.33218 9.999e-05 ***
Residuals      131   181.339   1.384         0.64937              
Total          134   279.256                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                  Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run     1     0.943   0.943   0.678 0.00338    0.5087    
Sex                1     4.211   4.211   3.029 0.01508    0.0470 *  
Site_General       1    92.763  92.763  66.731 0.33218 9.999e-05 ***
Sex:Site_General   1     0.626   0.626   0.450 0.00224    0.6451    
Residuals        130   180.713   1.390         0.64712              
Total            134   279.256                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal Colon
```R
Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    35.075  35.075  35.472 0.25510 9.999e-05 ***
Sex             1     0.387   0.387   0.391 0.00281    0.6833    
Site            2    39.741  19.871  20.095 0.28903 9.999e-05 ***
Residuals      63    62.296   0.989         0.45306              
Total          67   137.499                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    35.075  35.075  34.989 0.25510 9.999e-05 ***
Sex             1     0.387   0.387   0.386 0.00281    0.6828    
Site            2    39.741  19.871  19.822 0.28903 9.999e-05 ***
Sex:Site        2     1.146   0.573   0.572 0.00833    0.6843    
Residuals      61    61.150   1.002         0.44473              
Total          67   137.499                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1     1.545  1.5445  0.8068 0.01091 0.4535546    
Sex             1     3.093  3.0931  1.6157 0.02185 0.1976802    
Site            2    18.201  9.1004  4.7536 0.12860 0.0008999 ***
Residuals      62   118.694  1.9144         0.83864              
Total          66   141.533                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
## After CombatSeq
### Mucosal
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1     5.499   5.499   3.970 0.01953    0.0163 *  
Sex              1     4.435   4.435   3.201 0.01575    0.0387 *  
Site_General     1    88.816  88.816  64.110 0.31538 9.999e-05 ***
Residuals      132   182.868   1.385         0.64935              
Total          135   281.618                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    61.042  61.042  79.796 0.44347 9.999e-05 ***
Sex             1     0.065   0.065   0.085 0.00047    0.9293    
Site            2    28.345  14.172  18.526 0.20592 9.999e-05 ***
Residuals      63    48.194   0.765         0.35013              
Total          67   137.646                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)   
Sequencing_Run  1     1.986  1.9857  1.0291 0.01380 0.35976   
Sex             1     4.566  4.5661  2.3663 0.03172 0.09069 . 
Site            2    15.813  7.9063  4.0974 0.10986 0.00280 **
Residuals      63   121.564  1.9296         0.84461           
Total          67   143.928                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)   
Sequencing_Run  1     1.986  1.9857  1.0402 0.01380 0.36236   
Sex             1     4.566  4.5661  2.3920 0.03172 0.09289 . 
Site            2    15.813  7.9063  4.1418 0.10986 0.00200 **
Sex:Site        2     5.120  2.5600  1.3411 0.03557 0.25347   
Residuals      61   116.444  1.9089         0.80904           
Total          67   143.928                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

## PERMANOVA, repeated Measures
### Mucosal Dataset
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run   1     5.499   5.499   3.970 0.01953     0.1973    
Sex              1     4.435   4.435   3.201 0.01575     0.2704    
Site_General     1    88.816  88.816  64.110 0.31538 0.00009999 ***
Residuals      132   182.868   1.385         0.64935 0.00009999 ***
Total          135   281.618                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    61.042  61.042  79.796 0.44347 0.00009999 ***
Sex             1     0.065   0.065   0.085 0.00047     0.9906    
Site            2    28.345  14.172  18.526 0.20592 0.00009999 ***
Residuals      63    48.194   0.765         0.35013 0.00009999 ***
Total          67   137.646                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sequencing_Run  1     1.986  1.9857  1.0291 0.01380 0.6396    
Sex             1     4.566  4.5661  2.3663 0.03172 0.3515    
Site            2    15.813  7.9063  4.0974 0.10986 0.0003 ***
Residuals      63   121.564  1.9296         0.84461 0.0493 *  
Total          67   143.928                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```