## Repeated Measures PERMANOVA
### Luminal 
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1     6.720  6.7196  5.7328 0.06135   0.407059    
Sex             1    23.945 23.9452 20.4288 0.21861   0.005399 ** 
Site_General    1    22.607 22.6066 19.2869 0.20639 0.00009999 ***
Residuals      48    56.262  1.1721         0.51365 0.00009999 ***
Total          51   109.533                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal
```r
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sequencing_Run  2    20.855 10.4275 10.5299 0.23393 0.2345    
Sex             1    25.763 25.7630 26.0160 0.28898 0.0013 ** 
Site_General    1     5.892  5.8917  5.9496 0.06609 0.0048 ** 
Residuals      37    36.640  0.9903         0.41100 0.0003 ***
Total          41    89.150                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sequencing_Run  2    17.101  8.5505 11.0152 0.42416 0.1261  
Sex             1    10.826 10.8257 13.9462 0.26851 0.1110  
Site            2     2.300  1.1500  1.4815 0.05705 0.1475  
Residuals      13    10.091  0.7762         0.25029 0.0130 *
Total          18    40.318                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
Sequencing_Run  2    13.008  6.5042  6.5230 0.27423 0.186881   
Sex             1    16.336 16.3357 16.3829 0.34437 0.009199 **
Site            2     1.141  0.5706  0.5723 0.02406 0.832017   
Residuals      17    16.951  0.9971         0.35734 0.005000 **
Total          22    47.436                 1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1     0.280  0.2796  0.1825 0.00531     0.8884    
Sex             1    20.315 20.3151 13.2607 0.38604 0.00009999 ***
Site            2     1.390  0.6951  0.4537 0.02642     0.5645    
Residuals      20    30.640  1.5320         0.58223     0.0343 *  
Total          24    52.625                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sequencing_Run  1    19.638 19.6377 16.2834 0.35975 0.0273 *
Sex             1     8.234  8.2343  6.8278 0.15085 0.2977  
Site            2     0.184  0.0918  0.0761 0.00336 0.9387  
Residuals      22    26.532  1.2060         0.48604 0.0174 *
Total          26    54.587                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sequencing_Run  1    17.374 17.3737 13.3925 0.17087 0.0284 *
Sex             1    24.760 24.7604 19.0866 0.24352 0.0213 *
Site            2     0.601  0.3003  0.2315 0.00591 0.8999  
Type            1     1.864  1.8637  1.4366 0.01833 0.2295  
Residuals      44    57.080  1.2973         0.56138 0.0257 *
Total          49   101.678                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
Sequencing_Run  1     4.970  4.9701  3.0950 0.05295 0.615038   
Sex             1    21.564 21.5636 13.4283 0.22974 0.009499 **
Site            2     1.985  0.9925  0.6181 0.02115 0.493451   
Type            1     4.319  4.3191  2.6896 0.04602 0.091391 . 
Residuals      38    61.022  1.6058         0.65014 0.028797 * 
Total          43    93.860                 1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Duodenum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sequencing_Run  1    0.7056  0.7056 0.28880 0.02815 0.7656
Sex             1    3.7803  3.7803 1.54722 0.15084 0.1455
Type            1    1.0301  1.0301 0.42161 0.04110 0.6961
Residuals       8   19.5465  2.4433         0.77991 0.5771
Total          11   25.0626                 1.00000
```

### Jejunum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
Sequencing_Run  1     0.923  0.9226  0.5775 0.02667 0.827817   
Sex             1     8.897  8.8975  5.5701 0.25721 0.009099 **
Type            1     5.604  5.6039  3.5082 0.16200 0.009099 **
Residuals      12    19.168  1.5974         0.55412 0.012399 * 
Total          15    34.592                 1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Ileum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
Sequencing_Run  1     3.869  3.8686  2.0733 0.11098 0.38876  
Sex             1     8.082  8.0824  4.3317 0.23186 0.06209 .
Type            1     0.517  0.5169  0.2770 0.01483 0.84212  
Residuals      12    22.390  1.8659         0.64233 0.12619  
Total          15    34.858                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Cecum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
Sequencing_Run  1     3.869  3.8686  2.0733 0.11098 0.38966  
Sex             1     8.082  8.0824  4.3317 0.23186 0.05909 .
Type            1     0.517  0.5169  0.2770 0.01483 0.84032  
Residuals      12    22.390  1.8659         0.64233 0.11749  
Total          15    34.858                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Proximal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
Sequencing_Run  1     7.100  7.1004  4.3341 0.19930 0.15158  
Sex             1     7.048  7.0481  4.3022 0.19783 0.19018  
Type            1     0.181  0.1806  0.1102 0.00507 0.78722  
Residuals      13    21.297  1.6383         0.59780 0.09209 .
Total          16    35.627                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’
```

### Distal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)    
Sequencing_Run  2    9.0957  4.5479  4.8066 0.29963 0.165083    
Sex             1    9.7571  9.7571 10.3122 0.32141 0.005799 ** 
Type            1    2.0424  2.0424  2.1586 0.06728 0.325667    
Residuals      10    9.4617  0.9462         0.31168 0.000200 ***
Total          14   30.3569                 1.00000             
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```