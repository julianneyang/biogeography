
# Site RPCA 
## Repeated Measures PERMANOVA
### Luminal
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    43.405  43.405  85.969 0.35295 0.00009999 ***
Sex             1     6.247   6.247  12.372 0.05080     0.3706    
Site_General    1    45.051  45.051  89.228 0.36634 0.00009999 ***
Residuals      56    28.274   0.505         0.22991 0.00009999 ***
Total          59   122.977                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    36.380  36.380 109.969 0.37052     0.0177 *  
Sex             1     7.760   7.760  23.457 0.07903     0.3195    
Site_General    1    39.491  39.491 119.371 0.40220 0.00009999 ***
Residuals      44    14.556   0.331         0.14825 0.00009999 ***
Total          47    98.187                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    19.088 19.0877 14.3256 0.39171 0.0006999 ***
Sex             1     2.776  2.7755  2.0831 0.05696 0.6544346    
Site            2     1.550  0.7749  0.5815 0.03180 0.3372663    
Residuals      19    25.316  1.3324         0.51953 0.0245975 *  
Total          23    48.729                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    21.126 21.1262 24.3919 0.43939 0.00009999 ***
Sex             1     8.084  8.0839  9.3335 0.16813     0.2843    
Site            2     2.415  1.2074  1.3941 0.05022     0.0278 *  
Residuals      19    16.456  0.8661         0.34226     0.0025 ** 
Total          23    48.081                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
Sequencing_Run  1    20.614 20.6144 14.4106 0.33565 0.005499 **
Sex             1     3.701  3.7012  2.5874 0.06026 0.532147   
Site            2     1.338  0.6688  0.4676 0.02178 0.443556   
Residuals      25    35.762  1.4305         0.58230 0.013399 * 
Total          29    61.416                 1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    22.665 22.6646 17.2488 0.37332 0.00009999 ***
Sex             1     3.113  3.1131  2.3692 0.05128    0.63314    
Site            2     2.083  1.0415  0.7927 0.03431    0.07829 .  
Residuals      25    32.850  1.3140         0.54109    0.01240 *  
Total          29    60.710                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
----
## PERMANOVA
### Luminal Subset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    43.405  43.405  85.969 0.35295 9.999e-05 ***
Sex             1     6.247   6.247  12.372 0.05080 9.999e-05 ***
Site_General    1    45.051  45.051  89.228 0.36634 9.999e-05 ***
Residuals      56    28.274   0.505         0.22991              
Total          59   122.977                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    43.405  43.405  83.352 0.35295 9.999e-05 ***
Sex             1     6.247   6.247  11.996 0.05080 9.999e-05 ***
Site            5    46.246   9.249  17.762 0.37606 9.999e-05 ***
Residuals      52    27.079   0.521         0.22019              
Total          59   122.977                 1.00000
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    1    43.405  43.405  86.245 0.35295 9.999e-05 ***
Sex               1     6.247   6.247  12.412 0.05080 9.999e-05 ***
Site_General      1    45.051  45.051  89.515 0.36634 9.999e-05 ***
Sex:Site_General  1     0.594   0.594   1.180 0.00483    0.3188    
Residuals        55    27.680   0.503         0.22509              
Total            59   122.977                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal Subset 
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    36.380  36.380 109.969 0.37052 9.999e-05 ***
Sex             1     7.760   7.760  23.457 0.07903 9.999e-05 ***
Site_General    1    39.491  39.491 119.371 0.40220 9.999e-05 ***
Residuals      44    14.556   0.331         0.14825              
Total          47    98.187                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    36.380  36.380 106.336 0.37052 9.999e-05 ***
Sex             1     7.760   7.760  22.682 0.07903 9.999e-05 ***
Site            5    40.362   8.072  23.595 0.41107 9.999e-05 ***
Residuals      40    13.685   0.342         0.13938              
Total          47    98.187                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    1    36.380  36.380 110.046 0.37052 9.999e-05 ***
Sex               1     7.760   7.760  23.473 0.07903 9.999e-05 ***
Site_General      1    39.491  39.491 119.455 0.40220 9.999e-05 ***
Sex:Site_General  1     0.341   0.341   1.031 0.00347    0.3657    
Residuals        43    14.215   0.331         0.14478              
Total            47    98.187                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    22.665 22.6646 17.2488 0.37332 9.999e-05 ***
Sex             1     3.113  3.1131  2.3692 0.05128   0.09989 .  
Site            2     2.083  1.0415  0.7927 0.03431   0.53875    
Residuals      25    32.850  1.3140         0.54109              
Total          29    60.710                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    20.614 20.6144 14.4106 0.33565 9.999e-05 ***
Sex             1     3.701  3.7012  2.5874 0.06026   0.08509 .  
Site            2     1.338  0.6688  0.4676 0.02178   0.77822    
Residuals      25    35.762  1.4305         0.58230              
Total          29    61.416                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    21.126 21.1262 24.3919 0.43939 9.999e-05 ***
Sex             1     8.084  8.0839  9.3335 0.16813    0.0003 ***
Site            2     2.415  1.2074  1.3941 0.05022    0.2533    
Residuals      19    16.456  0.8661         0.34226              
Total          23    48.081                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    19.088 19.0877 14.3256 0.39171 9.999e-05 ***
Sex             1     2.776  2.7755  2.0831 0.05696     0.138    
Site            2     1.550  0.7749  0.5815 0.03180     0.683    
Residuals      19    25.316  1.3324         0.51953              
Total          23    48.729                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


# Type RPCA
## Repeated Measures PERMANOVA
### Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)   
Sequencing_Run  1    42.210  42.210  39.828 0.38819 0.00120 **
Sex             1    10.939  10.939  10.321 0.10060 0.38696   
Site            2     2.814   1.407   1.328 0.02588 0.04820 * 
Type            1     1.901   1.901   1.794 0.01748 0.07859 . 
Residuals      48    50.872   1.060         0.46785 0.00250 **
Total          53   108.735                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### SI
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    37.973  37.973 27.4174 0.34255 0.00009999 ***
Sex             1     5.620   5.620  4.0582 0.05070     0.6101    
Site            2     0.447   0.223  0.1613 0.00403     0.8208    
Type            1     0.334   0.334  0.2410 0.00301     0.7418    
Residuals      48    66.479   1.385         0.59971     0.0192 *  
Total          53   110.853                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
-----
## PERMANOVA
### Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    42.210  42.210  39.828 0.38819 9.999e-05 ***
Sex             1    10.939  10.939  10.321 0.10060    0.0002 ***
Site            2     2.814   1.407   1.328 0.02588    0.2622    
Type            1     1.901   1.901   1.794 0.01748    0.1694    
Residuals      48    50.872   1.060         0.46785              
Total          53   108.735                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    37.973  37.973 27.4174 0.34255 9.999e-05 ***
Sex             1     5.620   5.620  4.0582 0.05070    0.0199 *  
Site            2     0.447   0.223  0.1613 0.00403    0.9640    
Type            1     0.334   0.334  0.2410 0.00301    0.7951    
Residuals      48    66.479   1.385         0.59971              
Total          53   110.853                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Duodenum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    13.138 13.1380  9.3890 0.35761 0.00009999 ***
Sex             1     2.241  2.2409  1.6014 0.06100    0.62964    
Type            1     1.770  1.7697  1.2647 0.04817    0.05529 .  
Residuals      14    19.590  1.3993         0.53323    0.01660 *  
Total          17    36.739                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Jejunum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    12.707 12.7072  9.5530 0.34778 0.00009999 ***
Sex             1     5.203  5.2030  3.9115 0.14240     0.2610    
Type            1     0.005  0.0055  0.0041 0.00015     0.9977    
Residuals      14    18.623  1.3302         0.50967     0.0142 *  
Total          17    36.538                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Ileum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2   Pr(>F)   
Sequencing_Run  1    12.034 12.0342  9.6321 0.33012 0.009999 **
Sex             1     6.362  6.3622  5.0922 0.17452 0.191781   
Type            1     0.567  0.5666  0.4535 0.01554 0.480152   
Residuals      14    17.491  1.2494         0.47982 0.006999 **
Total          17    36.454                 1.00000            
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```


### Cecum
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1    12.023 12.0227 12.5737 0.32784 0.0076992 ** 
Sex             1     6.581  6.5808  6.8824 0.17945 0.1600840    
Type            1     4.682  4.6821  4.8967 0.12768 0.0064994 ** 
Residuals      14    13.386  0.9562         0.36503 0.0009999 ***
Total          17    36.672                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Proximal_Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    15.278 15.2780 14.7425 0.42314 0.00009999 ***
Sex             1     4.724  4.7243  4.5587 0.13084     0.3167    
Type            1     1.596  1.5957  1.5398 0.04419     0.1083    
Residuals      14    14.509  1.0363         0.40183     0.0016 ** 
Total          17    36.106                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Distal_Colon
```R
Permutation: free
Number of permutations: 0

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run  1    16.099 16.0991 11.9653 0.44394 0.00009999 ***
Sex             1     0.934  0.9338  0.6941 0.02575     0.7985    
Type            1     0.394  0.3941  0.2929 0.01087     0.6542    
Residuals      14    18.837  1.3455         0.51944     0.0105 *  
Total          17    36.264                 1.00000               
---
Signif. codes:  0 ‘
```