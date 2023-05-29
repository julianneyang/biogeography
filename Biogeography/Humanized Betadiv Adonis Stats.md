# Full Dataset
### Microbiota* Site General
```R
                       Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                         Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run            2     51.20  25.600  19.590 0.12358 9.999e-05 ***
Sex                       1      3.31   3.312   2.535 0.00799   0.07359 .  
Type                      1      1.86   1.865   1.427 0.00450   0.23648    
Microbiota                1     90.58  90.576  69.313 0.21863 9.999e-05 ***
Site_General              1      7.08   7.075   5.414 0.01708   0.00310 ** 
Microbiota:Site_General   1     22.43  22.435  17.168 0.05415 9.999e-05 ***
Residuals               182    237.83   1.307         0.57406              
Total                   189    414.29                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
```

### Microbiota* Site
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2     51.20  25.600  19.273 0.12358 9.999e-05 ***
Sex               1      3.31   3.312   2.493 0.00799   0.07489 .  
Type              1      1.86   1.865   1.404 0.00450   0.24058    
Microbiota        1     90.58  90.576  68.190 0.21863 9.999e-05 ***
Site              5     10.50   2.100   1.581 0.02534   0.10069    
Microbiota:Site   5     25.72   5.144   3.873 0.06208   0.00020 ***
Residuals       174    231.12   1.328         0.55787              
Total           189    414.29                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Microbiota * Type
```R
Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2     51.20  25.600  17.757 0.12358 9.999e-05 ***
Sex               1      3.31   3.312   2.297 0.00799   0.09469 .  
Site              5     11.58   2.317   1.607 0.02796   0.09129 .  
Type              1      1.54   1.541   1.069 0.00372   0.34707    
Microbiota        1     89.81  89.814  62.300 0.21679 9.999e-05 ***
Type:Microbiota   1      0.23   0.231   0.160 0.00056   0.89191    
Residuals       178    256.61   1.442         0.61939              
Total           189    414.29                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
---
# Site Subsets
### Luminal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2    14.247   7.124   6.591 0.12065   0.00020 ***
Sex              1     3.112   3.112   2.879 0.02635   0.05609 .  
Microbiota       1    51.047  51.047  47.230 0.43228 9.999e-05 ***
Site             2     0.397   0.199   0.184 0.00337   0.96570    
Microbiota:Site  2     0.648   0.324   0.300 0.00549   0.90221    
Residuals       45    48.637   1.081         0.41187              
Total           53   118.090                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Luminal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1     2.449  2.4490  1.4476 0.02344 0.2396760    
Sex              1    12.045 12.0452  7.1201 0.11531 0.0005999 ***
Microbiota       1    13.928 13.9277  8.2329 0.13333 0.0002000 ***
Site             2     1.384  0.6920  0.4091 0.01325 0.8271173    
Microbiota:Site  2     5.294  2.6471  1.5647 0.05068 0.1865813    
Residuals       41    69.361  1.6917         0.66399              
Total           48   104.461                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2    20.931  10.465  8.4494 0.19255 9.999e-05 ***
Sex              1     1.448   1.448  1.1693 0.01332    0.3079    
Microbiota       1    34.028  34.028 27.4732 0.31303 9.999e-05 ***
Site             2     1.550   0.775  0.6259 0.01426    0.6593    
Microbiota:Site  2     1.204   0.602  0.4862 0.01108    0.7584    
Residuals       40    49.543   1.239         0.45576              
Total           48   108.705                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal SI 
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2    25.922 12.9612  9.4652 0.29312 9.999e-05 ***
Sex              1     1.963  1.9635  1.4339 0.02220    0.2431    
Microbiota       1    15.905 15.9048 11.6149 0.17985 9.999e-05 ***
Site             2     4.189  2.0947  1.5297 0.04737    0.1982    
Microbiota:Site  2     0.745  0.3723  0.2719 0.00842    0.9344    
Residuals       29    39.711  1.3693         0.44904              
Total           37    88.436                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

# Type Subsets
### Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2    28.246  14.123  13.452 0.12614 9.999e-05 ***
Sex               1     5.054   5.054   4.814 0.02257  0.005799 ** 
Site              2     1.404   0.702   0.669 0.00627  0.629837    
Microbiota        1    86.140  86.140  82.050 0.38468 9.999e-05 ***
Type              1     0.588   0.588   0.560 0.00263  0.590541    
Microbiota:Type   1     3.810   3.810   3.629 0.01701  0.023098 *  
Residuals        94    98.685   1.050         0.44070              
Total           102   223.926                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2    10.983  5.4914  3.2066 0.05830   0.01160 *  
Sex              1    12.519 12.5186  7.3101 0.06646   0.00050 ***
Site             2     3.270  1.6352  0.9549 0.01736   0.44986    
Microbiota       1    23.074 23.0742 13.4739 0.12249 9.999e-05 ***
Type             1     0.927  0.9267  0.5411 0.00492   0.60854    
Microbiota:Type  1     4.022  4.0219  2.3486 0.02135   0.09169 .  
Residuals       78   133.575  1.7125         0.70911              
Total           86   188.370                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Distal Colon
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2    16.056  8.0281  6.6146 0.22446    0.0002 ***
Sex              1     1.233  1.2326  1.0156 0.01723    0.3774    
Microbiota       1    20.883 20.8825 17.2057 0.29193 9.999e-05 ***
Type             1     1.444  1.4441  1.1899 0.02019    0.3178    
Microbiota:Type  1     1.575  1.5749  1.2976 0.02202    0.2717    
Residuals       25    30.342  1.2137         0.42418              
Total           31    71.533                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### PC
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2     8.014   4.007  3.3701 0.10244    0.0112 *  
Sex              1     1.325   1.325  1.1142 0.01693    0.3424    
Microbiota       1    33.836  33.836 28.4585 0.43250 9.999e-05 ***
Type             1     0.052   0.052  0.0434 0.00066    0.9805    
Microbiota:Type  1     1.717   1.717  1.4442 0.02195    0.2401    
Residuals       28    33.291   1.189         0.42553              
Total           34    78.235                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Cec
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2     9.343   4.672  4.1015 0.11720   0.00350 ** 
Sex              1     2.799   2.799  2.4577 0.03511   0.09139 .  
Microbiota       1    33.557  33.557 29.4620 0.42094 9.999e-05 ***
Type             1     0.395   0.395  0.3469 0.00496   0.74613    
Microbiota:Type  1     0.594   0.594  0.5211 0.00745   0.61814    
Residuals       29    33.030   1.139         0.41434              
Total           35    79.718                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1```

### Ile
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sequencing_Run   2     4.144  2.0721  1.1951 0.06078 0.3155    
Sex              1     5.543  5.5428  3.1966 0.08128 0.0400 *  
Microbiota       1    13.678 13.6779  7.8884 0.20058 0.0004 ***
Type             1     2.947  2.9473  1.6998 0.04322 0.1823    
Microbiota:Type  1     0.264  0.2637  0.1521 0.00387 0.9149    
Residuals       24    41.614  1.7339         0.61027           
Total           30    68.190                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Jej
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sequencing_Run   2     6.941  3.4703  1.7184 0.10311 0.1447  
Sex              1     3.131  3.1312  1.5506 0.04652 0.2167  
Microbiota       1     6.522  6.5224  3.2299 0.09690 0.0291 *
Type             1     0.480  0.4802  0.2378 0.00713 0.8414  
Microbiota:Type  1     3.792  3.7919  1.8777 0.05633 0.1490  
Residuals       23    46.447  2.0194         0.69001         
Total           29    67.313                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Duo
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)  
Sequencing_Run   2     3.998  1.9990  1.0654 0.07166 0.38196  
Sex              1     7.654  7.6541  4.0792 0.13719 0.01850 *
Microbiota       1     4.490  4.4901  2.3929 0.08048 0.09269 .
Type             1     1.041  1.0413  0.5549 0.01866 0.59974  
Microbiota:Type  1     2.957  2.9571  1.5759 0.05300 0.22108  
Residuals       19    35.651  1.8764         0.63901          
Total           25    55.792                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

# Source Subset: Humanized vs SPF
### Humanized
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    41.429  20.715  18.786 0.20774 9.999e-05 ***
Sex             1    38.016  38.016  34.476 0.19062 9.999e-05 ***
Site            5    26.730   5.346   4.848 0.13403    0.0002 ***
Type            1     0.629   0.629   0.571 0.00316    0.5785    
Residuals      84    92.623   1.103         0.46444              
Total          93   199.427                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    41.429  20.715  18.786 0.20774 9.999e-05 ***
Sex             1    38.016  38.016  34.476 0.19062 9.999e-05 ***
Type            1     1.434   1.434   1.301 0.00719    0.2833    
Site            5    25.925   5.185   4.702 0.13000 9.999e-05 ***
Residuals      84    92.623   1.103         0.46444              
Total          93   199.427                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    41.429  20.715  19.124 0.20774 9.999e-05 ***
Sex             1    38.016  38.016  35.097 0.19062 9.999e-05 ***
Type            1     1.434   1.434   1.324 0.00719    0.2596    
Site_General    1    23.231  23.231  21.448 0.11649 9.999e-05 ***
Residuals      88    95.317   1.083         0.47796              
Total          93   199.427                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Cedars SPF 
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    29.672 14.8362  9.1307 0.14204 9.999e-05 ***
Sex             1     1.070  1.0704  0.6588 0.00512    0.5376    
Site            5    37.219  7.4439  4.5812 0.17816 9.999e-05 ***
Type            1     1.203  1.2034  0.7406 0.00576    0.5034    
Residuals      86   139.740  1.6249         0.66891              
Total          95   208.905                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    29.672 14.8362  9.1307 0.14204 9.999e-05 ***
Sex             1     1.070  1.0704  0.6588 0.00512    0.5455    
Type            1     1.309  1.3087  0.8054 0.00626    0.4632    
Site            5    37.114  7.4228  4.5682 0.17766 9.999e-05 ***
Residuals      86   139.740  1.6249         0.66891              
Total          95   208.905                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    29.672 14.8362  9.1679 0.14204 9.999e-05 ***
Sex             1     1.070  1.0704  0.6614 0.00512    0.5387    
Type            1     1.309  1.3087  0.8087 0.00626    0.4497    
Site_General    1    31.209 31.2091 19.2854 0.14939 9.999e-05 ***
Residuals      90   145.645  1.6183         0.69718              
Total          95   208.905                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


```

# Humanized Dataset broken down into Site and Type
### Luminal
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    12.418  6.2088  5.2546 0.11337 0.0009999 ***
Sex             1    23.945 23.9452 20.2651 0.21861 9.999e-05 ***
Site_General    1    17.636 17.6357 14.9253 0.16101 9.999e-05 ***
Residuals      47    55.535  1.1816         0.50701              
Total          51   109.533                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2    12.418  6.2088  5.6844 0.11337  0.000300 ***
Sex               1    23.945 23.9452 21.9226 0.21861 9.999e-05 ***
Site_General      1    17.636 17.6357 16.1461 0.16101 9.999e-05 ***
Sex:Site_General  1     5.291  5.2911  4.8442 0.04831  0.009099 ** 
Residuals        46    50.244  1.0923         0.45871              
Total            51   109.533                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    27.694 13.8469  13.974 0.31064 9.999e-05 ***
Sex             1    21.343 21.3427  21.538 0.23940 9.999e-05 ***
Site_General    1     3.449  3.4494   3.481 0.03869    0.0269 *  
Residuals      37    36.664  0.9909         0.41126              
Total          41    89.150                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2    27.694 13.8469 13.8483 0.31064 9.999e-05 ***
Sex               1    21.343 21.3427 21.3448 0.23940 9.999e-05 ***
Site_General      1     3.449  3.4494  3.4497 0.03869    0.0293 *  
Sex:Site_General  1     0.667  0.6675  0.6676 0.00749    0.5356    
Residuals        36    35.996  0.9999         0.40377              
Total            41    89.150                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Luminal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    19.716  9.8581  7.8168 0.36118 9.999e-05 ***
Sex             1     8.234  8.2343  6.5292 0.15085    0.0015 ** 
Site            2     0.153  0.0765  0.0607 0.00280    0.9945    
Residuals      21    26.484  1.2611         0.48517              
Total          26    54.587                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Luminal SI
```R
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  1     0.280  0.2796  0.1825 0.00531     0.863    
Sex             1    20.315 20.3151 13.2607 0.38604 9.999e-05 ***
Site            2     1.390  0.6951  0.4537 0.02642     0.798    
Residuals      20    30.640  1.5320         0.58223              
Total          24    52.625                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

### Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    20.600 10.2999 13.7538 0.43426    0.0003 ***
Sex             1    12.575 12.5746 16.7914 0.26509 9.999e-05 ***
Site            2     1.531  0.7655  1.0222 0.03227    0.4363    
Residuals      17    12.731  0.7489         0.26838              
Total          22    47.436                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

```

### Mucosal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sequencing_Run  2    17.604  8.8018  9.8228 0.43662 0.0002 ***
Sex             1     9.083  9.0833 10.1370 0.22529 0.0004 ***
Site            2     1.982  0.9912  1.1062 0.04917 0.3776    
Residuals      13    11.649  0.8961         0.28892           
Total          18    40.318                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Colon
```R
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    19.482  9.7408  8.2085 0.19160 9.999e-05 ***
Sex             1    24.998 24.9981 21.0657 0.24585 9.999e-05 ***
Site            2     1.212  0.6061  0.5108 0.01192    0.7453    
Type            1     4.959  4.9594  4.1792 0.04878    0.0204 *  
Residuals      43    51.027  1.1867         0.50185              
Total          49   101.678                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    14.288  7.1440  4.8815 0.15223    0.0012 ** 
Sex             1    21.494 21.4941 14.6871 0.22900 9.999e-05 ***
Site            2     2.860  1.4299  0.9770 0.03047    0.4397    
Type            1     1.070  1.0696  0.7309 0.01140    0.5000    
Residuals      37    54.148  1.4635         0.57691              
Total          43    93.860                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

# Cedars SPF broken down into Site and Type

### Luminal
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    15.001  7.5007  4.6652 0.13401    0.0013 ** 
Sex             1     1.308  1.3080  0.8135 0.01168    0.4552    
Site_General    1    21.676 21.6764 13.4821 0.19364 9.999e-05 ***
Residuals      46    73.959  1.6078         0.66067              
Total          50   111.945                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2    15.001  7.5007  4.6378 0.13401 0.0008999 ***
Sex               1     1.308  1.3080  0.8088 0.01168 0.4590541    
Site_General      1    21.676 21.6764 13.4029 0.19364 9.999e-05 ***
Sex:Site_General  1     1.180  1.1803  0.7298 0.01054 0.5020498    
Residuals        45    72.778  1.6173         0.65013              
Total            50   111.945                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    23.470 11.7349  7.6226 0.24038 9.999e-05 ***
Sex             1     0.513  0.5134  0.3335 0.00526    0.7622    
Site_General    1    12.072 12.0721  7.8417 0.12365    0.0005 ***
Residuals      40    61.579  1.5395         0.63071              
Total          44    97.635                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                 Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run    2    23.470 11.7349  8.5198 0.24038 9.999e-05 ***
Sex               1     0.513  0.5134  0.3728 0.00526    0.7315    
Site_General      1    12.072 12.0721  8.7646 0.12365    0.0003 ***
Sex:Site_General  1     7.862  7.8616  5.7077 0.08052    0.0019 ** 
Residuals        39    53.718  1.3774         0.55019              
Total            44    97.635                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Luminal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sequencing_Run  2     8.007  4.0033 2.03565 0.14584 0.1017
Sex             1     4.532  4.5320 2.30449 0.08255 0.1107
Site            2     1.062  0.5309 0.26995 0.01934 0.9062
Residuals      21    41.299  1.9666         0.75227       
Total          26    54.899                 1.00000  
```

### Luminal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)
Sequencing_Run  1     0.691  0.6913 0.31271 0.01319 0.8010
Sex             1     4.627  4.6269 2.09301 0.08826 0.1244
Site            2     5.103  2.5514 1.15415 0.09734 0.3445
Residuals      19    42.002  2.2106         0.80121       
Total          23    52.423                 1.00000  
```

### Mucosal Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Sequencing_Run  2    14.371  7.1855  4.0632 0.26829 0.0027 **
Sex             1     3.051  3.0513  1.7255 0.05696 0.1741   
Site            2     0.775  0.3874  0.2191 0.01446 0.9398   
Residuals      20    35.369  1.7684         0.66028          
Total          25    53.566                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Mucosal SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Sequencing_Run  2    15.472  7.7359  4.8610 0.37477 0.0014 **
Sex             1     3.025  3.0249  1.9008 0.07327 0.1544   
Site            2     2.098  1.0491  0.6592 0.05082 0.6718   
Residuals      13    20.689  1.5914         0.50113          
Total          18    41.284                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### Colon
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run  2    18.474  9.2369  6.1216 0.16977 9.999e-05 ***
Sex             1     6.774  6.7744  4.4897 0.06226   0.01250 *  
Site            2     6.182  3.0910  2.0485 0.05681   0.08719 .  
Type            1     7.978  7.9778  5.2872 0.07331   0.00420 ** 
Residuals      46    69.410  1.5089         0.63785              
Total          52   108.818                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```

### SI
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)  
Sequencing_Run  2    11.002  5.5012  2.8514 0.11673 0.0231 *
Sex             1     5.932  5.9320  3.0747 0.06294 0.0409 *
Site            2     5.763  2.8817  1.4936 0.06115 0.2022  
Type            1     2.102  2.1023  1.0897 0.02230 0.3394  
Residuals      36    69.455  1.9293         0.73688         
Total          42    94.255                 1.00000         
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```