# Pre- any sort of Batch Correction 
### Mucosal Dataset
```R
> data.adonis=adonis(data.dist ~ Sequencing_Run + Sex + Genotype*Site, data=metadata, permutations=10000)
'adonis' will be deprecated: use 'adonis2' instead
> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3    275.40  91.801 161.994 0.29152 9.999e-05 ***
Sex              1      3.73   3.735   6.590 0.00395    0.0011 ** 
Genotype         2     27.90  13.950  24.616 0.02953 9.999e-05 ***
Site             5    364.54  72.907 128.654 0.38587 9.999e-05 ***
Genotype:Site   10     49.28   4.928   8.697 0.05217 9.999e-05 ***
Residuals      395    223.84   0.567         0.23695              
Total          416    944.70                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2     Pr(>F)    
Sequencing_Run   3    275.40  91.801 136.125 0.29152 0.00009999 ***
Sex              1      3.73   3.735   5.538 0.00395      0.003 ** 
Site             5    366.04  73.208 108.554 0.38746 0.00009999 ***
Genotype         2     26.40  13.199  19.572 0.02794 0.00009999 ***
Residuals      405    273.13   0.674         0.28911               
Total          416    944.70                 1.00000               
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
# Post- Batch Correction: Intersection(ASV) across Seq_Run, followed by ComBatSeq  
### Mucosal Dataset
```R
 Terms added sequentially (first to last)

                       Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run          3     44.34  14.780  13.301 0.05731 9.999e-05 ***
Sex                     1      2.31   2.307   2.076 0.00298    0.1234    
Genotype                2     28.69  14.347  12.912 0.03709 9.999e-05 ***
Site_General            1    287.81 287.814 259.020 0.37202 9.999e-05 ***
Genotype:Site_General   2      2.71   1.354   1.219 0.00350    0.3008    
Residuals             367    407.80   1.111         0.52710              
Total                 376    773.66                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3     44.34  14.780  15.116 0.05731 9.999e-05 ***
Sex              1      2.31   2.307   2.359 0.00298   0.08589 .  
Genotype         2     28.69  14.347  14.674 0.03709 9.999e-05 ***
Site             5    339.40  67.880  69.424 0.43869 9.999e-05 ***
Genotype:Site   10     11.82   1.182   1.209 0.01528   0.23658    
Residuals      355    347.10   0.978         0.44865              
Total          376    773.66                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Mucosal SI Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3     39.54 13.1799  7.8746 0.09846 9.999e-05 ***
Sex              1     15.20 15.1969  9.0797 0.03784    0.0002 ***
Genotype         2      9.68  4.8424  2.8932 0.02412    0.0212 *  
Site             2     21.40 10.7018  6.3940 0.05330    0.0002 ***
Genotype:Site    4      9.48  2.3691  1.4155 0.02360    0.1811    
Residuals      183    306.29  1.6737         0.76269              
Total          195    401.59                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

> data.adonis$aov.tab
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3     39.54 13.1799  7.8052 0.09846 9.999e-05 ***
Sex              1     15.20 15.1969  8.9997 0.03784 9.999e-05 ***
Site             2     20.07 10.0338  5.9421 0.04997    0.0003 ***
Genotype         2     11.02  5.5105  3.2633 0.02744    0.0109 *  
Residuals      187    315.77  1.6886         0.78629              
Total          195    401.59                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Permutation: free
Number of permutations: 10000
```

### Mucosal Colon Dataset
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3    130.96  43.652  43.690 0.35856 9.999e-05 ***
Sex              1      1.45   1.453   1.455 0.00398    0.2272    
Genotype         2     11.91   5.955   5.960 0.03261    0.0002 ***
Site             2     51.51  25.753  25.775 0.14102 9.999e-05 ***
Genotype:Site    4      1.55   0.388   0.388 0.00425    0.9273    
Residuals      168    167.85   0.999         0.45958              
Total          180    365.23                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   3    130.96  43.652  44.320 0.35856 9.999e-05 ***
Sex              1      1.45   1.453   1.476 0.00398    0.2376    
Site             2     50.15  25.074  25.457 0.13730 9.999e-05 ***
Genotype         2     13.27   6.634   6.736 0.03633    0.0003 ***
Residuals      172    169.41   0.985         0.46383              
Total          180    365.23                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### WT Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1    10.714  10.714   8.292 0.03828 0.0006999 ***
Sex              1     5.335   5.335   4.129 0.01906 0.0150985 *  
Site_General     1    91.992  91.992  71.192 0.32866 9.999e-05 ***
Residuals      133   171.860   1.292         0.61400              
Total          136   279.901                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1    10.714 10.7142  9.3770 0.03828 9.999e-05 ***
Sex              1     5.335  5.3351  4.6693 0.01906  0.009799 ** 
Site             5   116.456 23.2913 20.3844 0.41606 9.999e-05 ***
Residuals      129   147.396  1.1426         0.52660              
Total          136   279.901                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1


```
###### WT Mucosal SI Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2  Pr(>F)   
Sequencing_Run  1     1.383  1.3828  0.7415 0.00974 0.48715   
Sex             1     5.403  5.4029  2.8971 0.03806 0.05649 . 
Site            2    15.820  7.9100  4.2414 0.11144 0.00170 **
Residuals      64   119.356  1.8649         0.84076           
Total          68   141.961                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### RAGROR Dataset
```R
Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2     12.84   6.419   5.646 0.03561    0.0002 ***
Sex              1      0.99   0.988   0.869 0.00274    0.4283    
Site_General     1    154.48 154.481 135.881 0.42859 9.999e-05 ***
Residuals      169    192.13   1.137         0.53305              
Total          173    360.44                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   2     12.84   6.419   6.256 0.03561 9.999e-05 ***
Sex              1      0.99   0.988   0.963 0.00274    0.3952    
Site             5    177.33  35.467  34.570 0.49200 9.999e-05 ***
Residuals      165    169.28   1.026         0.46965              
Total          173    360.44                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
###### RAGROR Mucosal SI Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

               Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)   
Sequencing_Run  2    17.631  8.8156  4.7485 0.09204 0.0011 **
Sex             1     0.365  0.3651  0.1966 0.00191 0.8334   
Site            2    10.183  5.0917  2.7426 0.05316 0.0233 * 
Residuals      88   163.371  1.8565         0.85289          
Total          93   191.551                 1.00000          
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### TCR KO Dataset
[//]: Note all of these were performed in the same sequencing run
```R
Terms added sequentially (first to last)

             Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex           1     2.548   2.548   1.872 0.01877    0.1585    
Site_General  1    47.416  47.416  34.832 0.34936 9.999e-05 ***
Residuals    63    85.759   1.361         0.63187              
Total        65   135.723                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sex        1     2.548  2.5481  2.0508 0.01877    0.1359    
Site       5    59.868 11.9735  9.6367 0.44110 9.999e-05 ***
Residuals 59    73.307  1.2425         0.54012              
Total     65   135.723                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
###### TCR_KO Mucosal SI Dataset
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

          Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
Sex        1    12.636 12.6365  7.4610 0.18500 0.0003 ***
Site       2     6.552  3.2762  1.9344 0.09593 0.1046    
Residuals 29    49.116  1.6937         0.71907           
Total     32    68.305                 1.00000           
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```