Prior to doing any sort of batch correction: Genotype* Site is significant, and Sequencing_Run accounts for 29% of the variation
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