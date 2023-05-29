# Omit two sequencing runs and luminal samples and perform combat seq
Then filter on sequencing depth of 10k or less, and prevalence filter to 10%
### Genotype* Site_General interaction 
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                       Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run          1     21.98  21.975  18.305 0.03730 9.999e-05 ***
Sex                     1      2.78   2.782   2.317 0.00472   0.09689 .  
Genotype                2     38.68  19.340  16.110 0.06566 9.999e-05 ***
Site_General            1    192.28 192.285 160.170 0.32640 9.999e-05 ***
Genotype:Site_General   2      2.06   1.028   0.856 0.00349   0.49525    
Residuals             276    331.34   1.201         0.56243              
Total                 283    589.11                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```
### Genotype* Site Interaction
```R
Permutation: free
Number of permutations: 10000

Terms added sequentially (first to last)

                Df SumsOfSqs MeanSqs F.Model      R2    Pr(>F)    
Sequencing_Run   1     21.98  21.975  21.635 0.03730 9.999e-05 ***
Sex              1      2.78   2.782   2.739 0.00472   0.06259 .  
Genotype         2     38.68  19.340  19.040 0.06566 9.999e-05 ***
Site             5    249.42  49.884  49.113 0.42339 9.999e-05 ***
Genotype:Site   10      8.11   0.811   0.798 0.01376   0.72883    
Residuals      264    268.15   1.016         0.45517              
Total          283    589.11                 1.00000              
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
```