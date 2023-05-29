### Interaction effect between Genotype* Site after controlling for covariates: Sequencing_Run, Sex 
Shannon 
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Genotype * Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                        3.173195 0.2482909 216 12.780151  0.0000
Sequencing_Run2014-September       0.023771 0.2006599  48  0.118467  0.9062
SexMale                            0.107806 0.1361990  48  0.791536  0.4325
GenotypeTCR_KO                     0.379745 0.3376719  48  1.124597  0.2664
GenotypeRAGROR                    -2.065236 0.3294287  48 -6.269147  0.0000
SiteJejunum                       -1.182253 0.2525130 216 -4.681950  0.0000
SiteIleum                         -1.631777 0.2506998 216 -6.508887  0.0000
SiteCecum                          1.016420 0.2517010 216  4.038206  0.0001
SiteProximal_Colon                 1.175978 0.2605406 216  4.513609  0.0000
SiteDistal_Colon                   0.606822 0.2626541 216  2.310349  0.0218
GenotypeTCR_KO:SiteJejunum        -0.197976 0.4258552 216 -0.464891  0.6425
GenotypeRAGROR:SiteJejunum         0.559088 0.3802707 216  1.470237  0.1430
GenotypeTCR_KO:SiteIleum          -1.193420 0.4247826 216 -2.809483  0.0054
GenotypeRAGROR:SiteIleum           0.811556 0.3829384 216  2.119286  0.0352
GenotypeTCR_KO:SiteCecum          -0.368368 0.4253742 216 -0.865986  0.3875
GenotypeRAGROR:SiteCecum           2.343749 0.4671022 216  5.017637  0.0000
GenotypeTCR_KO:SiteProximal_Colon  0.126494 0.4306637 216  0.293719  0.7693
GenotypeRAGROR:SiteProximal_Colon  2.509420 0.3895312 216  6.442153  0.0000
GenotypeTCR_KO:SiteDistal_Colon    0.514229 0.4319456 216  1.190496  0.2352
GenotypeRAGROR:SiteDistal_Colon    2.498302 0.3908554 216  6.391883  0.0000
```
Observed_OTUs
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Genotype * Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                        64.24217  6.335522 216 10.139995  0.0000
Sequencing_Run2014-September        8.79816  4.943090  48  1.779890  0.0814
SexMale                             5.57011  3.351505  48  1.661974  0.1030
GenotypeTCR_KO                     10.46689  8.684448  48  1.205245  0.2340
GenotypeRAGROR                    -38.66099  8.437004  48 -4.582312  0.0000
SiteJejunum                       -27.88389  6.614623 216 -4.215492  0.0000
SiteIleum                         -22.92188  6.566259 216 -3.490858  0.0006
SiteCecum                          62.06113  6.589287 216  9.418490  0.0000
SiteProximal_Colon                 66.98063  6.816041 216  9.826911  0.0000
SiteDistal_Colon                   43.73501  6.872660 216  6.363622  0.0000
GenotypeTCR_KO:SiteJejunum          0.52026 11.162492 216  0.046608  0.9629
GenotypeRAGROR:SiteJejunum         15.69639  9.966757 216  1.574875  0.1167
GenotypeTCR_KO:SiteIleum          -27.80540 11.133901 216 -2.497363  0.0133
GenotypeRAGROR:SiteIleum           13.76896 10.034924 216  1.372104  0.1715
GenotypeTCR_KO:SiteCecum          -10.51568 11.147497 216 -0.943322  0.3466
GenotypeRAGROR:SiteCecum           44.20475 12.215155 216  3.618845  0.0004
GenotypeTCR_KO:SiteProximal_Colon -11.79881 11.283014 216 -1.045714  0.2969
GenotypeRAGROR:SiteProximal_Colon  41.78418 10.202084 216  4.095651  0.0001
GenotypeTCR_KO:SiteDistal_Colon     2.26499 11.317308 216  0.200135  0.8416
GenotypeRAGROR:SiteDistal_Colon    41.56025 10.237639 216  4.059554  0.0001
```
Chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Genotype * Site 
                                      Value Std.Error  DF   t-value p-value
(Intercept)                        66.42565  7.042480 216  9.432139  0.0000
Sequencing_Run2014-September       10.73674  5.535572  48  1.939590  0.0583
SexMale                             5.93114  3.754062  48  1.579925  0.1207
GenotypeTCR_KO                     11.68195  9.637641  48  1.212118  0.2314
GenotypeRAGROR                    -31.12019  9.371157  48 -3.320848  0.0017
SiteJejunum                       -25.29293  7.313786 216 -3.458253  0.0007
SiteIleum                         -17.16840  7.260517 216 -2.364625  0.0189
SiteCecum                          72.36969  7.286720 216  9.931723  0.0000
SiteProximal_Colon                 73.35201  7.538577 216  9.730220  0.0000
SiteDistal_Colon                   48.95322  7.600892 216  6.440457  0.0000
GenotypeTCR_KO:SiteJejunum         -2.10358 12.340654 216 -0.170460  0.8648
GenotypeRAGROR:SiteJejunum          5.31175 11.018923 216  0.482057  0.6303
GenotypeTCR_KO:SiteIleum          -31.22819 12.309158 216 -2.536989  0.0119
GenotypeRAGROR:SiteIleum            2.99641 11.094705 216  0.270075  0.7874
GenotypeTCR_KO:SiteCecum          -11.41734 12.324633 216 -0.926383  0.3553
GenotypeRAGROR:SiteCecum           31.60235 13.511227 216  2.338969  0.0203
GenotypeTCR_KO:SiteProximal_Colon -13.94907 12.475192 216 -1.118144  0.2647
GenotypeRAGROR:SiteProximal_Colon  32.05632 11.280831 216  2.841663  0.0049
GenotypeTCR_KO:SiteDistal_Colon    -1.09810 12.512947 216 -0.087757  0.9302
GenotypeRAGROR:SiteDistal_Colon    32.48452 11.319947 216  2.869671  0.0045
```
Pielou's evenness
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Genotype * Site 
                                       Value  Std.Error  DF   t-value p-value
(Intercept)                        0.5329087 0.03796133 216 14.038201  0.0000
Sequencing_Run2014-September      -0.0124179 0.03106258  48 -0.399771  0.6911
SexMale                            0.0065394 0.02109223  48  0.310038  0.7579
GenotypeTCR_KO                     0.0419463 0.05148068  48  0.814796  0.4192
GenotypeRAGROR                    -0.3046166 0.05030078  48 -6.055903  0.0000
SiteJejunum                       -0.1479411 0.03822468 216 -3.870303  0.0001
SiteIleum                         -0.2465776 0.03795193 216 -6.497104  0.0000
SiteCecum                          0.0761761 0.03811004 216  1.998846  0.0469
SiteProximal_Colon                 0.0883601 0.03945770 216  2.239363  0.0262
SiteDistal_Colon                   0.0304258 0.03977502 216  0.764947  0.4451
GenotypeTCR_KO:SiteJejunum        -0.0430670 0.06445084 216 -0.668214  0.5047
GenotypeRAGROR:SiteJejunum         0.0456904 0.05755359 216  0.793877  0.4281
GenotypeTCR_KO:SiteIleum          -0.1877803 0.06428945 216 -2.920858  0.0039
GenotypeRAGROR:SiteIleum           0.0923323 0.05796077 216  1.593013  0.1126
GenotypeTCR_KO:SiteCecum          -0.0395750 0.06438291 216 -0.614682  0.5394
GenotypeRAGROR:SiteCecum           0.3307751 0.07074902 216  4.675331  0.0000
GenotypeTCR_KO:SiteProximal_Colon  0.0363271 0.06518968 216  0.557253  0.5779
GenotypeRAGROR:SiteProximal_Colon  0.3633178 0.05896996 216  6.161065  0.0000
GenotypeTCR_KO:SiteDistal_Colon    0.0785800 0.06538224 216  1.201856  0.2307
GenotypeRAGROR:SiteDistal_Colon    0.3658435 0.05916858 216  6.183070  0.0000
 
```
### Interaction between Genotype* Site_General after accounting for covariates: Sequencing_Run, Sex
shannon
```R
Fixed effects: shannon ~ Sequencing_Run + Sex + Genotype * Site_General 
                                  Value Std.Error  DF    t-value p-value
(Intercept)                    4.089201 0.1800606 228  22.710138  0.0000
Sequencing_Run2014-September   0.081259 0.2039879  48   0.398353  0.6921
SexMale                        0.128992 0.1379464  48   0.935089  0.3544
GenotypeTCR_KO                 0.418554 0.2425075  48   1.725942  0.0908
GenotypeRAGROR                 0.397051 0.2523671  48   1.573309  0.1222
Site_GeneralSI                -2.000282 0.1729204 228 -11.567641  0.0000
GenotypeTCR_KO:Site_GeneralSI -0.425386 0.2995713 228  -1.419983  0.1570
GenotypeRAGROR:Site_GeneralSI -1.866412 0.2804023 228  -6.656193  0.0000
```
otus
```R
Fixed effects: observed_otus ~ Sequencing_Run + Sex + Genotype * Site_General 
                                  Value Std.Error  DF    t-value p-value
(Intercept)                   121.75204  4.405083 228  27.638988  0.0000
Sequencing_Run2014-September    9.45868  5.006268  48   1.889367  0.0649
SexMale                         6.05025  3.386529  48   1.786565  0.0803
GenotypeTCR_KO                  2.94369  5.914669  48   0.497693  0.6210
GenotypeRAGROR                  2.05630  6.157171  48   0.333968  0.7399
Site_GeneralSI                -76.80761  4.172633 228 -18.407467  0.0000
GenotypeTCR_KO:Site_GeneralSI  -0.13179  7.225380 228  -0.018239  0.9855
GenotypeRAGROR:Site_GeneralSI -28.78445  6.764337 228  -4.255324  0.0000
```
chao1
```R
Fixed effects: chao1 ~ Sequencing_Run + Sex + Genotype * Site_General 
                                  Value Std.Error  DF    t-value p-value
(Intercept)                   131.39606  4.901097 228  26.809520  0.0000
Sequencing_Run2014-September   11.26165  5.600473  48   2.010839  0.0500
SexMale                         6.49745  3.790620  48   1.714087  0.0930
GenotypeTCR_KO                  1.94789  6.545318  48   0.297600  0.7673
GenotypeRAGROR                 -1.29593  6.817642  48  -0.190085  0.8500
Site_GeneralSI                -81.13026  4.527888 228 -17.917904  0.0000
GenotypeTCR_KO:Site_GeneralSI  -0.20425  7.834093 228  -0.026071  0.9792
GenotypeRAGROR:Site_GeneralSI -25.34660  7.336637 228  -3.454798  0.0007
```
pielou
```R
Fixed effects: pielou_e ~ Sequencing_Run + Sex + Genotype * Site_General 
                                   Value  Std.Error  DF   t-value p-value
(Intercept)                    0.5952772 0.02768024 228 21.505494  0.0000
Sequencing_Run2014-September  -0.0035729 0.03146869  48 -0.113538  0.9101
SexMale                        0.0095988 0.02128797  48  0.450903  0.6541
GenotypeTCR_KO                 0.0591620 0.03715361  48  1.592362  0.1179
GenotypeRAGROR                 0.0527165 0.03867831  48  1.362948  0.1793
Site_GeneralSI                -0.2152252 0.02617976 228 -8.221055  0.0000
GenotypeTCR_KO:Site_GeneralSI -0.0833282 0.04533085 228 -1.838222  0.0673
GenotypeRAGROR:Site_GeneralSI -0.2906784 0.04243923 228 -6.849285  0.0000
 
```