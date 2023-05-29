```R
testing <- site_heatmap$feature
testing2<-unique(data$feature)
c(setdiff(testing, testing2), setdiff(testing2, testing))
```