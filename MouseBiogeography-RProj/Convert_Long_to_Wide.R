library(reshape2)

setwd("C:/Users/Jacobs Laboratory/Desktop/Mouse_Biogeography_Julianne/Regional-Mouse-Biogeography-Analysis/2021-8-Pathway-Batch-Correction/")

data_long= read.csv ("pathwaynames_long.csv", header=TRUE)

data_wide <- dcast(data_long, L1...Terminal.Pathway ~ Splay, value.var="L2")
data_wide

write.csv(data_wide, "pathwaynames_wide_L2.csv")

data_wide2 <- dcast(data_long, L1...Terminal.Pathway ~ Splay, value.var="L3")
data_wide2

write.csv(data_wide2, "pathwaynames_wide_L3.csv")

data_wide3 <- dcast(data_long, L1...Terminal.Pathway ~ Splay, value.var="L4")
data_wide3

write.csv(data_wide3, "pathwaynames_wide_L4.csv")

library(dplyr)
library(tidyr)

data_wide1 <- data_long%>%
  pivot_wider(names_from = Splay,
              values_from = L2)

#find duplicates
data_long %>%
  dplyr::group_by(L1...Terminal.Pathway, L2, L4, Splay) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 