library(renv)
devtools::install_github("biobakery/melonnpan")
renv::restore()
renv::snapshot()
library(melonnpan)
library(here)

here::i_am("MouseBiogeography-RProj/MelonnPan.R")

###Running of MelonnPan
df<-read.tsv(here("Shotgun/relab_normalized/"))
predict_metabolites()
