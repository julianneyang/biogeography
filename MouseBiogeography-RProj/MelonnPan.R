library(renv)
devtools::install_github("biobakery/melonnpan")
renv::restore()
renv::snapshot()
library(melonnpan)
library(here)
library(stringr)

here::i_am("MouseBiogeography-RProj/MelonnPan.R")

###Running of MelonnPan
df<-read.delim(here("Shotgun/relab_normalized/merged_humann_genefamilies_relab_unstratified.tsv"),row.names=1)
shotgun_dat <- as.data.frame(t(df))
shotgun_dat <- shotgun_dat[,-1]

# Mean relative abundance of 0.01% and prevalent in 10% of samples
threshold <- 0.01/100
sample_fraction <- 0.10
binary_abundance <- shotgun_dat > threshold
proportion_samples <- colMeans(binary_abundance)
filtered_data <- shotgun_dat[,proportion_samples>=sample_fraction]

melonnpan::melonnpan.predict(metag = filtered_data,
                            output = here("Shotgun/melonnpan"))

# Prevalent in 15 % of samples 
threshold <- 0.15
binary_presence <- shotgun_dat > 0
proportion_present <- colMeans(binary_presence)
filtered_dat2 <- shotgun_dat[,proportion_present>=threshold]

melonnpan::melonnpan.predict(metag = filtered_dat2,
                             output = here("Shotgun/melonnpan/"))

