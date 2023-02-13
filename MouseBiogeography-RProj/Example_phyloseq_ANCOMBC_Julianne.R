#Author: Carra Simpson
#Project: For Julianne
#Script: Import into Phyloseq / ANCOM-BC example
#Date (c): 1 Nov 2021
#Dated last ed:

library(phyloseq)
library(qiime2R)
library(ANCOMBC)

setwd("~/Documents/ETC")

# Import sample metadata and make into a phyloseq data type
SAM <- read.csv("samplemetadata.csv", row = 1)
SAM1 <- sample_data(SAM)

# Import OTU table and make into a phyloseq data type
OTU <- read.csv("otu_table.csv", row = 1)
OTU1 = otu_table(OTU, taxa_are_rows = TRUE)

# Import taxonomy and make it into a matrix for phyloseq
TAX <- read.csv("taxonomy.csv", row =1)
TAX1 <- as.matrix(TAX)

# Make taxonomy into a phyloseq data type
TAX2 <- tax_table(TAX1)

# Make the final phyloseq object with the factor sample data
DADA2_biom2 <- merge_phyloseq(OTU1, SAM1, TAX2)


########################################################################

# Collapse phyloseq object at specific level of the taxonomy
DADA2_biom2species <- tax_glom(DADA2_biom2, taxrank = "Species")

# Run ANCOMBC (have to collapse phyloseq object to specific groups if >3)
# Species level ANCOM-BC (cannot currently do random effects in ANCOM-BC)
# Global = TRUE gives omnibus results, might need to collapse into groups
out = ancombc(phyloseq = DADA2_biom2species, formula = "Group + Age + ETC + ETC", 
              p_adj_method = "fdr", zero_cut = 0.90, lib_cut = 0, 
              group = "Group", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.15, global = TRUE)

res = out$res
tab_diff = res$diff_abn

res$q_val
res$p_val
res$beta
res$W

write.csv(res, "FINALRESULTS.csv")
