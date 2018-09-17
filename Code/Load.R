#!/usr/bin/env Rscript

d <- read.csv(file = "./raw_data//gene_count_matrix.csv")

saveRDS(d, "./rds/Freichel.rds")
