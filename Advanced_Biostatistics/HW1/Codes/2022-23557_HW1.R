library(tidyverse)  # dplyr, ggplot2
library(DESeq2)
library(edgeR)


count_data <- read.csv("../Dataset/LIHC_counts.csv", row.names = 1)
count_data <- count_data[grepl("ENSG", rownames(count_data)),]
colnames(count_data) <- gsub('\\.', '-', colnames(count_data))  # Genes : 60483, Samples : 424

count_data <- count_data[rowSums(count_data) > 0, ]  # Genes : 56451, Samples : 424
pseudo_count <- log2(count_data + 1)

count_primary <- count_data[, 1]  # peek a primary tumor sample
count_recurrent <- count_data[, 194]  # peek a recurrent tumor sample
count_normal <- count_data[, 33]  # peek a normal tissue sample

pseudo_primary <- pseudo_count[, 1]
pseudo_recurrent <- pseudo_count[, 194]
pseudo_normal <- pseudo_count[, 33]

hist(count_primary, 
     freq = FALSE, breaks = 30, col = "pink", 
     xlim = c(0, 2e+06), ylim = c(0, 2e-06),
     xlab = "gene counts", main = "Histogram of Primary Tumor")

hist(count_recurrent, 
     freq = FALSE, breaks = 30, col = "pink", 
     xlim = c(0, 1e+06), ylim = c(0, 2e-05),
     xlab = "gene counts", main = "Histogram of Recurrent Tumor")

hist(count_normal,
     freq = FALSE, breaks = 500, col = "pink", 
     xlim = c(0, 1e+06), ylim = c(0, 2e-05),
     xlab = "gene counts", main = "Histogram of Normal Tissue")

hist(pseudo_primary, 
     freq = FALSE, breaks = 50, col = "pink", 
     xlim = c(0, 25), ylim = c(0, 1.2),
     xlab = "log2(gene counts)", main = "Histogram of Primary Tumor")

hist(pseudo_recurrent, 
     freq = FALSE, breaks = 50, col = "pink", 
     xlim = c(0, 25), ylim = c(0, 1.1),
     xlab = "log2(gene counts)", main = "Histogram of Recurrent Tumor")

hist(pseudo_normal, 
     freq = FALSE, breaks = 50, col = "pink", 
     xlim = c(0, 25), ylim = c(0, 1.2),
     xlab = "log2(gene counts)", main = "Histogram of Normal Tissue")


pheno_data <- read.csv("../Dataset/LIHC_phenotypes.csv", row.names = 1)
