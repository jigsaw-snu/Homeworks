library(tidyverse)  # dplyr, ggplot2
library(DESeq2)
library(edgeR)


count_data <- read.csv("../Dataset/LIHC_counts.csv", row.names = 1)
count_data <- count_data[grepl("ENSG", rownames(count_data)),]
colnames(count_data) <- gsub('\\.', '-', colnames(count_data))  # Genes : 60483, Samples : 424

count_data <- count_data[rowSums(count_data) > 0, ]  # Genes : 56451, Samples : 424
pseudo_count <- log2(count_data + 1)



pheno_data <- read.csv("../Dataset/LIHC_phenotypes.csv", row.names = 1)
