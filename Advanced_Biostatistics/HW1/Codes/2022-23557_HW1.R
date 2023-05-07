library(tidyverse)  # dplyr, ggplot2
library(DESeq2)
library(edgeR)
library(reshape)
library(reshape2)


pheno_data <- read.csv("../Dataset/LIHC_phenotypes.csv", row.names = 1)
colnames(pheno_data) <- c("Sample", "Condition")

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
     freq = FALSE, breaks = 200, col = "pink", 
     xlim = c(0, 1e+06), ylim = c(0, 1e-06),
     xlab = "gene counts", main = "Histogram of Primary Tumor")

hist(count_recurrent, 
     freq = FALSE, breaks = 200, col = "pink", 
     xlim = c(0, 1e+06), ylim = c(0, 1e-06),
     xlab = "gene counts", main = "Histogram of Recurrent Tumor")

hist(count_normal,
     freq = FALSE, breaks = 200, col = "pink", 
     xlim = c(0, 1e+06), ylim = c(0, 1e-06),
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


primary_mean_var <- data.frame(mean = apply(count_data[, pheno_data$Condition == "Primary Tumor"], 1, mean),
                               var = apply(count_data[, pheno_data$Condition == "Primary Tumor"], 1, var))

primary_mean_var2 <- primary_mean_var[primary_mean_var$mean < 1e+05, ]

recurrent_mean_var <- data.frame(mean = apply(count_data[, pheno_data$Condition == "Recurrent Tumor"], 1, mean),
                                 var = apply(count_data[, pheno_data$Condition == "Recurrent Tumor"], 1, var))

recurrent_mean_var2 <- recurrent_mean_var[recurrent_mean_var$mean < 1e+05, ]

normal_mean_var <- data.frame(mean = apply(count_data[, pheno_data$Condition == "Solid Tissue Normal"], 1, mean),
                      var = apply(count_data[, pheno_data$Condition == "Solid Tissue Normal"], 1, var))

normal_mean_var2 <- normal_mean_var[normal_mean_var$mean < 1e+05, ]


plot(primary_mean_var2$mean, primary_mean_var2$var,
     xlab = "Mean", ylab = "Variance",
     col = "orange2", main = "Mean-Variance plot")

plot(recurrent_mean_var2$mean, recurrent_mean_var2$var,
     xlab = "Mean", ylab = "Variance",
     col = "orange2", main = "Mean-Variance plot")

plot(normal_mean_var2$mean, normal_mean_var2$var,
     xlab = "Mean", ylab = "Variance",
     col = "orange2", main = "Mean-Variance plot")


long_RAW <- reshape2::melt(as.matrix(pseudo_count), id = rownames(count_data))
names(long_RAW)[1:2] <- c("gene", "sample")
long_RAW$method <- rep("RAW", nrow(long_RAW))


###### Total Count Normalization
DGE <- edgeR::DGEList(count_data)
norm_TC <- edgeR::cpm(DGE)
pseudo_TC <-log2(norm_TC + 1) 
long_TC <- reshape2::melt(pseudo_TC, id = rownames(count_data))
colnames(long_TC)[1:2] <- c("gene", "sample")
long_TC$method <- rep("TC", nrow(long_TC))


###### Upper Quartile Normalization
DGE <- edgeR::calcNormFactors(DGE, method = "upperquartile")
norm_UQ <- edgeR::cpm(DGE)
pseudo_UQ <- log2(norm_UQ + 1)
long_UQ <- reshape2::melt(pseudo_UQ, id = rownames(count_data))
colnames(long_UQ)[1:2] <- c("gene", "sample")
long_UQ$method <- rep("UQ", nrow(long_UQ))


###### Trimmed Mean of M-values Normalization
DGE <- edgeR::calcNormFactors(DGE, method = "TMM")
norm_TMM <- edgeR::cpm(DGE)
pseudo_TMM <- log2(norm_TMM + 1)
long_TMM <- reshape2::melt(pseudo_TMM, id = rownames(count_data))
colnames(long_TMM)[1:2] <- c("gene", "sample")
long_TMM$method <- rep("TMM", nrow(long_TMM))


###### Relative Log Expression Normalization
meta_data <- pheno_data[match(colnames(count_data), pheno_data$Sample), ]
rownames(meta_data) <- NULL
meta_data <- meta_data %>% tibble::column_to_rownames(var = "Sample")
meta_data$Condition <- as.factor(meta_data$Condition)
levels(meta_data$Condition) <- c("Solid Normal Tissue", "Recurrent Tumor", "Primary Tumor")

dds <- DESeq2::DESeqDataSetFromMatrix(countData = round(count_data),
                                      colData = meta_data,
                                      design = ~Condition)

dds <- DESeq2::estimateSizeFactors(dds)

norm_RLE <- DESeq2::counts(dds, normalized = TRUE)
pseudo_RLE <- log2(norm_RLE + 1)
long_RLE <- reshape2::melt(pseudo_RLE, id = rownames(count_data))
names(long_RLE)[1:2] <- c("gene", "sample")
long_RLE$method <- rep("RLE", nrow(long_RLE))



ggplot(long_RAW, aes(x = sample, y = value, fill = method)) +
    geom_boxplot() +
    facet_wrap(~method, scales = "free")

all_norms <- rbind(long_RAW, long_TC, long_UQ, long_TMM, long_RLE)






