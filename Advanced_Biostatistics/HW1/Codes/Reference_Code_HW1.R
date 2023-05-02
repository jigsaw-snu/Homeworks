#Necessary Packages#
install.packages("BiocManager")
BiocManager::install("edgeR")
BiocManager::install("DESeq2")

library(edgeR)
library(DESeq2)


##Input and data preparations [count data]##
BiocManager::install("pasilla")
library(pasilla)

counts <- system.file("extdata", "pasilla_gene_counts.tsv", package = "pasilla", mustWork = TRUE)
dat <- as.matrix(read.csv(counts, sep = "\t", row.names = "gene_id"))
head(dat)
dim(dat) # 14599 genes 7 samples

##Input and data preparations [sample information]##
pasAnno <- system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla", mustWork = TRUE)
coldata <- read.csv(pasAnno, row.names = 1)
head(coldata)
coldata <- coldata[, c("condition", "type")]
head(coldata)

## Basic exploratory analysis ##
#Gene filtering
raw_counts_wn <- dat[ rowSums(dat)>0,] # remove all zero counts
dim(raw_counts_wn)  
pseudo_counts <- log2(raw_counts_wn + 1)

colnames(raw_counts_wn)
raw_counts_wn1<-raw_counts_wn[, ]
hist(raw_counts_wn1, probability = T, breaks=30, col="blue", xlim = c(0, 200000), 
     ylim = c(0, 0.00012), xlab = "counts", main = "Histogram of treated1") 
pseudo_counts1<-pseudo_counts[,5]
hist(pseudo_counts1, probability = T, breaks=30, col="blue", 
     ylim = c(0, 0.20), xlab = "log2(counts+1)", main = "Histogram of treated1") 
     
                      
library(MASS) 
library(reshape2) 
library(reshape)

df_raw <- melt(pseudo_counts, id=rownames(raw_counts_wn)) 
names(df_raw)[1:2] <- c("id", "sample")  
df_raw$method <- rep("Raw counts", nrow(df_raw))
head(df_raw)

#Relationship between mean and variance#
df <- data.frame(mean = apply(raw_counts_wn[, coldata$condition == "untreated"], 1 , mean), var= apply(raw_counts_wn[, coldata$condition == "untreated"], 1 , var))
df<- df[df$mean<=5000,]
head(df)
plot(df$mean, df$var, xlab="Mean", ylab="Variance", col="green4", main="Plot of variance versus mean in counts")

########################Total count normalization ###################################
library(edgeR)
dge2 <- DGEList(raw_counts_wn)
dge2

# count per million read (normalized count)#
norm_counts <- cpm(dge2)
head(norm_counts)

pseudo_TC <- log2(cpm(dge2) + 1)
df_TC <- melt(pseudo_TC, id=rownames(raw_counts_wn))
names(df_TC)[1:2] <- c("id", "sample") 
df_TC$method <- rep("TC", nrow(df_TC))
head(df_TC)

##############Upper quartile normalization is obtained with the function calcNormFactors#####################
# normalize for library size by calculating scaling factor#
dge2 <- calcNormFactors(dge2, method="upperquartile")
# normalization factors for each library#
dge2$samples

# count per million read (normalized count)#
norm_counts_uq <- cpm(dge2)
head(norm_counts_uq)

#Pseudo counts are obtained and stored in a data frame#
pseudo_UQ <- log2(cpm(dge2) + 1)
head(pseudo_UQ)
df_UQ<- melt(pseudo_UQ, id=rownames(raw_counts_wn))
names(df_UQ)[1:2] <- c("id", "sample") 
df_UQ$method <- rep("UQ", nrow(df_UQ))
head(df_UQ)


###############################Trimmed mean (TMM)####################################
#TMM normalization works similarly as UQ and updates the value of the variable norm.factors#
dge2 <- calcNormFactors(dge2, method="TMM")
dge2$samples

# count per million read (normalized count)#
norm_counts_tmm <- cpm(dge2)
head(norm_counts_tmm)

#Pseudo counts are obtained with the function cpm and stored in a data frame#
pseudo_TMM <- log2(cpm(dge2) + 1)
df_TMM <- melt(pseudo_TMM, id=rownames(raw_counts_wn))
names(df_TMM)[1:2] <- c("id", "sample") 
df_TMM$method <- rep("TMM", nrow(df_TMM))
head(df_TMM)


###########Relative log expectation (RLE)#####################
library(DESeq2)
meta = data.frame(Group=c(rep("untreated", 4), rep("treated", 3)))
rownames(meta)<-colnames(raw_counts_wn)
dds <- DESeqDataSetFromMatrix(countData = raw_counts_wn, colData = meta, design = ~ Group)
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
 #DESeq2 normalization counts#
deseq_normcount<- counts(dds, normalized=T)
pseudo_deseq<-log2(deseq_normcount+1)
df_deseq <- melt(pseudo_deseq, id=rownames(raw_counts_wn))
names(df_deseq)[1:2] <- c("id", "sample") 
df_deseq$method <- rep("DESeq2(RLE)", nrow(df_raw))
head(df_deseq)

##############################Comparison Methods################################################

df_allnorm <- rbind(df_raw, df_deseq, df_TC, df_UQ, df_TMM)

df_allnorm$method <- factor(df_allnorm$method, levels = c("Raw counts", "DESeq2(RLE)", "TC", "UQ", "TMM"))


library(ggplot2)
boxplots <- ggplot(df_allnorm, aes(x=sample, y = value, fill=method)) + 
  geom_boxplot() +
  facet_wrap(~method, scale="free")+
  theme(axis.text.x=element_blank() 
  )
boxplots


####################################Differential Analysis#####################################

Anno <- system.file("extdata", "pasilla_sample_annotation.csv", package = "pasilla", mustWork = TRUE)
sampleAnno <- read.csv(Anno)
sname <- sub("fb$", "", sampleAnno$file)
rownames(sampleAnno) <- sname
sampleAnno <- sampleAnno[, -1]
head(sampleAnno)

sampleAnno<- sampleAnno[match(colnames(dat), rownames(sampleAnno)),]
head(sampleAnno)
all.equal(colnames(dat), rownames(sampleAnno))

sampleAnno <- sampleAnno[, c("condition", "type")]
sampleAnno$condition <- factor(sampleAnno[, "condition"], levels = c("untreated", "treated")) 
sampleAnno$type <- factor(sampleAnno[, "type"], levels = c("single-read", "paired-end"))
table(sampleAnno)

#######################################Testing with DESeq2###################################################

library(DESeq2)
  #Wald Test#
dds <- DESeqDataSetFromMatrix(countData = dat, colData = sampleAnno, design = ~ condition)
dds <- DESeq(dds, test = "Wald", quiet = FALSE)

res <- results(dds, alpha = 0.05) 
head(res)
summary(res)
 
 #LRT#

dds1 <- DESeq(dds, test = "LRT", reduced = ~ 1)
res1 <- results(dds1) 
head(res1)
summary(res1)

 #Remove NA values and update the results #

table(is.na(res$padj))
res_update <- res[complete.cases(res),]
summary(res_update)

res_update_ordered <- res_update[order(res_update$padj), ]
head(res_update_ordered)

  #Save results as a fie#

save1 <- merge.data.frame(as.data.frame(res_update_ordered), as.data.frame(counts(dds, normalized = T)), by="row.names", 
                          sort = FALSE)
colnames(save1)[1]  <- "gene"
write.csv(save1, "treated_vs_untreated.csv")

  #Principal Component Analysis#

dds_rlog <- rlogTransformation(dds)
plotPCA(dds_rlog, intgroup = c("condition", "type"))


###########################################################Testing with edgeR################################################
library(edgeR)
condition <- sampleAnno$condition
dat_1 <- DGEList(counts =  dat, group = condition)
dat_1
dat_norm <- calcNormFactors(dat_1, method="TMM")
dat_norm_d1 <- estimateCommonDisp(dat_norm)

dat_norm_d1 <- estimateTagwiseDisp(dat_norm_d1)
dat_norm_d1
plotBCV(dat_norm_d1, )
head(dat_norm_d1$counts[order(dat_norm_d1$tagwise.dispersion),], 3)
tail(dat_norm_d1$counts[order(dat_norm_d1$tagwise.dispersion),], 3)

