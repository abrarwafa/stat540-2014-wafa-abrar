# Seminar 07
# Created by: Abrar Wafa
# Date: Feb.27, 2013

library(edgeR)
library(DESeq)
library(limma)
library(gplots)

dat <- read.table("../Data_RNAseq/bottomly_count_table.tsv", header = TRUE, row.names = 1)
des <- read.table("../Data_RNAseq/bottomly_phenodata.tsv", header = TRUE, row.names = 1)

############################################################################################
# edgeR
############################################################################################
# Mini exercise:
all(rownames(des) == colnames(dat))
# remove any gene that has count equal to zero across all samples
filDat <- dat[(rowSums(dat) !=0),]
# remove any gene that has count equal to zero in at least one sample in each genotype group
selrow <- apply(filDat, 1, function(row) all(row[1:10]!=0) | all(row[11:21]!=0))
fdat <- filDat [selrow, ]

group <- factor(c(rep("1", 10), rep("2", 11)))
group
dge.glm <- DGEList(counts = fdat, group = group)
str(dge.glm)
names(dge.glm)
dge.glm[["samples"]]
nrow(dge.glm[[1]])
ncol(dge.glm[[2]])

design <- model.matrix(~group)
design

dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)

dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
# plot the tagwise dispersion against log2-CPM (counts per million)
plotBCV(dge.glm.tag.disp)
fit <- glmFit(dge.glm.tag.disp, design)
colnames(coef(fit))
lrt <- glmLRT(fit, coef = 2)
topTags(lrt)
tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(tt.glm$table[tt.glm$table$FDR < 0.01, ])
interestingSamples <- rownames(tt.glm$table[tt.glm$table$FDR < 1e-50, ])
cpm(dge.glm.tag.disp)[interestingSamples, ]
summary(de.glm <- decideTestsDGE(lrt, p = 0.05, adjust = "BH"))
tags.glm <- rownames(dge.glm.tag.disp)[as.logical(de.glm)]
plotSmear(lrt, de.tags = tags.glm)
abline(h = c(-2, 2), col = "blue")

############################################################################################
# DESeq
############################################################################################
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)
deSeqDat <- estimateDispersions(deSeqDat)
# plotting the estimated dispersions against the mean normalized counts
plotDispEsts(deSeqDat)
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
str(results)
nrow(results[results$padj<1e-4,])
plotMA(results)

############################################################################################
# Voom & limma
############################################################################################
norm.factor <- calcNormFactors(dat)
dat.voomed <- voom(dat, design, plot = TRUE, lib.size = colSums(dat) * norm.factor)
vfit <- lmFit(dat.voomed, design)
vfit <- eBayes(vfit)
vhits <- topTable(vfit)
str(vhits)

############################################################################################
# Take Home Problem
############################################################################################

# edgR
dge.glm <- DGEList(counts = dat, group = group)
dge.glm[["samples"]]
design <- model.matrix(~group)
dge.glm.com.disp <- estimateGLMCommonDisp(dge.glm, design, verbose = TRUE)
dge.glm.trend.disp <- estimateGLMTrendedDisp(dge.glm.com.disp)
dge.glm.tag.disp <- estimateGLMTagwiseDisp(dge.glm.trend.disp, design)
edgfit <- glmFit(dge.glm.tag.disp, design)
lrt <- glmLRT(edgfit, coef = 2)
topTags(lrt,  adjust.method="BH")
tt.glm <- topTags(lrt, n = Inf)
class(tt.glm)
nrow(ehits <- tt.glm$table[tt.glm$table$FDR < 1e-4, ])

#deSeq
deSeqDat <- newCountDataSet(dat, group)
head(counts(deSeqDat))
deSeqDat <- estimateSizeFactors(deSeqDat)
sizeFactors(deSeqDat)
deSeqDat <- estimateDispersions(deSeqDat)
results <- nbinomTest(deSeqDat, levels(group)[1], levels(group)[2])
dhits <- (results[results$padj<1e-4,])
dhits <- na.exclude(dhits)
rownames(dhits) <- dhits$id

# Voom & limma
vhits <- topTable(vfit, coef= 2, p.value = 1e-4, number = nrow(dat))
nrow(vhits)

diagram <- list(rownames(ehits), rownames(dhits), rownames(vhits))
venn(diagram)
