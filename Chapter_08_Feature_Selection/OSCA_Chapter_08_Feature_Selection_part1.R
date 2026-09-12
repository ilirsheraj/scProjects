# Part 2: Basics - Feature Selection Dataset 1

library(DropletUtils)
# Data downloaded directly from terminal (If you want manually download it, but
# OSCA's current package does not download it at all)
# wget wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
# tar -xvf pbmc4k_raw_gene_bc_matrices.tar.gz
sce <- read10xCounts("raw_gene_bc_matrices/GRCh38", col.names = TRUE)
sce
dim(sce)

# About 750k cells (mostly empty droplets of course)
totals <- colSums(counts(sce))
summary(totals)
table(totals == 0)
head(rowData(sce))

# Gene annotation
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

library(EnsDb.Hsapiens.v86)
# Get chromosome locations
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")
head(location)

# Cell detection
# This is droplet-based 10x Genomics, so we detect empty drops and remove them
set.seed(100)
empty_drops <- emptyDrops(counts(sce))
table(empty_drops$FDR <= 0.001, useNA = "ifany")

keep <- !is.na(empty_drops$FDR) & empty_drops$FDR <= 0.001

# Plot barcode ranks: Right upper part is what you want
bar_out <- barcodeRanks(counts(sce))

plot(bar_out$rank, bar_out$total, log = "xy", 
     xlab = "Rank",
     ylab = "Total UMI count")

abline(h = metadata(bar_out)$knee, lty = 2)
abline(h = metadata(bar_out)$inflection, lty = 2)

plot(bar_out$rank, bar_out$total, log = "xy", 
     xlab = "Rank",
     ylab = "Total UMI count")

is.cell <- !is.na(empty_drops$FDR) & empty_drops$FDR <= 0.001

points(bar_out$rank[is.cell],
       bar_out$total[is.cell],
       pch = 16,
       cex = 0.5,
       col="red")

# Filter out empty droplets
sce <- sce[, keep]

# Now we have a data structure ready for QC
sce
# 3302 cells passed the filter

# ------------------------------------------------------------------------------
# Quality Control (QC)
# ------------------------------------------------------------------------------
# Quality Control: Old Method
# stats <- perCellQCMetrics(sce, subsets=list(Mito=which(location=="MT")))
# high.mito <- isOutlier(stats$subsets_Mito_percent, type="higher")
# sce.pbmc <- sce[,!high.mito]
# sce.pbmc

# More modern way: Quick and shorter, In this case more stringent
is_mito <- grepl("^MT", location)
sum(is_mito)
sce <- scrapper::quickRnaQc.se(
  sce,
  subsets = list(Mito = is_mito))

table(sce$keep)

# Check in more details whhat's going on
qc.metrics <- scrapper::computeRnaQcMetrics(
  counts(sce),
  subsets = list(Mito = is_mito))

qc.thresholds <- scrapper::suggestRnaQcThresholds(qc.metrics)

# The cut-offs
qc.thresholds

# Remove low-quality cells
sce <- sce[, sce$keep]
sce
# 3937 passed QC filters

# ------------------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------------------
library(scran)
set.seed(1000)

# First clustering
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, cluster=clusters)
summary(sizeFactors(sce))

# Normalization Step
sce <- scrapper::normalizeRnaCounts.se(
  sce, 
  size.factors = sizeFactors(sce))

# ------------------------------------------------------------------------------
# Quantifying per-gene variation
# ------------------------------------------------------------------------------
dec <- scrapper::modelGeneVariances(logcounts(sce))
head(dec)
dec_formatted <- scrapper::formatModelGeneVariancesResult(dec)
head(dec_formatted)

# Instead of bio in scran, we have residuals
dec_formatted[order(dec_formatted$residuals, decreasing=TRUE),]

ord <- order(dec_formatted$means)

plot(dec_formatted$means,
     dec_formatted$variances,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression",
     pch = 16)

lines(dec_formatted$means[ord],
      dec_formatted$fitted[ord],
      col = "dodgerblue",
      lwd = 2)

# # Older Method
# dec.pbmc <- modelGeneVar(sce)
# dec.pbmc[order(dec.pbmc$bio, decreasing=TRUE),]
# fit.pbmc <- metadata(dec.pbmc)
# plot(fit.pbmc$mean, fit.pbmc$var, xlab="Mean of log-expression",
#      ylab="Variance of log-expression")
# curve(fit.pbmc$trend(x), col="dodgerblue", add=TRUE, lwd=2)

# ------------------------------------------------------------------------------
# Quantifying technical noise
# ------------------------------------------------------------------------------
set.seed(0010101)
dec_poisson <- modelGeneVarByPoisson(sce)
# Get a warning, but ignore it

head(dec_poisson)

dec_poisson <- dec_poisson[order(dec_poisson$bio, decreasing=TRUE),]
head(dec_poisson)

plot(dec_poisson$mean, dec_poisson$total, 
     pch=16, 
     xlab="Mean of log-expression",
     ylab="Variance of log-expression")
curve(metadata(dec_poisson)$trend(x), col="dodgerblue", add=TRUE)

# ------------------------------------------------------------------------------
# Selecting highly variable genes
# ------------------------------------------------------------------------------
n_genes <- 0.1*dim(sce)[1]
hvg <- scrapper::chooseHighlyVariableGenes(
  dec_formatted$residuals,
  top = n_genes)

length(hvg)
str(hvg)

# Keep only Highly Variable Genes
sce_hvg <- sce[hvg,]
dim(sce_hvg)

# Keep the entire dataset, but compute dimensionality reduction on HVG
library(scater)
sce <- runPCA(sce, subset_row=hvg)
reducedDimNames(sce)

# Add more options
rowSubset(sce) <- hvg

# Bigger subset of genes
n_genes_2 <- 0.2*dim(sce)[1]
rowSubset(sce, "HVGs.more") <- scrapper::chooseHighlyVariableGenes(
  dec_formatted$residuals,
  top = n_genes_2)

# Check their availability
colnames(rowData(sce))

# Recycling the class above.
altExp(sce_hvg, "original") <- sce
altExpNames(sce_hvg)
dim(sce_hvg)

# No need for explicit subset_row= specification in downstream operations.
sce_hvg <- runPCA(sce_hvg)

# Recover original data:
sce_hgv_original <- altExp(sce_hvg, "original", withColData=TRUE)
sce_hgv_original

# EOF