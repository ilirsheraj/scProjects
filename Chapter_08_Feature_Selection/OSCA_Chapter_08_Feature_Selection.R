library(DropletUtils)

# Data downloaded directly from terminal (If you want manually download it, but
# OSCA's current package does not download it at all)
# wget wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
# tar -xvf pbmc4k_raw_gene_bc_matrices.tar.gz
sce <- read10xCounts("raw_gene_bc_matrices/GRCh38", col.names = TRUE)
sce
dim(sce)

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
     ylab = "Variance of log-expression")

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
# Another Dataset
# ------------------------------------------------------------------------------
library(scRNAseq)
sce.416b <- LunSpikeInData(which="416b") 

