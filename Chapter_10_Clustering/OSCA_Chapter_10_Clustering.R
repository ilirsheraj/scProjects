library(DropletUtils)
# Data downloaded directly from terminal (If you want manually download it, but
# OSCA's current package does not download it at all)
# wget wget https://cf.10xgenomics.com/samples/cell-exp/2.1.0/pbmc4k/pbmc4k_raw_gene_bc_matrices.tar.gz
# tar -xvf pbmc4k_raw_gene_bc_matrices.tar.gz
# ------------------------------------------------------------------------------
# Part 1: Load the Data
# ------------------------------------------------------------------------------
sce <- read10xCounts("raw_gene_bc_matrices/GRCh38", col.names = TRUE)
sce
dim(sce)

# This is the same as Chapter_08
# About 750k cells (mostly empty droplets of course)
totals <- colSums(counts(sce))
summary(totals)
table(totals == 0)
head(rowData(sce))

# ------------------------------------------------------------------------------
# Part 2: Gene annotation
# ------------------------------------------------------------------------------
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ID, rowData(sce)$Symbol)
head(rownames(sce))

library(EnsDb.Hsapiens.v86)
# Get chromosome locations: Thisis for mitochondrial gene annotation
location <- mapIds(EnsDb.Hsapiens.v86, keys=rowData(sce)$ID, 
                   column="SEQNAME", keytype="GENEID")
table(location)
# 13 Mitochondrial Genes

# ------------------------------------------------------------------------------
# Part 3: Cell Detection
# ------------------------------------------------------------------------------
# This is droplet-based 10x Genomics, so we detect empty drops and remove them
set.seed(100)
empty_drops <- emptyDrops(counts(sce))
table(empty_drops$FDR <= 0.001, useNA = "ifany")

keep <- !is.na(empty_drops$FDR) & empty_drops$FDR <= 0.001

# # Plot barcode ranks: Right upper part is what you want
# bar_out <- barcodeRanks(counts(sce))
# 
# plot(bar_out$rank, bar_out$total, log = "xy", 
#      xlab = "Rank",
#      ylab = "Total UMI count")
# 
# abline(h = metadata(bar_out)$knee, lty = 2)
# abline(h = metadata(bar_out)$inflection, lty = 2)
# 
# plot(bar_out$rank, bar_out$total, log = "xy", 
#      xlab = "Rank",
#      ylab = "Total UMI count")
# 
# is.cell <- !is.na(empty_drops$FDR) & empty_drops$FDR <= 0.001
# 
# points(bar_out$rank[is.cell],
#        bar_out$total[is.cell],
#        pch = 16,
#        cex = 0.5,
#        col="red")

# Filter out empty droplets
sce <- sce[, keep]

# Now we have a data structure ready for QC
sce
# 4402 cells passed the filter

# ------------------------------------------------------------------------------
# Part 4: Quality Control (QC)
# ------------------------------------------------------------------------------
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
# Part 5: Normalization
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
# Part 6: Quantifying per-gene variation (Variance Modeling)
# ------------------------------------------------------------------------------
# Poisson Modeling
set.seed(1001)
dec_pbmc <- modelGeneVarByPoisson(sce)
top_pbmc <- getTopHVGs(dec_pbmc, prop=0.1)

# Have a look at what you got
dec_pbmc[top_pbmc,]

library(scater)
# First, ill run conventional PCA with 50 PCs
sce <- runPCA(sce, subset_row=top_pbmc, ncomponents = 50)

percentVar <- attr(reducedDim(sce, "PCA"), "percentVar")

pdf("Scree_Plot.pdf", width = 6, height = 5)
plot(percentVar,
     type = "b",
     pch = 16,
     xlab = "Principal component",
     ylab = "Variance explained (%)")
dev.off()

barplot(percentVar[1:10],
        names.arg = seq_along(percentVar[1:10]),
        xlab = "Principal component",
        ylab = "Variance explained (%)",
        main = "PCA Scree Plot")

cumulative_sum <- cumsum(percentVar)

pdf("PC_Cummulative_Plot.pdf", width = 6, height = 5)
plot(seq_along(cumulative_sum),
     cumulative_sum,
     type = "b", 
     pch = 16,
     xlab = "Principal component",
     ylab = "Variance explained (%)",
     main = "Cummulative Sum")
dev.off()
# 5-6 PCs, the rest are useless

# Now lets use variance modeling to support it
sce <- denoisePCA(sce, subset.row=top_pbmc, technical=dec_pbmc)
ncol(reducedDim(sce, "PCA"))

set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")

set.seed(1000000)
sce <- runUMAP(sce, dimred="PCA")

reducedDimNames(sce)

plotReducedDim(sce, "TSNE")

# ------------------------------------------------------------------------------
# Graph-Based Clustering
# ------------------------------------------------------------------------------
# clustering is done by scran library
nn_clusters <- clusterCells(sce, use.dimred="PCA")
table(nn_clusters)

# Assign clusters to the data for visualization
colLabels(sce) <- nn_clusters

pdf("Cluster_Labeled_TSNE.pdf", width = 6, height = 6)
plotReducedDim(sce, "TSNE", colour_by="label")
dev.off()

# Change default parameters
library(bluster)
nn_clusters2 <- clusterCells(sce, use.dimred="PCA", 
                             BLUSPARAM=SNNGraphParam(k=10, type="rank", 
                                                     cluster.fun="walktrap"))
table(nn_clusters2)

# Get cluster info
nn_clust_info <- clusterCells(sce, use.dimred="PCA", full=TRUE)
nn_clust_info$objects$graph

set.seed(11000)
reducedDim(sce, "force") <- igraph::layout_with_fr(nn_clust_info$objects$graph)
plotReducedDim(sce, colour_by="label", dimred="force")




