# Part 2: Basics - Feature Selection Dataset 2
library(scRNAseq)
sce_416b <- LunSpikeInData(which="416b")
sce_416b

head(colData(sce_416b))
head(rowData(sce_416b))

# ------------------------------------------------------------------------------
# Explore the Data
# ------------------------------------------------------------------------------
# We explored this before
for (nm in colnames(colData(sce_416b))) {
  cat("---", nm, "---")
  print(table(colData(sce_416b)[[nm]], useNA = "ifany"))
}

# Phenotype and block have equal numbers, check whether they are confounded
table(phenotype = colData(sce_416b)$phenotype,
      block = colData(sce_416b)$block)

# Convert block into factor
sce_416b$block <- factor(sce_416b$block)

# Annotate the genes: Ensemble to Gene Symbol
library(AnnotationHub)
ens_mm_v97 <- AnnotationHub()[["AH73905"]]
rowData(sce_416b)$ENSEMBL <- rownames(sce_416b)

rowData(sce_416b)$SYMBOL <- mapIds(ens_mm_v97, keys=rownames(sce_416b),
                                   keytype="GENEID", column="SYMBOL")

rowData(sce_416b)$SEQNAME <- mapIds(ens_mm_v97, keys=rownames(sce_416b),
                                    keytype="GENEID", column="SEQNAME")

head(rowData(sce_416b))



library(scater)
rownames(sce_416b) <- uniquifyFeatureNames(rowData(sce_416b)$ENSEMBL, 
                                           rowData(sce_416b)$SYMBOL)

# ------------------------------------------------------------------------------
# Quality Control
# ------------------------------------------------------------------------------
# There is Spike-in but we use only mitochondrial percentage in this case
mito <- which(rowData(sce_416b)$SEQNAME=="MT")

sce_416b <- scrapper::quickRnaQc.se(
  sce_416b,
  subsets = list(Mito = mito))

table(sce_416b$keep)

qc.metrics <- scrapper::computeRnaQcMetrics(
  counts(sce_416b),
  subsets = list(Mito = mito))

qc.thresholds <- scrapper::suggestRnaQcThresholds(qc.metrics)

# The cut-offs
qc.thresholds

# Remove low-quality cells
sce_416b <- sce_416b[, sce_416b$keep]
sce_416b


# ------------------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------------------
library(scran)
sce_416b <- computeSumFactors(sce_416b)
sce_416b <- scrapper::normalizeRnaCounts.se(
  sce_416b, 
  size.factors = sizeFactors(sce_416b))

sce_416b

# ------------------------------------------------------------------------------
# Quantify Technical Noise
# ------------------------------------------------------------------------------
ercc <- altExp(sce_416b, "ERCC")

# ERCC library sizes
ercc_sf <- colSums(counts(ercc))

# Center the ERCC size factors
ercc_sf <- scrapper::centerSizeFactors(ercc_sf)

# Normalize ERCC counts
logcounts(ercc) <- scrapper::normalizeCounts(
  counts(ercc),
  size.factors = ercc_sf
)

# Put it back
altExp(sce_416b, "ERCC") <- ercc

assayNames(altExp(sce_416b, "ERCC"))

# Endogenous genes
dec_gene <- scrapper::modelGeneVariances(logcounts(sce_416b))
dec_gene <- scrapper::formatModelGeneVariancesResult(dec_gene)


# ERCC spike-ins
dec_ercc <- scrapper::modelGeneVariances(logcounts(altExp(sce_416b, "ERCC")))
dec_ercc <- scrapper::formatModelGeneVariancesResult(dec_ercc)

# Fit technical variance trend to ERCCs
fit_ercc <- scrapper::fitVarianceTrend(dec_ercc$means,
                                       dec_ercc$variances)

# Make interpolation function for the ERCC trend
trend_ercc <- approxfun(dec_ercc$means,
                        fit_ercc$fitted,
                        rule = 2)

# Estimate technical variance for every endogenous gene
technical <- trend_ercc(dec_gene$means)

# Biological/excess variance
dec_spike <- dec_gene
dec_spike$fitted <- technical
dec_spike$residuals <- dec_spike$variances - technical

# Rank genes by biological variance
dec_spike[order(dec_spike$residuals, decreasing = TRUE),]

# Plot the whole thing
plot(dec_spike$means,
     dec_spike$variances,
     xlab = "Mean of log-expression",
     ylab = "Variance of log-expression")

points(dec_ercc$means,
       dec_ercc$variances,
       col = "red",
       pch = 16)

# ERCC-derived technical trend
ord <- order(dec_ercc$means)

lines(dec_ercc$means[ord],
      fit_ercc$fitted[ord],
      col = "dodgerblue",
      lwd = 2)

# x <- seq(min(dec_ercc$means),
#          max(dec_ercc$means),
#          length.out = 1000)
# 
# lines(x, trend_ercc(x), col = "dodgerblue", lwd = 2)

# ------------------------------------------------------------------------------
# Include Batch Effects
# ------------------------------------------------------------------------------
# Model endogenous genes within each block
gene_block_raw <- scrapper::modelGeneVariances(
  logcounts(sce_416b),
  block = sce_416b$block)

# Model ERCCs within each block
ercc_block_raw <- scrapper::modelGeneVariances(
  logcounts(altExp(sce_416b, "ERCC")),
  block = sce_416b$block)

block_ids <- gene_block_raw$block.ids

block_results <- vector("list", length(block_ids))
names(block_results) <- block_ids

block_ids <- gene_block_raw$block.ids

block_results <- vector("list", length(block_ids))
names(block_results) <- block_ids

for (i in seq_along(block_ids)) {
  
  # Already formatted
  g <- gene_block_raw$per.block[[i]]
  e <- ercc_block_raw$per.block[[i]]
  
  # Fit technical trend using ERCCs in this block
  fit <- scrapper::fitVarianceTrend(
    e$means,
    e$variances
  )
  
  # Turn fitted ERCC trend into a function
  trend <- approxfun(
    e$means,
    fit$fitted,
    rule = 2
  )
  
  # Evaluate ERCC-derived technical variance
  # at endogenous gene means
  tech <- trend(g$means)
  
  # Replace gene-wise fitted trend with ERCC-derived trend
  g$fitted <- tech
  
  # Equivalent of old "bio"
  g$residuals <- g$variances - tech
  
  block_results[[i]] <- g
}

head(block_results[[1]])
head(block_results[[2]])

dec_block <- block_results[[1]]

for (nm in c("means", "variances", "fitted", "residuals")) {
  dec_block[[nm]] <- Reduce(
    "+",
    lapply(block_results, function(x) x[[nm]])
  ) / length(block_results)
}


head(dec_block)

head(dec_block[order(dec_block$residuals, decreasing = TRUE),])

par(mfrow = c(1, 2))

for (i in seq_along(block_results)) {
  
  current <- block_results[[i]]
  ercc_current <- ercc_block_raw$per.block[[i]]
  
  fit <- scrapper::fitVarianceTrend(
    ercc_current$means,
    ercc_current$variances)
  
  trend <- approxfun(
    ercc_current$means,
    fit$fitted,
    rule = 2)
  
  plot(current$means,
       current$variances,
       main = gene_block_raw$block.ids[i],
       pch = 16,
       cex = 0.5,
       xlab = "Mean of log-expression",
       ylab = "Variance of log-expression")
  
  points(ercc_current$means,
         ercc_current$variances,
         col = "red",
         pch = 16)
  
  x <- seq(min(ercc_current$means),
           max(ercc_current$means),
           length.out = 500)
  
  lines(x, trend(x), col = "dodgerblue", lwd = 2)
}

par(mfrow = c(1, 1))

# ------------------------------------------------------------------------------
# Feature Selection
# ------------------------------------------------------------------------------
library(scater)
# Select HVGs using biological/residual variance

# Keep genes with positive biological variance
keep_var <- dec_block$residuals > 0

# Rank them by biological variance
ranked <- order(dec_block$residuals,decreasing = TRUE)

ranked <- ranked[keep_var[ranked]]

# Take top 10% of genes as HVGs
n_hvg <- round(0.1 * nrow(dec_block))

n_hvg <- min(n_hvg, length(ranked))

hvg <- rownames(dec_block)[ranked[seq_len(n_hvg)]]

length(hvg)
head(hvg)

# 2. PCA on normalized expression of HVGs
set.seed(100)

sce_416b <- runPCA(sce_416b, subset_row = hvg, ncomponents = 25,
                   exprs_values = "logcounts")

reducedDimNames(sce_416b)

# PCA by batch
pca_colors <- c("#0072B2", "#D55E00")

plotPCA(sce_416b,colour_by = "block") + 
  scale_colour_manual(values = pca_colors)

plotPCA(sce_416b, colour_by = "phenotype") +
  scale_colour_manual(values = pca_colors)

plotPCA(sce_416b, colour_by = "phenotype", shape_by = "block") +
  scale_colour_manual(values = pca_colors)

# Variance explained by each PC
percent_var <- attr(reducedDim(sce_416b, "PCA"), "percentVar")

# Scree plot
plot(seq_along(percent_var),
     percent_var,
     type = "b", 
     pch = 16,
     xlab = "Principal component",
     ylab = "Variance explained (%)",
     main = "PCA Scree Plot")

barplot(
  percent_var,
  names.arg = seq_along(percent_var),
  xlab = "Principal component",
  ylab = "Variance explained (%)",
  main = "PCA Scree Plot"
)

# EOF