# Part 2: Basics - Dimensionality Reduction
# The sama dataset was used to run a similar analysis before as well
library(scRNAseq)
library(scater)
library(org.Mm.eg.db)

# ------------------------------------------------------------------------------
# Load the Data
# ------------------------------------------------------------------------------
sce_zeisel <- ZeiselBrainData()
sce_zeisel
head(rownames(sce_zeisel))

sce_zeisel <- aggregateAcrossFeatures(sce_zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce_zeisel)))

# Annotate: Convert gene symbols into ENSEMBLE IDs
# This is mouse dataset: look at gene symbols above
rowData(sce_zeisel)$Ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce_zeisel), 
                                      keytype="SYMBOL", column="ENSEMBL")

rowData(sce_zeisel)
table(rowData(sce_zeisel)$featureType)
altExpNames(sce_zeisel)

# ------------------------------------------------------------------------------
# QC
# ------------------------------------------------------------------------------
# Quality Control: use both mitochondria and Spike-In
mito_genes <- rowData(sce_zeisel)$featureType == "mito"
sum(mito_genes)

sce_zeisel <- scrapper::quickRnaQc.se(
  sce_zeisel,
  subsets = list(Mito = mito_genes),
  altexp.proportions = "ERCC")

colData(sce_zeisel)
table(sce_zeisel$keep)
sce_zeisel <- sce_zeisel[, sce_zeisel$keep]
dim(sce_zeisel)

# ------------------------------------------------------------------------------
# Normalization
# ------------------------------------------------------------------------------
library(scran)
set.seed(100)
clust_zeisel <- quickCluster(sce_zeisel) 
table(clust_zeisel)

# Deconvolution
sce_zeisel <- computeSumFactors(sce_zeisel, cluster=clust_zeisel)
summary(sizeFactors(sce_zeisel))

sce_zeisel <- scrapper::normalizeRnaCounts.se(
  sce_zeisel, 
  size.factors = sizeFactors(sce_zeisel))

sce_zeisel
assayNames(sce_zeisel)

# Suse spike-in to remove technical noise
ercc <- altExp(sce_zeisel, "ERCC")

# ERCC library sizes
ercc_sf <- colSums(counts(ercc))

# Center the ERCC size factors
ercc_sf <- scrapper::centerSizeFactors(ercc_sf)

# Normalize ERCC counts
logcounts(ercc) <- scrapper::normalizeCounts(
  counts(ercc),
  size.factors = ercc_sf)

# Put it back in the sc class
altExp(sce_zeisel, "ERCC") <- ercc

assayNames(altExp(sce_zeisel, "ERCC"))

# Endogenous genes
dec_gene <- scrapper::modelGeneVariances(logcounts(sce_zeisel))
dec_gene <- scrapper::formatModelGeneVariancesResult(dec_gene)

# ERCC spike-ins
dec_ercc <- scrapper::modelGeneVariances(logcounts(altExp(sce_zeisel, "ERCC")))
dec_ercc <- scrapper::formatModelGeneVariancesResult(dec_ercc)

# Fit technical variance trend to ERCCs
fit_ercc <- scrapper::fitVarianceTrend(dec_ercc$means,
                                       dec_ercc$variances)

# Make interpolation function for the ERCC trend
trend_ercc <- approxfun(dec_ercc$means,
                        fit_ercc$fitted,
                        rule = 2)

# Estimate technical variance for every endogenous gene
dec_gene$tech <- trend_ercc(dec_gene$means)

# Biological/excess variance
dec_gene$bio <- dec_gene$variances - dec_gene$tech
head(dec_gene)

# Rank genes by biological variance
dec_gene[order(dec_gene$residuals, decreasing = TRUE),]

# Choose top 2000 Genes
top_zeisel <- rownames(sce_zeisel)[order(dec_gene$bio, decreasing = TRUE)[1:2000]]

# Run the PCA on these genes without removing the rest from the matrix
set.seed(100)
sce_zeisel <- runPCA(sce_zeisel, subset_row = top_zeisel, ncomponents = 50)

dim(reducedDim(sce_zeisel, "PCA"))

reducedDimNames(sce_zeisel)

# Plot the PCA
cell_colors <- c(
  "astrocytes_ependymal" = "#0066CC",
  "endothelial_mural"    = "#FF8C00",
  "interneurons"         = "#00A651",
  "microglia"            = "#E31A1C",
  "oligodendrocytes"     = "#6A00A8",
  "pyramidal CA1"        = "#7F4F24",
  "pyramidal SS"         = "#FF1493")

plotReducedDim(
  sce_zeisel,
  dimred = "PCA",
  colour_by = "level1class",
  point_size = 1.5
) +
  scale_colour_manual(values = cell_colors)

percent_var <- attr(reducedDim(sce_zeisel, "PCA"), "percentVar")

# Scree plot
plot(seq_along(percent_var),
     percent_var,
     type = "b", 
     pch = 16,
     xlab = "Principal component",
     ylab = "Variance explained (%)",
     main = "PCA Scree Plot")

barplot(
  percent_var[1:20],
  names.arg = seq_along(percent_var[1:20]),
  xlab = "Principal component",
  ylab = "Variance explained (%)",
  main = "PCA Scree Plot")

# Plot PCAs for the first 4 components
plotReducedDim(sce_zeisel, 
               dimred="PCA", 
               ncomponents=4,
               colour_by="level1class") +
  scale_colour_manual(values = cell_colors)
  
# TCB