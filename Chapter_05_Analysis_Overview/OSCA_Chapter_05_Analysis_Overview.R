# Chapter 5: Analysis overview
## A quick "end-to-end" analysis to get a general idea
# Here, we use the a droplet-based retina dataset from Macosko et al that comes
# with scRNAseq package
# BiocManager::install("scRNAseq")
library(scRNAseq)
library(scater)
library(bluster)
library(ggplot2)
library(batchelor) # for batch correction
# ------------------------------------------------------------------------------
# Quick Start - Simple
# ------------------------------------------------------------------------------
# Load the data
sce <- MacoskoRetinaData()
sce

dim(sce)
# 24658 49300

colData(sce)
rowData(sce)

# parse matrix
assay(sce)[1:5, 1:5]

# We can carry out Quality control using mitochondrial genes.
is_mito <- grepl("^MT-", rownames(sce))
# There are 31 mitochondrial genes
sum(is_mito)
# 31 mitochondrial genes

# This is deprecated
# qcstats <- perCellQCMetrics(sce, subsets=list(Mito=is_mito))
# use scrapper instead
# The long, step-by-step way of doing this
qcstats <- scrapper::computeRnaQcMetrics(
  counts(sce),
  subsets = list(Mito = is_mito))
head(qcstats)

# define the threshold
thresholds <- scrapper::suggestRnaQcThresholds(qcstats)

# Create a TRUE/FALSE vector
filtered <- scrapper::filterRnaQcMetrics(
  thresholds,
  qcstats)

# Apply and remove filtereed out cells
sce_filtered <- sce[, filtered]
dim(sce_filtered)


# Quick and shorter way
sce.qc <- scrapper::quickRnaQc.se(
  sce,
  subsets = list(Mito = is_mito))

table(sce.qc$keep)

sce.filtered <- sce.qc[, sce.qc$keep]

dim(sce.filtered)
# 24658 45877

# Normalization: updated to scrapper
sce.filtered <- scrapper::normalizeRnaCounts.se(
  sce.filtered,
  size.factors = colSums(counts(sce.filtered)))

assayNames(sce.filtered)

# Feature selection.
# Deprecated too
# library(scran)
# dec <- modelGeneVar(sce.filtered)
# hvg <- getTopHVGs(dec, prop=0.1)

# Use scrapper instead
dec2 <- scrapper::modelGeneVariances(logcounts(sce.filtered))
head(dec2)

# If you want it formatted (check docummentation)
dec2.formatted <- scrapper::formatModelGeneVariancesResult(dec2)
head(dec2.formatted)

hvg <- scrapper::chooseHighlyVariableGenes(
  dec2.formatted$residuals,
  top = 2000)

# PCA
set.seed(1234)
sce.filtered <- runPCA(sce.filtered, ncomponents=25, subset_row=hvg)
reducedDimNames(sce.filtered)

# 25 was chosen by authors, so lets diagnose it a bit
var_exp <- attr(reducedDim(sce.filtered, "PCA"), "percentVar")

# Generate a scree plot and observe where it flattens
scree_df <- data.frame(
  PC = seq_along(var_exp),
  Variance = var_exp)

ggplot(scree_df, aes(x = PC, y = Variance)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq_along(var_exp)) +
  labs(
    x = "Principal Component",
    y = "Variance Explained (%)",
    title = "PCA Scree Plot") +
  theme_classic()
# Flattened by 8-9th component

# PLot cummulative variance
scree_df$Cumulative <- cumsum(scree_df$Variance)

ggplot(scree_df, aes(PC, Cumulative)) +
  geom_line() +
  geom_point() +
  scale_x_continuous(breaks = seq_along(var_exp)) +
  labs(
    x = "Principal Component",
    y = "Cumulative Variance Explained (%)",
    title = "Cumulative PCA Variance") +
  theme_classic()

# Visually inspect the 25 ones
var_exp <- attr(reducedDim(sce.filtered, "PCA"), "percentVar")

data.frame(
  PC = seq_along(var_exp),
  Variance = round(var_exp, 3),
  Cumulative = round(cumsum(var_exp), 3))

# Clustering: Also updated
colLabels(sce.filtered) <- clusterRows(reducedDim(sce.filtered, "PCA"),
                                       BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))

table(colLabels(sce.filtered))

# Visualization
set.seed(1234)
sce.filtered <- runUMAP(sce.filtered, dimred = "PCA")

plotUMAP(sce.filtered, colour_by = "label")

# Marker detection: Old way
# markers <- findMarkers(sce.filtered, test.type="wilcox", direction="up", lfc=1)

# New Way
markers <- scrapper::scoreMarkers.se(sce.filtered, colLabels(sce.filtered),
                                     more.marker.args = list(threshold = 1))

head(markers[[1]])

# Top 5 markers per each cluster
top5 <- lapply(markers, function(x) {
  rownames(x)[seq_len(min(5, nrow(x)))]})

top_genes <- unique(unlist(top5))
top_genes

expr <- as.matrix(logcounts(sce.filtered)[top.genes, ])
clusters <- colLabels(sce.filtered)

dot_data <- do.call(rbind, lapply(seq_along(top_genes), function(i) {
  gene <- top_genes[i]
  do.call(rbind, lapply(levels(factor(clusters)), function(cl) {
    cells <- clusters == cl
    
    data.frame(
      gene = gene,
      cluster = cl,
      avg_expr = mean(expr[i, cells]),
      pct_expr = mean(expr[i, cells] > 0) * 100
    )
  }))
}))

dot_data$cluster <- factor(dot_data$cluster, levels = as.character(1:11))

ggplot(dot_data, aes(x = cluster, y = gene, size = pct_expr, colour = avg_expr)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  labs(
    x = "Cluster",
    y = "Marker gene",
    size = "% expressed",
    colour = "Average\nexpression") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

dot_data$scaled_expr <- ave(dot_data$avg_expr, dot_data$gene,
                            FUN = function(x) as.numeric(scale(x)))

head(dot_data)

ggplot(dot_data,aes(x = cluster, y = gene, size = pct_expr, colour = scaled_expr)) +
  geom_point() +
  scale_size(range = c(1, 8)) +
  labs(
    x = "Cluster",
    y = "Marker gene",
    size = "% expressed",
    colour = "Scaled\nexpression"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1))

markers[["9"]][1:20, ]
markers[["10"]][1:20, ]
markers[["11"]][1:20, ]

sce.filtered$cluster_9_11 <- ifelse(
  colLabels(sce.filtered) %in% c("9", "10", "11"),
  as.character(colLabels(sce.filtered)),
  "Other")

plotUMAP(
  sce.filtered,
  colour_by = "cluster_9_11")

m9  <- rownames(markers[["9"]])[1:10]
m10 <- rownames(markers[["10"]])[1:10]
m11 <- rownames(markers[["11"]])[1:10]

genes.9.11 <- unique(c(m9, m10, m11))

# ------------------------------------------------------------------------------
# Part 2: Batches
# ------------------------------------------------------------------------------
# Pancreas Smart-seq2 dataset from Segerstolpe Data
sce <- SegerstolpePancreasData()
dim(sce)
sce

# Have a look at metadata
colData(sce)
table(colData(sce)$disease)
# Gene annotation
rowData(sce)

# This is deprecated
# # Quality control (using ERCCs).
# qcstats <- perCellQCMetrics(sce)
# head(qcstats)
# 
# filtered <- quickPerCellQC(qcstats, percent_subsets="altexps_ERCC_percent")
# sce2 <- sce[, !filtered$discard]
# 
# dim(sce2)

# Use this instead
sce <- scrapper::quickRnaQc.se(sce, subsets = list(), 
                               altexp.proportions = "ERCC")

table(sce$keep)

sce <- sce[, sce$keep]

dim(sce)

# Normalization: updated to scrapper
sce <- scrapper::normalizeRnaCounts.se(
  sce,
  size.factors = colSums(counts(sce)))

assayNames(sce)

# Feature selection, blocking on the individual of origin.
table(colData(sce)$individual)
dec <- scrapper::modelGeneVariances(logcounts(sce), block=sce$individual)
head(dec)
dec.formatted <- scrapper::formatModelGeneVariancesResult(dec)
head(dec.formatted)

n_hvg <- ceiling(0.1 * nrow(sce))
hvg <- scrapper::chooseHighlyVariableGenes(
  dec.formatted$residuals,
  top = n_hvg)

# Batch correction.
set.seed(1234)
sce <- correctExperiments(sce, batch=sce$individual, 
                          subset.row=hvg, 
                          correct.all=TRUE)

sce
reducedDimNames(sce)

colLabels(sce) <- clusterRows(
  reducedDim(sce, "corrected"),
  BLUSPARAM = NNGraphParam(cluster.fun = "louvain"))

table(colLabels(sce))

# Visualization
set.seed(1234)
sce <- runUMAP(sce, dimred = 'corrected')

gridExtra::grid.arrange(
  plotUMAP(sce, colour_by="label"),
  plotUMAP(sce, colour_by="individual"),
  ncol=2)

# Marker detection, blocking on the individual of origin.
markers <- scrapper::scoreMarkers.se(sce, colLabels(sce),
                                     more.marker.args = list(threshold = 1))

head(markers[[1]])

# Top 5 markers per each cluster
top5 <- lapply(markers, function(x) {
  rownames(x)[seq_len(min(5, nrow(x)))]})

top_genes <- unique(unlist(top5))
top_genes
