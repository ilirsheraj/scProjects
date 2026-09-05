# Chapter 4: The SingleCellExperiment class
library(SingleCellExperiment)
library(BiocFileCache)

# To construct a rudimentary SingleCellExperiment object, we only need to fill the assays slot.
# Let's download a data set from ArrayExpress servers
bfc <- BiocFileCache("raw_data", ask = FALSE)
calero_counts <- bfcrpath(bfc, file.path("https://www.ebi.ac.uk/biostudies", 
                          "files/E-MTAB-5522/counts_Calero_20160113.tsv"))

mat <- read.delim(calero_counts, header=TRUE, row.names=1, check.names=FALSE)
mat[1:5, 1:5]

# Keep only endogenous genes for the moment, Set the spike matrix aside
spike_mat <- mat[grepl("^ERCC-", rownames(mat)),]
# Keep genes only
mat <- mat[grepl("^ENSMUSG", rownames(mat)),]
dim(mat)

# Remove gene length column to keep read counts only
gene_length <- mat[,1]
mat <- as.matrix(mat[,-1])
dim(mat)

# ------------------------------------------------------------------------------
# Now construct singlecellexperiment object
# ------------------------------------------------------------------------------
# # We provide our data as a named list where each entry of the list is a matrix
# sce <- SingleCellExperiment(assays = list(counts = mat))
# sce
# 
# # We can get the counts back:use getter functions assay() or counts()
# mat2 <- counts(sce)
# mat2[1:5, 1:5]
# 
# # Or using assay() function
# assay(sce, "counts")
# 
# # Now we can add more asssays
# # We can log-normalize the counts and they will be stored in assay automatically
# sce <- scuttle::logNormCounts(sce)
# sce
# 
# dim(logcounts(sce))
# 
# # We can add custom assays: Create a new matrix where you add 100 to all values
# counts_100 <- counts(sce) + 100
# assay(sce, "counts_100") <- counts_100
# 
# # Check all the assays available
# assays(sce)
# 
# # Remove the third assay
# assays(sce) <- assays(sce)[1:2]
# assays(sce)
# 
# # To get the assay names as a vector rather than list
# assayNames(sce)

# Since some functions of scuttle have been deprecated, im using scrapper here which
# is the current API and test to show the results are identical
# Make it more up to date using scrapper
library(scrapper)
sce2 <- SingleCellExperiment(assays = list(counts = mat))
sce2

# Calculate library-size-based size factors
sizeFactors(sce2) <- scrapper::centerSizeFactors(colSums(counts(sce2)))
# sizeFactors(sce2)

# Normalize and log-transform counts
logcounts(sce2) <- scrapper::normalizeCounts(
  counts(sce2),
  size.factors = sizeFactors(sce2),
  log = TRUE)

# identical(rownames(sce), rownames(sce2))
# identical(colnames(sce), colnames(sce2))
# 
# x <- as.matrix(logcounts(sce))
# y <- as.matrix(logcounts(sce2))
# 
# cor(as.vector(x), as.vector(y))
# 
# max(abs(x - y))
# mean(abs(x - y))
# 
# cor(sizeFactors(sce), sizeFactors(sce2))
# 
# summary(sizeFactors(sce))
# summary(sizeFactors(sce2))
# 
# all.equal(
#   as.matrix(logcounts(sce)),
#   as.matrix(logcounts(sce2))
# )
# The two methods are essentially the same, so time to move with scrapper

# We can add custom assays: Create a new matrix where you add 100 to all values
counts_100 <- counts(sce2) + 100
# assign a new entry to assays slot
assay(sce2, "counts_100") <- counts_100
# new assay has now been added.
assays(sce2)

# Remove the third assay
assays(sce2) <- assays(sce2)[1:2]
assays(sce2)

# To get the assay names as a vector rather than list
assayNames(sce2)

# ------------------------------------------------------------------------------
# Handling Metadata
# ------------------------------------------------------------------------------
# Download the SDRF file containing the metadata.
# dir.create("raw_data", showWarnings = FALSE)
download.file(
  url = "https://www.ebi.ac.uk/biostudies/files/E-MTAB-5522/E-MTAB-5522.sdrf.txt",
  destfile = "raw_data/E-MTAB-5522.sdrf.txt",
  mode = "wb"
)

coldata <- read.delim("raw_data/E-MTAB-5522.sdrf.txt", check.names = FALSE)
head(coldata)

# Filter to keep interesting columns only
coldata <- coldata[coldata[,"Derived Array Data File"]=="counts_Calero_20160113.tsv",]

# Only keeping interesting columns, and setting the library names as the row names.
coldata <- DataFrame(
  genotype=coldata[,"Characteristics[genotype]"],
  phenotype=coldata[,"Characteristics[phenotype]"],
  spike_in=coldata[,"Factor Value[spike-in addition]"],
  row.names=coldata[,"Source Name"]
)

# Different Options to add metadata
# 1 - Create a new singlecellexperiment object from scratch
sce <- SingleCellExperiment(assays = list(counts=mat), colData=coldata)
# Now it has assay() and colData()

# Accessing metadata is similar to accessing assay
colData(sce)

head(sce$spike_in)
head(sce$phenotype)
head(sce$genotype)

# 2 - Add to an existing object en masse
sce <- SingleCellExperiment(list(counts=mat))
colData(sce) <- coldata
sce

# 3 - Add to an existing object piece by piece
sce <- SingleCellExperiment(list(counts=mat))
sce$phenotype <- coldata$phenotype
colData(sce)

# Check that names in expression matrix and metadata are the same
stopifnot(identical(rownames(coldata), colnames(mat)))

# This is deprecated too
sce <- scuttle::addPerCellQC(sce)
colData(sce)

# Since scuttle is deprecated, use scrapper
sce2 <- SingleCellExperiment(assays = list(counts = mat))
colData(sce2) <- coldata

# Use this instead
sce2 <- scrapper::quickRnaQc.se(
  sce2,
  subsets = list()
)
colData(sce2)

# They are identical
# x <- as.matrix(colData(sce))
# y <- as.matrix(colData(sce2))
# 
# cor(as.vector(as.numeric((x[,"detected"]))), as.vector(as.numeric((y[,"detected"]))))
# cor(as.vector(as.numeric((x[,"sum"]))), as.vector(as.numeric((y[,"sum"]))))

# Add metadata to the assay (expression matrix): Manually set gene length to the assay matrix
rowData(sce2)$Length <- gene_length
rowData(sce2)

# We can use this function to get per feature quality control metrics
sce2 <- scuttle::addPerFeatureQC(sce2)
rowData(sce2)

# To get genomic coordinates: At the moment its empty, but it can be added
rowRanges(sce2)

# The paper uses ensembl 82 release
# Get the GTF file containing the Ensembl annotation used in this dataset
mm10_gtf <- bfcrpath(bfc, file.path("http://ftp.ensembl.org/pub/release-82",
                                    "gtf/mus_musculus/Mus_musculus.GRCm38.82.gtf.gz"))
# Now get the genomic ranges
gene_data <- rtracklayer::import(mm10_gtf)

# Keep genes only
mcols(gene_data)
# Check the column names: standard functions for working with GRanges class
colnames(mcols(gene_data))
# Keep gene only
gene_data <- gene_data[gene_data$type=="gene"]
# Ensemble
names(gene_data) <- gene_data$gene_id
# Keep only gene-related columns
is.gene.related <- grep("gene_", colnames(mcols(gene_data)))
mcols(gene_data) <- mcols(gene_data)[,is.gene.related]
gene_data

# Add it to sce object
rowRanges(sce2) <- gene_data[rownames(sce2)]

# Now we have row ranges
rowRanges(sce2)[1:10,]
rowData(sce2)

# We can add other metadata according to our wishes too: This is not synchronized
# with rows and columns. For synchronization use rowData() and colData()
my_genes <- c("gene_1", "gene_5")
metadata(sce2) <- list(favorite_genes = my_genes)
metadata(sce2)

your_genes <- c("gene_4", "gene_8")
metadata(sce2)$your_genes <- your_genes
metadata(sce2)

# ------------------------------------------------------------------------------
# Subsetting and Combining
# ------------------------------------------------------------------------------
# Everything is integrated
# Lets pick only the first 10 cells (columns)
first_10 <- sce2[,1:10]
ncol(counts(first_10))
colData(first_10)

table(sce2$phenotype)

# Subset only wild-type cells
wt_only <- sce2[, sce2$phenotype == "wild type phenotype"]
ncol(counts(wt_only))
colData(wt_only)

# We can do the same to features: Say we only want to keep protein-coding genes:
# First, lets see the identity of the "genes" starting from the most common
sort(table(rowData(sce2)$gene_biotype), decreasing = TRUE)

# Now ill keep protein-coding only
coding_only <- sce2[rowData(sce2)$gene_biotype == "protein_coding",]
nrow(counts(coding_only))
rowData(coding_only)

# Let's combine two sces to create a double of everything
sce3 <- cbind(sce2, sce2)
ncol(counts(sce3))
colData(sce3)

# Do the same for rows: double number of genes
sce3 <- rbind(sce2, sce2)
nrow(counts(sce3))

dim(assay(sce3))

# ------------------------------------------------------------------------------
# Single-cell-specific fields
# ------------------------------------------------------------------------------
# Dimensionality Reduction
# sce <- scater::logNormCounts(sce)

sce2 <- scrapper::normalizeRnaCounts.se(
  sce2,
  size.factors = colSums(counts(sce2)))

sce2 <- scater::runPCA(sce2)

dim(reducedDim(sce2, "PCA"))

sce2 <- scater::runTSNE(sce2, perplexity = 0.1)
head(reducedDim(sce2, "TSNE"))

# We can see the list of results in reducedDim
reducedDims(sce2)

# Manually add UMAP: Change this a bit 
u <- uwot::umap(
  t(as.matrix(logcounts(sce2))),
  n_neighbors = 2)
# Add it manually to reducedDim slot
reducedDim(sce2, "UMAP_uwot") <- u

reducedDims(sce2)

head(reducedDim(sce2, "UMAP_uwot"))

# ------------------------------------------------------------------------------
# Alternative Experiments: spike-in
# ------------------------------------------------------------------------------

# spike-in matrix contains a column with gene length by default in this case, so
# let's remove it first
spike_se <- SummarizedExperiment(list(counts=spike_mat[,-1]))
spike_se

# Than add it in alternative experiment slot
altExp(sce2, "spike") <- spike_se
altExps(sce2)

# It's in synch with the endogenous genes dataset
sub <- sce2[,1:2]
altExp(sub, "spike")

# Size factors: stored in colData
sce2 <- scran::computeSumFactors(sce2)
summary(sizeFactors(sce2))

# Add them manually
# sizeFactors(sce2) <- scater::librarySizeFactors(sce2)
# summary(sizeFactors(sce))
# Size factors already calculated in the beginning, but just in case
sizeFactors(sce2) <- scrapper::centerSizeFactors(
  colSums(counts(sce2)))
sizeFactors(sce2)

# Column Labels: Here i have changed the old code too
# colLabels(sce) <- scran::clusterCells(sce, use.dimred="PCA")
# table(colLabels(sce))
colLabels(sce2) <- bluster::clusterRows(
  reducedDim(sce2, "PCA"),
  BLUSPARAM = bluster::NNGraphParam())

table(colLabels(sce2))
