# Reading Counts into R from tabular data
# BiocManager::install("BiocFileCache")

###################################################
# Part 1: Tabular Formats
###################################################
# To use simple read.del
library(BiocFileCache)
bfc <- BiocFileCache(ask=FALSE)
# Righ-click on ftp in geo and get the following link:
url <- file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
                 "GSE85nnn/GSE85241/suppl",
                 "GSE85241%5Fcellsystems%5Fdataset%5F4donors%5Fupdated%2Ecsv%2Egz")

# Download it in the local directory
muraro.fname <- bfcrpath(bfc, url)
local.name <- URLdecode(basename(url))
unlink(local.name)
file.symlink(muraro.fname, local.name)

# Read it into R
mat <- as.matrix(read.delim("GSE85241_cellsystems_dataset_4donors_updated.csv.gz"))
dim(mat)
mat[1:5, 1:10]

# Memory usage of matrix
object.size(mat)

# ------------------------------------------------------------------------------
# Read it in sparse format way: Doesn't store zeros
# ------------------------------------------------------------------------------
# BiocManager::install("scuttle")
library(scuttle)
sparse_mat <- readSparseCounts("GSE85241_cellsystems_dataset_4donors_updated.csv.gz")
dim(sparse_mat)
# See some entries like a normal matrix
as.matrix(sparse_mat[1:5, 1:10, drop = FALSE])

class(sparse_mat)

# Memory of sparse matrix
object.size(sparse_mat)

# ------------------------------------------------------------------------------
# When data is stored in Excel Format
# ------------------------------------------------------------------------------
bfc <- BiocFileCache("raw_data", ask=FALSE)
wilson_fname <- bfcrpath(bfc, file.path("ftp://ftp.ncbi.nlm.nih.gov/geo/series",
                "GSE61nnn/GSE61533/suppl/GSE61533_HTSEQ_count_results.xls.gz"))

# Unzip it
library(R.utils)
wilson_name2 <- "GSE61533_HTSEQ_count_results.xls"
gunzip(wilson_fname, destname=wilson_name2, remove=FALSE, overwrite=TRUE)

all_counts <- readxl::read_excel("GSE61533_HTSEQ_count_results.xls")
head(all_counts)

gene.names <- all_counts$ID
all_counts <- as.matrix(all_counts[,-1])
rownames(all_counts) <- gene.names
dim(all_counts)

################################################
# Part 2: Cellranger Output
################################################
# Reading data from Cellranger Output
# BiocManager::install("DropletTestFiles")
library(DropletTestFiles)
# cached <- getTestFile("tenx-2.1.0-pbmc4k/1.0.0/filtered.tar.gz")
# fpath <- "tenx-2.1.0-pbmc4k"
# untar(cached, exdir=fpath)

# Download filtered data directly from 10x repository
fpath <- "pbmc4k_filtered_gene_bc_matrices.tar.gz"
untar(fpath)
# BiocManager::install("DropletUtils")
library(DropletUtils)
sce <- read10xCounts("filtered_gene_bc_matrices/GRCh38")
sce

# We can also read from different experiments: In this case same twice
dirA <- "filtered_gene_bc_matrices/GRCh38"
dirB <- "filtered_gene_bc_matrices/GRCh38"
sce <- read10xCounts(c(dirA, dirB))
sce
################################################################################
# From HDF5-based formats (Hierarchical Data Format version 5 (HDF5))
# BiocManager::install("zellkonverter")
library(zellkonverter)
demo <- system.file("extdata", "krumsiek11.h5ad", package = "zellkonverter")
sce <- readH5AD(demo)
sce

# BiocManager::install("LoomExperiment")
library(LoomExperiment)
demo <- system.file("extdata", "L1_DRG_20_example.loom", package = "LoomExperiment")
scle <- import(demo, type="SingleCellLoomExperiment")
scle

