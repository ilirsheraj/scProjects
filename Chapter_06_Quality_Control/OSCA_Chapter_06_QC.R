# Part 2 - Basics: Quality Control
# scRNA-seq dataset from Lun et al. (2017), provided with no prior QC.
library(scRNAseq)
library(SingleCellExperiment)

sce_416b <- LunSpikeInData(which="416b")
sce_416b
head(rowData(sce_416b))
head(colData(sce_416b))

# Get a quick summary of all metadata
for (nm in colnames(colData(sce_416b))) {
  cat("---", nm, "---")
  print(table(colData(sce_416b)[[nm]], useNA = "ifany"))
}

# Phenotype and block have equal numbers, check whether they are confounded
table(phenotype = colData(sce_416b)$phenotype,
      block = colData(sce_416b)$block)
# All Balanced

# Include spike-in as well
with(
  as.data.frame(colData(sce_416b)),
  table(block, phenotype, spike.in.addition))
# Perfectly balanced
# technical/experimental block is block itself :)

sce_416b$block <- factor(sce_416b$block)
sce_416b

# Save for now, no need to download again, the load it back
# saveRDS(sce_416b, file = "sce_416b.rds")

# ------------------------------------------------------------------------------
# Different kinds of QC Metrics
# ------------------------------------------------------------------------------
# Identifying the mitochondrial transcripts in our SingleCellExperiment.
location <- rowRanges(sce_416b)
head(location)

# Detect Mitochondrial genes
is_mito <- any(seqnames(location)=="MT")
sum(is_mito)

# # This has been deprecated
# library(scuttle)
# df <- perCellQCMetrics(sce_416b, subsets=list(Mito=is.mito))
# summary(df$sum)

# Use scrapper instead
sce_416b <- scrapper::quickRnaQc.se(
  sce_416b,
  subsets = list(Mito = is_mito),
  altexp.proportions = "ERCC",
  block = sce_416b$block)

colnames(colData(sce_416b))

summary(sce_416b$sum)
summary(sce_416b$detected)
summary(sce_416b$subset.proportion.Mito)
summary(sce_416b$subset.proportion.ERCC)

# Further exploration
altExp(sce_416b, "ERCC")
# 92 ERCC
altExp(sce_416b, "SIRV")
# 7 SIRV
head(rownames(altExp(sce_416b, "ERCC")))
rownames(altExp(sce_416b, "SIRV"))

# ------------------------------------------------------------------------------
# Identifying low-quality cells
# ------------------------------------------------------------------------------
# Method 1: Fixed Thresholds
qc.lib <- sce_416b$sum < 1e5
qc.nexprs <- sce_416b$detected < 5e3
qc.spike <- sce_416b$subset.proportion.ERCC > 0.1
qc.mito <- sce_416b$subset.proportion.Mito > 0.1
discard <- qc.lib | qc.nexprs | qc.spike | qc.mito
sum(discard)

# Summarize the number of cells removed for each reason as well as the final result
DataFrame(LibSize=sum(qc.lib), NExprs=sum(qc.nexprs),
          SpikeProp=sum(qc.spike), MitoProp=sum(qc.mito), 
          Total=sum(discard))
# 39 removed in total

# Method 2: Adaptive Thresholds
# # This is deprecated
# reasons <- perCellQCFilters(df,
#                             sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
# colSums(as.matrix(reasons))

# scrapper has already provided the decision points, but if you wanna see them, here they are
table(sce_416b$keep)

# Look at all rejected sample metrics
colData(sce_416b)[!sce_416b$keep, ]

# If you want to see specific columns
colData(sce_416b)[!sce_416b$keep,
                  c("sum", "detected", "subset.proportion.Mito", "subset.proportion.ERCC")]

# Or if you want to see one column
colData(sce_416b)[!sce_416b$keep, "detected"]

# Method 3: Robust Methods
stats <- cbind(log10(sce_416b$sum), log10(sce_416b$detected),
               sce_416b$subset.proportion.Mito, sce_416b$subset.proportion.ERCC)

library(robustbase)
outlying <- adjOutlyingness(stats, only.outlyingness = TRUE)
multi.outlier <- isOutlier(outlying, type = "higher")
summary(multi.outlier)

# ------------------------------------------------------------------------------
# Checking diagnostic plots
# ------------------------------------------------------------------------------
# To make the names cleaner
sce_416b$phenotype <- ifelse(grepl("induced", sce_416b$phenotype), "induced", "wild type")
library(scater)

# Make it visible, the online one is dull
qc_colors <- scale_colour_manual(values = c(
  "TRUE"  = "#0072B2",
  "FALSE" = "#D55E00"),
  labels = c(
    "TRUE"  = "Keep",
    "FALSE" = "Discard")
  )

gridExtra::grid.arrange(
  plotColData(
    sce_416b, x="block", y="sum", colour_by="keep",
    other_fields="phenotype"
  ) +
    facet_wrap(~phenotype) +
    scale_y_log10() +
    qc_colors +
    ggtitle("Total count"),
  
  plotColData(
    sce_416b, x="block", y="detected", colour_by="keep", 
    other_fields="phenotype"
  ) +
    facet_wrap(~phenotype) +
    scale_y_log10() +
    qc_colors +
    ggtitle("Detected features"),
  
  plotColData(
    sce_416b, x="block", y="subset.proportion.Mito", 
    colour_by="keep", other_fields="phenotype"
  ) +
    facet_wrap(~phenotype) +
    qc_colors +
    ggtitle("Mito proportion"),
  
  plotColData(
    sce_416b, x="block", y="subset.proportion.ERCC", 
    colour_by="keep", other_fields="phenotype"
  ) +
    facet_wrap(~phenotype) +
    qc_colors +
    ggtitle("ERCC proportion"),
  
  ncol = 1
)


plotColData(sce_416b, x="sum", y="subset.proportion.Mito", colour_by="keep")

