# Part 2: Basics - Nornalization
library(scRNAseq)
library(scater)
library(org.Mm.eg.db)

sce_zeisel <- ZeiselBrainData()
sce_zeisel
head(rownames(sce_zeisel))

# Although it issues a warning and suggestion to use scrapper::... it DNE!
sce_zeisel <- aggregateAcrossFeatures(sce_zeisel, 
                                      id=sub("_loc[0-9]+$", "", rownames(sce_zeisel)))

# Annotate
rowData(sce_zeisel)$Ensembl <- mapIds(org.Mm.eg.db, keys=rownames(sce_zeisel), 
                                      keytype="SYMBOL", column="ENSEMBL")

rowData(sce_zeisel)
table(rowData(sce_zeisel)$featureType)
altExpNames(sce_zeisel)

# Quality Control
mito_genes <- rowData(sce_zeisel)$featureType == "mito"
sum(mito_genes)

sce_zeisel <- scrapper::quickRnaQc.se(
  sce_zeisel,
  subsets = list(Mito = mito_genes),
  altexp.proportions = "ERCC")

colData(sce_zeisel)
table(sce_zeisel$keep)
sce_zeisel_new <- sce_zeisel[, sce_zeisel$keep]
dim(sce_zeisel_new)

# Quality Control the older way
stats <- perCellQCMetrics(sce_zeisel, subsets=list(
  Mt=rowData(sce_zeisel)$featureType=="mito"))
qc <- quickPerCellQC(stats, percent_subsets=c("altexps_ERCC_percent", 
                                              "subsets_Mt_percent"))
sce.zeisel <- sce_zeisel[,!qc$discard]
dim(sce.zeisel)

# Old and New QC decisions
old_keep <- !qc$discard
new_keep <- sce_zeisel$keep

table(old_keep, new_keep)
different <- which(!old_keep & new_keep)
length(different)

# Look at the discrepancy between the two: High ERC Percent
stats[different, ]
different_metrics <- as.data.frame(stats[different,])
ggplot(different_metrics, aes(altexps_ERCC_percent, subsets_Mt_percent)) +
  geom_point()

summary(stats[!qc$discard, "altexps_ERCC_percent"])
summary(stats[sce_zeisel$keep, "altexps_ERCC_percent"])
different <- which(qc$discard & sce_zeisel$keep)

stats[different,
      c("sum",
        "detected",
        "subsets_Mt_percent",
        "altexps_ERCC_percent")]
metadata(sce_zeisel)$qc$thresholds
# For the time being i'll go with the new scrapper API

# ------------------------------------------------------------------------------
# Library size normalization
# ------------------------------------------------------------------------------
# The older method from OSCA
# library(scater)
# lib.sf.zeisel <- librarySizeFactors(sce_zeisel)
# summary(lib.sf.zeisel)

# Modern Method: Same as old
lib_sf_zeisel <- scrapper::centerSizeFactors(
  colSums(counts(sce_zeisel_new)))
summary(lib_sf_zeisel)

# Visualize it
hist(log10(lib_sf_zeisel), xlab="Log10[Size factor]", 
    main = "Library Size Factor Histogram",  col='grey80')

# ------------------------------------------------------------------------------
# Normalization by Deconvolution
# ------------------------------------------------------------------------------
library(scran)
set.seed(100)
clust_zeisel <- quickCluster(sce_zeisel_new) 
table(clust_zeisel)

# Deconvolution
deconv_sf_zeisel <- calculateSumFactors(sce_zeisel_new, cluster=clust_zeisel)
summary(deconv_sf_zeisel)

plot(lib_sf_zeisel, deconv_sf_zeisel, xlab="Library size factor",
     ylab="Deconvolution size factor", log='xy', pch=16,
     col=as.integer(factor(sce_zeisel_new$level1class)))
abline(a=0, b=1, col="red")

# ------------------------------------------------------------------------------
# Normalization by spike-ins
# ------------------------------------------------------------------------------
# Another DataSet
sce_richard <- RichardTCellData()
sce_richard
colData(sce_richard)
dim(sce_richard)

# Quality already determined
sce_richard <- sce_richard[,sce_richard$`single cell quality`=="OK"]
dim(sce_richard)
altExpNames(sce_richard)

# Compute spike-in factor
# ## Old method
# sce.richard <- computeSpikeFactors(sce.richard, "ERCC")
# summary(sizeFactors(sce.richard))

# New Method with scrapper
ercc_counts <- counts(altExp(sce_richard, "ERCC"))
spike_sf <- scrapper::centerSizeFactors(
  colSums(ercc_counts))

summary(spike_sf)
sizeFactors(sce_richard) <- spike_sf

sce_richard
head(colData(sce_richard))

to_plot <- data.frame(
  DeconvFactor=calculateSumFactors(sce_richard),
  SpikeFactor=spike_sf,
  Stimulus=sce_richard$stimulus, 
  Time=sce_richard$time
)

ggplot(to_plot, aes(x=DeconvFactor, y=SpikeFactor, color=Time)) +
  geom_point() + 
  facet_wrap(~Stimulus) + 
  scale_x_log10() + 
  scale_y_log10() + 
  geom_abline(intercept=0, slope=1, color="red")

# Differences between SPIKE-In and deconvolution seize factors
sce_richard_deconv <- sce_richard

logcounts(sce_richard_deconv) <- scrapper::normalizeCounts(
  counts(sce_richard_deconv),
  size.factors = to_plot$DeconvFactor)

sce_richard_spike <- sce_richard

logcounts(sce_richard_spike) <- scrapper::normalizeCounts(
  counts(sce_richard_spike),
  size.factors = to_plot$SpikeFactor)

gridExtra::grid.arrange(
  plotExpression(sce_richard_deconv, x="stimulus", 
                 colour_by="time", features="ENSMUSG00000092341") + 
    theme(axis.text.x = element_text(angle = 90)) + 
    ggtitle("After deconvolution"),
  
  plotExpression(sce_richard_spike, x="stimulus", 
                 colour_by="time", features="ENSMUSG00000092341") + 
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After spike-in normalization"),
  ncol=2
)

# ------------------------------------------------------------------------------
# Scaling and log-transforming
# ------------------------------------------------------------------------------
set.seed(100)
clust_zeisel <- quickCluster(sce_zeisel_new) 
sce_zeisel <- computeSumFactors(sce_zeisel_new, cluster=clust_zeisel, min.mean=0.1)
summary(sizeFactors(sce_zeisel))

sce_zeisel <- scrapper::normalizeRnaCounts.se(
  sce_zeisel, 
  size.factors = sizeFactors(sce_zeisel))

sce_zeisel
assayNames(sce_zeisel)

g6pd <- grep("^G6pdx$", rownames(sce_zeisel),
             ignore.case = TRUE, value = TRUE)

g6pd_expr <- logcounts(sce_zeisel)[g6pd, ]

df <- data.frame(
  expression = g6pd_expr,
  tissue = sce_zeisel$tissue)

# Fraction expressing in each tissue
tapply(df$expression > 0, df$tissue, mean)

ggplot(subset(df, expression > 0),
       aes(x = tissue, y = expression)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  labs(
    x = "Tissue",
    y = "Log-normalized G6PD expression",
    title = "G6PD expression among expressing cells") +
  scale_x_discrete(labels = c(
    ca1hippocampus = "CA1 hippocampus",
    sscortex = "SS cortex")) +
  theme_classic(base_size = 13)

# EOF