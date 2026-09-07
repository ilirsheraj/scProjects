# Part 2 - Basics: Quality Control
# scRNA-seq dataset from Lun et al. (2017), provided with no prior QC.
library(scRNAseq)

sce_416b <- LunSpikeInData(which="416b")
sce_416b
head(rowData(sce_416b))
head(colData(sce_416b))

# Get a quick summary of all metadata
for (nm in colnames(colData(sce_416b))) {
  cat("---", nm, "---")
  print(table(cd[[nm]], useNA = "ifany"))
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

# Save for now, no need to download again
saveRDS(sce_416b, file = "sce_416b.rds")


