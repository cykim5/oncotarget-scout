# 01_download_tcga.R
# Download TCGA-LUAD RNA-seq counts, somatic mutations, and clinical data.
# Run from the repository root:  Rscript scripts/01_download_tcga.R
# Outputs land in data/processed/ (git ignored).

library(TCGAbiolinks)
library(SummarizedExperiment)

dir.create("data/processed", showWarnings = FALSE, recursive = TRUE)

# --- RNA-seq counts (tumor + solid tissue normal) ---------------------------
query_exp <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Transcriptome Profiling",
  data.type     = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_exp)
se <- GDCprepare(query_exp)

# Keep the unstranded count assay; DESeq2 runs in 02_differential_expression.R
saveRDS(se, "data/processed/tcga_luad_se.rds")

# --- Somatic mutations (for STK11 / KEAP1 / KRAS status) ---------------------
query_mut <- GDCquery(
  project       = "TCGA-LUAD",
  data.category = "Simple Nucleotide Variation",
  data.type     = "Masked Somatic Mutation"
)
GDCdownload(query_mut)
maf <- GDCprepare(query_mut)
saveRDS(maf, "data/processed/tcga_luad_maf.rds")

# --- Clinical ----------------------------------------------------------------
clinical <- GDCquery_clinic("TCGA-LUAD", type = "clinical")
saveRDS(clinical, "data/processed/tcga_luad_clinical.rds")

message("Done. Samples: ", ncol(se),
        " | Tumor: ", sum(colData(se)$sample_type == "Primary Tumor"),
        " | Normal: ", sum(colData(se)$sample_type == "Solid Tissue Normal"))
