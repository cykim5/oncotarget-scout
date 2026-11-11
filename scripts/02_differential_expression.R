# 02_differential_expression.R
# DESeq2 analyses + GSEA.
#   A) Tumor vs normal                -> tumor specificity layer
#   B) STK11/KEAP1 mutant vs wild type -> subgroup enrichment layer
# Run after 01_download_tcga.R:  Rscript scripts/02_differential_expression.R

library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(msigdbr)

dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

se  <- readRDS("data/processed/tcga_luad_se.rds")
maf <- readRDS("data/processed/tcga_luad_maf.rds")

counts <- assay(se, "unstranded")
mode(counts) <- "integer"

# Patient barcode = first 12 characters of the TCGA sample barcode
patient <- substr(colnames(se), 1, 12)

# Mutation status per patient for the subgroup genes
mut_genes <- c("STK11", "KEAP1", "KRAS")
mut_patients <- unique(substr(maf$Tumor_Sample_Barcode, 1, 12)[
  maf$Hugo_Symbol %in% mut_genes])
status <- data.frame(
  patient = patient,
  stk11_keap1_mut = patient %in% unique(substr(maf$Tumor_Sample_Barcode, 1, 12)[
    maf$Hugo_Symbol %in% c("STK11", "KEAP1")])
)

sample_type <- colData(se)$sample_type

run_deseq <- function(keep, group, group_name) {
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, keep],
    colData   = data.frame(group = factor(group)),
    design    = ~ group
  )
  dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]  # basic low count filter
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, alpha = 0.05))
  res$gene <- rownames(res)
  out <- res[order(res$padj), ]
  write.csv(out, paste0("results/tables/de_", group_name, ".csv"),
            row.names = FALSE)
  out
}

# A) Tumor vs normal
keep_a <- sample_type %in% c("Primary Tumor", "Solid Tissue Normal")
de_tumor <- run_deseq(keep_a,
                      ifelse(sample_type[keep_a] == "Primary Tumor",
                             "tumor", "normal"),
                      "tumor_vs_normal")

# B) Mutant vs wild type, tumors only
keep_b <- sample_type == "Primary Tumor"
grp_b  <- ifelse(status$stk11_keap1_mut[keep_b], "mutant", "wildtype")
de_sub <- run_deseq(keep_b, grp_b, "stk11_keap1_mut_vs_wt")

# --- GSEA on the subgroup contrast (Hallmark) --------------------------------
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
term2gene <- hallmark[, c("gs_name", "gene_symbol")]

ranked <- de_sub$stat
names(ranked) <- de_sub$gene
ranked <- sort(ranked[!is.na(ranked)], decreasing = TRUE)

gsea_res <- GSEA(ranked, TERM2GENE = term2gene, pvalueCutoff = 0.25)
write.csv(as.data.frame(gsea_res),
          "results/tables/gsea_hallmark_subgroup.csv", row.names = FALSE)

message("DE tables and GSEA written to results/.")
message("Answer key positions are computed in 03 from the score table; ",
        "report them in results/answer_key_check.md after the first run.")
