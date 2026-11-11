# Pipeline scripts

Run in number order. R steps handle TCGA and pathway analysis. Python steps handle single cell, DepMap joins, and score assembly. Paths are relative to the repository root.

| Script | Step | Output |
| --- | --- | --- |
| `01_download_tcga.R` | Download TCGA-LUAD counts, mutations, clinical data | `data/processed/tcga_luad_*.rds` |
| `02_differential_expression.R` | DESeq2: tumor vs normal, then STK11/KEAP1 mutant vs wild type; GSEA with clusterProfiler | `results/tables/de_*.csv`, `results/figures/` |
| `03_survival_and_scoring.py` | Cox survival layer, DepMap and GTEx joins, final score from docs/target-score.md | `results/tables/target_scores.csv` |
| `04_singlecell_layer.py` (Phase 3, to be added) | Scanpy processing of the CELLxGENE dataset, malignant cell share per gene | `results/tables/singlecell_layer.csv` |
| `05_reason_cards.py` (Phase 3, to be added) | One markdown reason card per top 20 gene, generated from the result tables | `results/reason_cards/` |

The numbered scripts are working templates: the structure, inputs, outputs, and analysis choices are final, and the dataset specific lines are marked where a first run may need adjustment to the current GDC and DepMap releases. Any adjustment gets logged in `results/run_log.md`.
