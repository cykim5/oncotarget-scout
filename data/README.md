# Data sources and download steps

No raw or processed data files are committed to this repository. This file is the record of where every input comes from, so the analysis can be rebuilt.

## 1. TCGA lung adenocarcinoma (TCGA-LUAD)

- Source: Genomic Data Commons (GDC), accessed with the R package TCGAbiolinks.
- Data: STAR counts RNA-seq (tumor and solid tissue normal), masked somatic mutations (for STK11, KEAP1, KRAS status), and clinical data (stage, age, survival).
- Download: run `scripts/01_download_tcga.R`. Processed tables are written to `data/processed/` (git ignored).
- Cohort note: about 500 tumors and about 59 normal samples in recent GDC releases. Exact counts are printed by the download script and recorded in the run log, rather than hard coded here.

## 2. DepMap CRISPR dependency

- Source: DepMap portal data downloads (depmap.org/portal/data_page), a current quarterly public release. Record the release name (for example 25Q2) in the run log when downloaded; results depend on it.
- Files: `CRISPRGeneEffect.csv` (Chronos gene effect; more negative means more essential; around minus 1 is the common essential range, around minus 0.5 is a likely dependency), `Model.csv` (cell line annotations), and the somatic mutation table for STK11/KEAP1/KRAS status per line.
- Download is manual through the portal page; place the files under `data/raw/depmap/` (git ignored).
- Used in this project: Public 24Q4 (figshare article 27993248) for the main dependency layer, and Public 20Q2 (figshare article 12280541, CERES scores) for the Phase 3 time split check. Both were downloaded from the Broad's figshare deposit because the portal download API returns a human-verification page to scripts; see the run log.

## 3. Single cell, CZ CELLxGENE Discover

- Source: cellxgene.cziscience.com collections; choose one human lung adenocarcinoma dataset with author provided cell type annotations and malignant cell labels.
- Download the H5AD file into `data/raw/cellxgene/` (git ignored) and record the dataset title and version in the run log.
- Used in this project: dataset 4ca422ae-e4ad-4884-9073-7c3799cd9143 (Chromium 10x, LUAD and squamous cell carcinoma, 4 donors, bioRxiv 2024.11.01.621259). Only the LUAD cells (2,558 of 35,954) were used for the single cell layer; see the run log.

## 4. Normal tissue, GTEx

- Source: GTEx portal median gene expression (gtexportal.org), or the recount3 processed GTEx data from R.
- Data: median TPM per tissue for all scored genes; heart, brain, liver, and lung are the tissues used in the safety window.

## 5. Pathway gene sets

- Source: MSigDB Hallmark collection and Reactome, accessed from R via the msigdbr package, or downloaded GMT files recorded in the run log.

## Run log

Every analysis run appends one entry to `results/run_log.md`: date, script versions, dataset releases, cohort counts, and any deviation from this file.
