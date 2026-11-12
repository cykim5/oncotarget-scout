# Run log

## 2026-11-11, Phase 1 v1 (full LUAD cohort)

**Scripts.** `scripts/phase1_python/run_phase1.py` (stages 1 to 3 as executed,
plus the stages 4 to 6 definitions) and
`scripts/phase1_python/run_phase1_stages4to6.py` (stages 4 to 6 as executed).
See the execution note at the end of this entry for why two scripts were used.

**Data source.** UCSC Xena GDC hub (https://gdc-hub.s3.us-east-1.amazonaws.com/download/),
downloaded 2026-11-11:
- `TCGA-LUAD.star_counts.tsv.gz` (66,495,260 bytes): GDC STAR counts for
  TCGA-LUAD as served by Xena. The file stores log2(count + 1) values.
- `TCGA-LUAD.clinical.tsv.gz` (149,532 bytes)
- `TCGA-LUAD.survival.tsv.gz` (5,253 bytes)
- `gencode.v36.annotation.gtf.gene.probemap` (3,206,945 bytes) for
  Ensembl to gene symbol mapping (GENCODE v36).
Xena does not stamp a GDC release number inside these files; the download
date above is the version record for this run.

**Cohort counts (actual).**
- Downloaded matrix: 60,660 Ensembl rows x 589 samples.
- Differential expression: 528 primary tumor samples (barcode type 01) and
  59 solid tissue normal samples (type 11). 2 recurrent tumor samples
  (type 02) were excluded.
- After symbol mapping (duplicate symbols summed): 59,427 unique symbols.
  After the prefilter (count >= 10 in at least 10 samples): 35,273 genes
  carried into DESeq2 and the ranking.
- Survival cohort: 485 patients (one primary tumor sample per patient,
  complete pathologic stage, age, and overall survival), 177 deaths.

**Packages.** Python 3.12.3; pydeseq2 0.5.4; gseapy 1.3.2 (prerank,
MSigDB_Hallmark_2020 library fetched from Enrichr, 1,000 permutations,
seed 42); lifelines 0.30.3; pandas 2.3.3; numpy 2.5.4; scipy 1.18.1;
matplotlib 3.6.3. Installed into the project venv (`.venv/`).

**Deviations from docs/project-spec.md, all deliberate and recorded here.**
1. Python route instead of R. The R / Bioconductor stack (TCGAbiolinks,
   DESeq2) is not installed in this environment. pydeseq2 is the Python
   port of DESeq2, so the differential expression method is the designed
   one; the code path differs.
2. Data access via the UCSC Xena GDC hub instead of a TCGAbiolinks GDC
   download. Same underlying GDC STAR counts. Because Xena publishes
   log2(count + 1), integer counts were recovered as round(2**x - 1)
   before running pydeseq2. The transform is invertible at the precision
   Xena publishes; any residual rounding sits in the lowest count digits.
3. GSEA was run on the Hallmark collection only in v1. The spec Phase 1
   text also names Reactome; it is deferred, not dropped.
4. No PCA figure was produced. The spec lists a PCA figure next to the
   volcano plot; v1 ships the volcano, the top 20 bar chart, and the
   Kaplan-Meier figure.
5. The survival (outcome) layer was computed for the top 300 genes by
   tumor specificity plus the 3 answer key genes, 303 Cox models total,
   following the spec wording "survival association for top candidates."
   All other genes carry coverage 1 of 5 in the partial table, the 303
   carry coverage 2 of 5, and scores are renormalized over covered layers
   as docs/target-score.md specifies. Rank comparisons across coverage
   groups are therefore not like for like, and the report says so.
6. Tumor specificity uses unshrunken DESeq2 log fold changes, as designed
   (padj as a gate, fold change as the ranking input). Consequence: the
   scale maximum is set by MAGEC1 (log2 fold change 43.3), a gene that is
   nearly silent in normal lung, which compresses every other gene's
   specificity value. No post hoc adjustment was made; the effect is
   analyzed openly in the Phase 1 report.
7. Execution history. The first run of run_phase1.py crashed after the
   DESeq2 fit on a results parsing bug (gene names taken from the wrong
   column; fixed). The second run wrote the DE and GSEA tables and was
   then terminated by this environment before stages 4 to 6 finished.
   Stages 4 to 6 were completed by run_phase1_stages4to6.py reading the
   two saved stage outputs, with identical analysis code, so every number
   traces to a committed script and the saved intermediate tables.
