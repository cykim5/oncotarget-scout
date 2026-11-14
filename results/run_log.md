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

## 2026-11-12, Phase 2 (STK11/KEAP1 subgroup, dependency, safety)

**Scripts.** `scripts/phase2_python/01_download_phase2.py`,
`02_de_shrinkage_subgroup.py`, `02b_subgroup_de.py` (companion, see
execution note 3), `03_depmap.py`, `04_gtex.py`, `05_scoring.py`,
`06_dependency_shortlist.py`.

**Data sources (all downloaded 2026-11-12).**
- TCGA-LUAD STAR counts, clinical, survival, GENCODE v36 probemap:
  same UCSC Xena GDC hub files as Phase 1, reused from
  `data/raw/xena/` (no re-download).
- Somatic mutations for STK11, KEAP1, KRAS: cBioPortal public API,
  study luad_tcga_pan_can_atlas_2018 (83 STK11, 103 KEAP1, 172 KRAS
  mutation records; sample list 566 samples). See deviation 1.
- DepMap Public 24Q4 (figshare article 27993248, DOI
  10.25452/figshare.plus.27993248.v1): CRISPRGeneEffect.csv
  (428,678,699 bytes), Model.csv (645,696), OmicsSomaticMutations.csv
  (338,945,382), AchillesCommonEssentialControls.csv (17,015). All
  size and md5 verified against the figshare-supplied values. See
  deviations 2 and 3.
- GTEx v8 median gene TPM by tissue
  (GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz,
  6,952,331 bytes), 56,200 genes x 54 tissues.

**Cohort and group counts (actual).**
- Mutation strata (522 of 528 primary tumors with mutation data; 6
  without, excluded from subgroup DE): wild type 379, STK11-only 48,
  KEAP1-only 70, both 25; combined mutant 143.
- Subgroup DE: 34,835 genes in the tumor-only matrix. Combined
  mutant vs WT: 19,922 genes at padj < 0.05 (7,222 up, 12,700 down),
  max |log2FC| 9.01. STK11-only vs WT: 13,083 significant.
  KEAP1-only vs WT: 15,947. Both vs WT: 13,299.
- DepMap: 2,105 models, 260 lung lineage; 67 STK11/KEAP1-mutant lung
  lines (19 with both mutations), 40 with CRISPR data; 86 other lung
  lines with CRISPR data. CRISPRGeneEffect matrix: 1,178 models x
  17,916 genes. Common essential controls list: 1,247 genes.
- Safety layer: 54,592 unique symbols; median safety window 0.979.
- Survival: mutant subgroup cohort 132 patients, 59 deaths. Cox
  models run for 503 genes (top 500 preliminary + 3 answer key), all
  converged; 52 with p < 0.05.
- Ranking: 28,328 genes ranked (dependency coverage 16,711, safety
  24,131, outcome 503); 6,945 excluded as unscoreable (no DepMap
  coverage and no expression signal); 848 genes with dependency > 0.

**Packages.** Python 3.12.3; pydeseq2 0.5.4 (incl. `lfc_shrink`,
apeglm Cauchy prior); lifelines 0.30.3; pandas 2.3.3; numpy 2.5.4;
scipy 1.18.1; matplotlib 3.6.3. Same project venv as Phase 1.

**Deviations from docs/project-spec.md and the Phase 2 task, all
deliberate and recorded here.**
1. Mutation calls from the cBioPortal API instead of the Xena GDC
   hub mutect2 file. The Xena file
   (TCGA-LUAD.mutect2_snv.tsv.gz) and the legacy Xena hub mutation
   files both return HTTP 403. cBioPortal serves the same underlying
   TCGA mutation calls for the PanCancer Atlas LUAD study. Only
   non-silent mutations count toward mutant status.
2. DepMap from the Broad's figshare deposit (24Q4, the newest
   release on figshare; hosting stopped after 24Q4) instead of the
   portal download API, which returns an HTML human-verification
   page to scripted requests. Release is pinned and citable.
3. First download of OmicsSomaticMutations.csv arrived truncated
   (333,381,632 of 338,945,382 bytes) and was caught by the md5
   check; re-downloaded in full and verified. Script 01 now verifies
   size and md5 on the temp file before it replaces the destination.
4. Mutant vs wild type fold changes are unshrunk (shrinkage applied
   only to the tumor vs normal comparison). Both groups are tumors,
   so the near-zero baseline artifact is far less likely, and the
   combined table confirms it (max |log2FC| 9.01, zero genes above
   10). The per-gene tables do contain a few extreme unshrunk values
   from low-count genes in small groups (PSG5 at 79.2 in KEAP1-only,
   padj 0.86; GCG at 49.2 in STK11-only, padj 0.95); they are
   non-significant, are not used in scoring, and are disclosed in
   the Phase 2 report.
5. Dependency layer rule, fixed during the run: the first version
   of script 05 set genes failing the essentiality gate to NaN,
   which wrongly excused them from the highest-weighted layer via
   renormalization (dependency coverage showed as 848 genes and
   CD274 fell out of the ranking entirely). The corrected rule,
   committed in script 05: a gene with DepMap coverage that fails
   the gate scores 0 on this layer (evidence against it); only
   genes with no DepMap coverage at all are missing. The buggy
   run's outputs were deleted and the answer key file truncated
   back to the Phase 1 section before re-running; no number from
   the buggy run survives in any committed file.
6. Outcome layer scaling weakness, disclosed not fixed: the 0 to 1
   scaling of Cox log hazard ratios is stretched by one extreme
   outlier (RN7SL564P, log HR -7.42, p = 0.076), so a gene with no
   survival effect can show an outcome value near 0.93. The
   relative ordering is preserved and the layer weight is 0.10,
   which bounds the effect; the report directs readers to the
   absolute hazard ratios and p values in cox_phase2.csv. A future
   version should clip or winsorize the scale.
7. Optional druggability flag skipped: the DGIdb monthly TSV
   download URLs now return the site's web page HTML instead of
   data files, so no clean public source was available in this run.
8. Human Cell Atlas normal lung component of the safety layer not
   done in Phase 2 (GTEx heart/brain/liver only); deferred to
   Phase 3 alongside the single cell work.

**Execution history.**
- Script 02 completed Part A (shrinkage) and then crashed at the
  start of Part B on a pandas set-indexer bug
  (`counts.loc[:, covered]` with `covered` a set; fixed in script
  02 since). Part B was completed by 02b_subgroup_de.py with
  identical analysis code.
- The first run of 02b was interrupted by a runtime restart with no
  confirmed result; verified afterward that no output files had
  been written and the process was gone, so it was re-run from
  scratch. The checkpoint commit 4cb06e8 (scripts, Part A outputs,
  dependency and safety layers, Phase 1 correction note) was made
  between the two runs.
- Phase 2 headline numbers: shrinkage moves MAGEC1 from rank 1 to
  rank 145 in the rebuilt Phase 1 ranking (new rank 1 PSG1);
  subgroup top 5 by fold change AC010941.1, PCSK2, CGA, AL353586.1,
  CALCA; dependency-supported shortlist led by INTS6, NFE2L2,
  SPRR1B, UXS1, NCAPG2 (13 of 20 robust); answer key v2 KRAS
  #4993, EGFR #6166, CD274 #5646 of 28,328.
