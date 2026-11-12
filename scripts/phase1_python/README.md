# Phase 1, Python route

`run_phase1.py` is the executed Phase 1 pipeline. It covers the same Phase 1
deliverables as the R templates in `scripts/` (TCGA-LUAD download, DESeq2
differential expression, GSEA, survival, partial ranking, answer key check)
using Python equivalents: UCSC Xena GDC-hub STAR counts instead of a
TCGAbiolinks GDC download, pydeseq2 instead of DESeq2 in R, gseapy for GSEA,
and lifelines for the Cox and Kaplan-Meier analyses. Every substitution is
recorded in `results/run_log.md`. The R templates remain in place, unmodified,
as the reference implementation of the original design.

Data downloads (public, no credentials), from the repository root:

```
mkdir -p data/raw/xena
for f in TCGA-LUAD.star_counts.tsv.gz TCGA-LUAD.clinical.tsv.gz \
         TCGA-LUAD.survival.tsv.gz gencode.v36.annotation.gtf.gene.probemap; do
  curl -L -o data/raw/xena/$f https://gdc-hub.s3.us-east-1.amazonaws.com/download/$f
done
```

Note: the Xena STAR counts file stores log2(count + 1) values. The script
recovers integer counts as round(2**x - 1) before running pydeseq2, and says
so in its own docstring and in the run log.

Run:

```
.venv/bin/python scripts/phase1_python/run_phase1.py
```
