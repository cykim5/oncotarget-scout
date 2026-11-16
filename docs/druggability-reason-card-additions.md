# Druggability lines for the top 5 reason cards

Prepared 2026-11-15, from docs/druggability-top5.md. Not applied yet.
Claire's open choice from chat 2026-11-15: add this as a column to the
shortlist and reason cards, or keep druggability as the standalone doc.
These are the exact lines to paste if she picks the column option.
No em dashes used. Nothing here changes the Phase 3 ranking.

## One line per reason card, under "## Watch out" or as a new "## Druggability" section

- INTS6 (rank 1): Druggability: hard for a classic small molecule inhibitor, no active site of its own. The route would be a degrader or a protein-protein interaction disruptor, not an active-site inhibitor.
- NFE2L2 (rank 2): Druggability: the only top 5 gene with a compound in hand. ML385, a selective NRF2 inhibitor (IC50 1.9 uM), works in KEAP1-deficient lung lines including A549. It is a research probe, not an approved drug, and it is already the orthogonal arm in docs/validation-plan.md.
- UXS1 (rank 3): Druggability: most actionable of the top 5 on structure. Enzyme with a defined active site (catalytic triad Thr118, Tyr147, Lys151 (construct numbering; Thr202, Tyr231, Lys235 in 4GLL)) and crystal structures PDB 4GLL and 2B69, one bound with NAD+ and UDP. No inhibitor found yet, so this is a screening starting point, not a ready drug. Independent support: Gebru et al., Cancer Res 2025, PMID 40966319, found UXS1 is a selective dependency in KEAP1-mutant lung cancer.
- NCAPG2 (rank 4): Druggability: hard for a small molecule. Structural condensin II subunit, no enzyme active site and no small-molecule pocket reported. A degrader or RNA-based knockdown would be the route, not an inhibitor.
- SPRR1B (rank 5): Druggability: hardest of the five as a drug target. Structural envelope protein with no pocket. It reads more as a biomarker than as something to inhibit with a small molecule. It also has no single-cell value in the Phase 3 ranking (4 of 5 layers only).

## Shortlist CSV column (results/tables/shortlist_dependency_supported_phase3.csv)

Add a `druggability` column for the top 5 rows only, values:
INTS6 = hard, degrader or PPI disruptor; NFE2L2 = probe available, ML385;
UXS1 = most actionable, active site solved, no inhibitor yet;
NCAPG2 = hard, degrader or RNA knockdown; SPRR1B = hardest, biomarker.
Leave the other 15 rows blank, druggability was only checked for the top 5.

## Full-20 note (for any writeup, not a reason card)

IMPDH2 is rank 20 of 20, not top 5. It is the most druggable gene in the
full 20: mycophenolic acid (CellCept / Myfortic) is a potent IMPDH
inhibitor, approved since 1995. State the rank correctly if it is mentioned.

## Limits, carried over from docs/druggability-top5.md

- No Open Targets tractability score was retrieved for any of the five.
- "No known inhibitor found" means none turned up in the Nov 15 searches,
  not proof none exists. A per-gene ChEMBL target search would settle it.
- Druggability does not change the Phase 3 ranking. It is a separate
  column for deciding what to do after CRISPR validation.
