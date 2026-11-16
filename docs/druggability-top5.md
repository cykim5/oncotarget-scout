# Druggability check, top 5 shortlist

Prepared 2026-11-15, after the Phase 3 ranking. This answers the question asked in chat that day: how do you check druggability, and what does it say about the top 5.

Top 5, from results/tables/shortlist_dependency_supported_phase3.csv: 1 INTS6, 2 NFE2L2, 3 UXS1, 4 NCAPG2, 5 SPRR1B.

Method, the three checks from chat: (1) does a known compound bind the protein, and is any of them an inhibitor, (2) does the protein have a pocket a drug could fit into, shown by a crystal structure with a small molecule bound or a defined active site, (3) what is the track record of the protein family, enzymes with active sites are usually druggable, transcription factors and structural proteins are usually hard.

## Table

| Rank | Gene | Protein class | Known inhibitor | Pocket / structure | Verdict |
| --- | --- | --- | --- | --- | --- |
| 1 | INTS6 | Scaffold subunit of the Integrator complex, positions PP2A | None found | No active site of its own. INTS6 blocks the PP2A active site in the closed Integrator-PP2A conformation, it does not have a pocket to block. | Hard for a classic small molecule inhibitor. A degrader (PROTAC) or a protein-protein interaction disruptor would be the route, not an active-site inhibitor. |
| 2 | NFE2L2 | Transcription factor, encodes NRF2 | Yes, ML385, a selective NRF2 inhibitor, IC50 1.9 uM | ML385 binds the Neh1 / CNC-bZIP domain and stops the NRF2-MAFG complex binding DNA. Works in KEAP1-deficient NSCLC lines (H460, A549, H1299), including A549, one of the validation test lines. | Druggable as a probe, with a caveat. ML385 is a research tool, not an approved drug, and transcription factors are still a hard class. Already named in docs/validation-plan.md as the orthogonal arm for NFE2L2. |
| 3 | UXS1 | Enzyme, UDP-xylose synthase, a lyase / decarboxylase in the short-chain dehydrogenase/reductase family | None found | Yes. Crystal structures exist: PDB 4GLL (2.5 A) and 2B69, plus a 1.26 A structure bound with NAD+ and UDP. Defined active site with catalytic triad Thr118, Tyr147, Lys151 (construct numbering, as in 2B69; Thr202, Tyr231, Lys235 in 4GLL). | Most druggable of the top 5 on structure. No compound yet, so this is a screening starting point, not a ready drug. Independent support: a Dec 2025 Cancer Research paper (Gebru et al., PMID 40966319) found UXS1 is a selective dependency in KEAP1-mutant lung cancer and that losing it sensitizes cells to clinical checkpoint inhibitors. |
| 4 | NCAPG2 | Structural subunit of the condensin II complex, chromosome condensation and segregation in mitosis | None found | No enzyme active site and no small-molecule pocket reported. Evidence so far is knockdown only. | Hard for a small molecule. A degrader or an RNA-based knockdown would be the route, not an inhibitor. |
| 5 | SPRR1B | Structural envelope protein, cornifin-B, crosslinked in keratinocytes | None found | No pocket. It is a proline-rich structural protein that transglutaminase crosslinks into the cell envelope. | Hardest of the five as a drug target. It reads more as a biomarker than as something to inhibit with a small molecule. Note it also has no single-cell value in the Phase 3 ranking (4 of 5 layers only). |

## What this changes

The ranking order and the druggability order are different, and that is the point of the check.

- UXS1 (rank 3) moves up as the most actionable target: an enzyme with a solved active site and independent 2025 validation in exactly this KEAP1-mutant setting.
- NFE2L2 (rank 2) is the only top 5 gene with a compound already in hand, ML385, which is why the validation plan already pairs knockout with ML385.
- INTS6 (rank 1) stays first on the biology, but it is the hardest to drug directly. If it confirms in the lab, the next step is a degrader or PPI approach, not a standard inhibitor screen.
- NCAPG2 and SPRR1B (ranks 4 and 5) are the weakest on druggability. They are still worth the CRISPR test, because the test is about dependency, not about having a drug today, but they should not be presented as near-term drug targets.

## Correction to an earlier chat claim

In chat on Nov 15, IMPDH2 was described as one of the top 5 with an approved inhibitor. That was wrong on the rank. IMPDH2 is rank 20 of 20 in the Phase 3 shortlist, not top 5. The inhibitor claim itself is right: mycophenolic acid (CellCept / Myfortic) is a potent, reversible IMPDH inhibitor, approved since 1995 for transplant immunosuppression. So IMPDH2 is the most druggable gene in the full 20, just not in the top 5. It is worth one line in any writeup, with the rank stated correctly.

## Limits of this check

- No Open Targets tractability score was retrieved for any of the five. Web searches on Nov 15 returned papers and GeneCards pages, not the Open Targets target pages, so the verdicts above rest on protein class, published structures, and known compounds, not on a scored tractability field. Pulling the Open Targets scores would be the next verification step before this table goes into a preprint or application.
- "No known inhibitor found" means none turned up in the searches run Nov 15 (ChEMBL / Open Targets / literature). It is not proof that none exists. A ChEMBL target search per gene is the check that would settle it.
- Druggability does not change the Phase 3 ranking. It is a separate column to add to the shortlist, for deciding what to do after the CRISPR validation, not a reason to reorder the ranking itself.

## Sources checked Nov 15

- Shortlist and ranks: results/tables/shortlist_dependency_supported_phase3.csv
- Validation plan, including the ML385 arm: docs/validation-plan.md
- NFE2L2 / ML385: Singh et al. 2016 via InvivoChem and MDPI reviews, IC50 1.9 uM, KEAP1-deficient NSCLC activity
- UXS1 structure: PDB 4GLL, PDB 2B69, J Biol Chem structure and mechanism paper (1.26 A, NAD+ and UDP bound, Tyr147 catalytic base; Tyr231 in 4GLL)
- UXS1 dependency: Gebru et al., Cancer Res 2025 Dec 1, PMID 40966319
- INTS6 / Integrator-PP2A: Nature 2024 structural basis paper, INTS6 blocks the PP2A active site in the closed conformation; Cell Rep 2023, INTS5/8 plus INTS6 recruit PP2A
- NCAPG2: Front Immunol 2023 pan-cancer analysis, condensin II subunit, knockdown evidence only
- SPRR1B: GeneCards / NCBI, cornifin-B envelope protein, transglutaminase crosslinking, no pocket
- IMPDH2 / mycophenolic acid: FDA DailyMed label, IMPDH inhibitor, approved 1995
