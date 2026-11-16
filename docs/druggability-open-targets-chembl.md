# Druggability, Open Targets and ChEMBL check

Run 2026-11-16. This closes the two gaps the Limits section of
[druggability-top5.md](druggability-top5.md) left open:

1. No Open Targets tractability score had been retrieved for any of the top 5.
2. "No known inhibitor found" meant the Nov 15 searches turned nothing up. It
   was not a database answer.

Both are now answered from the APIs rather than from reading. The script is
`scripts/druggability_checks.py` and the raw response is
`results/tables/druggability_open_targets_chembl.json`.

## What was queried

- **Open Targets Platform**, GraphQL API v4. Each gene resolved to an Ensembl
  id, then the `tractability` field, which carries 28 flags per target across
  four modalities: SM (small molecule), AB (antibody), PR (PROTAC or other
  degrader), OC (other clinical). A flag is either true or false, so "2 of 28
  true" is the whole answer for that target, not a sample.
- **ChEMBL**, REST API. Targets were looked up by UniProt accession rather than
  by gene symbol, because a symbol search for UXS1 returns nothing while the
  accession search is exact. For every ChEMBL target found, the total activity
  count and every annotated mechanism of action were pulled.

## Results

| Gene | UniProt | SM flags true | ChEMBL targets | ChEMBL activities | Annotated mechanisms |
| --- | --- | --- | --- | --- | --- |
| INTS6 | Q9UL03 | none | 1 | 6 | none |
| NFE2L2 | Q16236 | Approved Drug, High-Quality Ligand, Small Molecule Binder | 3 | 101,242 | 2 |
| UXS1 | Q8NBZ7 | Structure with Ligand, Med-Quality Pocket | 0 | 0 | none |
| NCAPG2 | Q86XI2 | none | 0 | 0 | none |
| SPRR1B | P22528 | none | 0 | 0 | none |

Full flag list per gene, true flags only:

| Gene | Flags returning true |
| --- | --- |
| INTS6 | PR Database Ubiquitination, PR Half-life Data |
| NFE2L2 | SM Approved Drug, SM High-Quality Ligand, SM Small Molecule Binder, AB GO CC high conf, PR UniProt Ubiquitination, PR Database Ubiquitination |
| UXS1 | SM Structure with Ligand, SM Med-Quality Pocket, AB UniProt SigP or TMHMM, PR Half-life Data |
| NCAPG2 | PR Database Ubiquitination |
| SPRR1B | AB GO CC high conf |

## What this confirms

**The reasoned verdicts hold.** Every verdict in druggability-top5.md survives
the check, and now rests on a retrieved field rather than on protein class
reasoning.

- **UXS1** is the only gene besides NFE2L2 with any small molecule flag true,
  and its two flags are exactly the ones that matter for a screening start:
  a structure with a ligand bound, and a pocket Open Targets scores as medium
  quality. ChEMBL has no target entry for it at all, so "no known inhibitor" is
  now a database answer, not a failed search. An enzyme with a scored pocket and
  zero compounds is a clean starting point, which is what the docking and the
  100-drug screen were built on.
- **INTS6, NCAPG2 and SPRR1B** have no small molecule flag true between them.
  Their only true flags are ubiquitination and half-life data, which are the
  degrader flags, and one antibody localization flag. This is the structured
  version of what the table already said: a degrader or knockdown route, not an
  active site inhibitor.
- **NCAPG2 and SPRR1B** have no ChEMBL target entry, so there is no compound
  history for either.

## What this changes, and it matters

**The approved NRF2 drug goes the wrong way.**

NFE2L2 carries the SM "Approved Drug" tractability flag, which at a glance looks
like the strongest result in the table. The two annotated mechanisms in ChEMBL
are what that flag rests on:

| Compound | ChEMBL | Action | Phase | First approval |
| --- | --- | --- | --- | --- |
| Omaveloxolone | CHEMBL4303525 | NRF2 activator | 4, approved | 2023 |
| Bardoxolone methyl | CHEMBL1762621 | KEAP1/NRF2 inhibitor, which releases NRF2 | 3 | not approved |

Both push NRF2 activity **up**. Omaveloxolone is approved for Friedreich ataxia,
where raising NRF2 is the point.

In STK11/KEAP1-mutant lung adenocarcinoma the goal is the opposite. KEAP1 loss
already leaves NRF2 constitutively active, and the hypothesis is that these
tumors depend on it. The direction needed is inhibition.

So the "Approved Drug" flag on NFE2L2 should not be read as a repurposing
opportunity for this subgroup. ML385 remains the only inhibitor in hand, it is
still a research tool with no clinical phase, and the validation plan's pairing
of knockout with ML385 is still the right design. If anything this sharpens the
case for UXS1, which is the one target where a scored pocket and an empty
compound landscape point the same way.

## Limits of this check

- Tractability flags describe a protein's general amenability to a modality.
  They say nothing about whether modulating it helps in this subgroup. That is
  what the CRISPR validation is for.
- A ChEMBL activity count is a measure of how much a target has been studied,
  not of how good a target it is. NFE2L2's 101,242 activities reflect decades of
  NRF2 chemistry, much of it screening against activators.
- "Med-Quality Pocket" for UXS1 is Open Targets' own grading. It is a reason to
  screen, not a prediction that a selective inhibitor exists.
- Queried 2026-11-16. Both databases are updated on a release cycle, so these
  numbers carry a date.
