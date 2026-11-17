"""The tools the agents call.

Each tool returns plain numbers, never a judgment. The judgment is the
backend's job. Keeping them apart is what lets the same agents run on a
rule backend and on a model backend and be compared.
"""
from scipy.stats import binomtest

from .context import ANSWER_KEY, LAYERS, numbers_in


def layer_depletion(ctx, k=20):
    """Is the top of the list missing a layer the population has?

    For each layer, compare how often it is present in the top k with
    how often it is present across the whole ranked list, and give the
    one sided binomial p value for the top being depleted.

    This is a general check. It does not know what any layer means.
    """
    out = []
    top = ctx.top(k)
    for layer in LAYERS:
        rate = ctx.layer_presence_rate(layer)
        present = int(top[layer].notna().sum())
        if rate == 0.0:
            p = 1.0
        else:
            p = float(binomtest(present, k, rate, alternative="less").pvalue)
        out.append({"layer": layer, "present_in_top": present, "k": k,
                    "population_rate": rate, "p_depleted": p})
    return out


def answer_key_position(ctx):
    """Where the known drug targets land in the ranking."""
    n = len(ctx.ranking)
    out = []
    for gene in ANSWER_KEY:
        row = ctx.gene_row(gene)
        if row is None:
            out.append({"gene": gene, "rank": None, "n": n, "percentile": None})
            continue
        rank = int(row["rank"])
        out.append({"gene": gene, "rank": rank, "n": n,
                    "percentile": 100.0 * (1.0 - (rank - 1) / n)})
    return out


def card_numbers_match(ctx, fields):
    """Re-read each reason card's numbers from the ranking table.

    fields maps a label printed in the card to the column it came from.
    Returns one record per mismatch, and the number of cards checked.
    """
    mismatches = []
    checked = 0
    for gene, text in ctx.cards():
        row = ctx.gene_row(gene)
        if row is None:
            mismatches.append({"gene": gene, "field": "(whole card)",
                               "card": None, "table": None,
                               "note": "gene is not in the ranking at all"})
            continue
        checked = checked + 1
        for label, column in fields.items():
            said = numbers_in(text, label)
            if said is None:
                continue
            if column not in row.index:
                continue
            actual = row[column]
            if actual != actual:        # NaN
                continue
            if abs(said - float(actual)) > 0.0015:
                mismatches.append({"gene": gene, "field": label,
                                   "card": said, "table": float(actual),
                                   "note": "card does not match the table"})
    return {"checked": checked, "mismatches": mismatches}


def shortlist_dependency_support(ctx):
    """Does every shortlisted gene actually carry the dependency layer?"""
    if ctx.shortlist is None:
        return {"checked": 0, "unsupported": []}
    bad = ctx.shortlist[ctx.shortlist["dependency"].isna() == True]
    return {"checked": int(len(ctx.shortlist)),
            "unsupported": list(bad["gene"])}
