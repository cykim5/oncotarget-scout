"""Everything the agents are allowed to read.

The agents never reach into the filesystem themselves. They get a
Context, which loads the committed tables once and hands out copies.
That is what makes a run reproducible and what makes it possible to
point the same agents at the pre-fix and post-fix rankings.
"""
import re
from pathlib import Path

import pandas as pd

LAYERS = ["tumor_specificity", "dependency", "safety_window",
          "malignant_share", "outcome_link"]
ANSWER_KEY = ["EGFR", "KRAS", "CD274"]


class Context:
    def __init__(self, root, ranking_file, cards_dir=None, shortlist_file=None):
        self.root = Path(root)
        self.ranking_file = str(ranking_file)
        self.ranking = pd.read_csv(self.root / ranking_file)
        self.cards_dir = self.root / cards_dir if cards_dir else None
        self.shortlist = None
        if shortlist_file is not None:
            self.shortlist = pd.read_csv(self.root / shortlist_file)

    def top(self, k):
        return self.ranking.head(k).copy()

    def layer_presence_rate(self, layer):
        return float(self.ranking[layer].notna().mean())

    def gene_row(self, gene):
        hit = self.ranking[self.ranking["gene"] == gene]
        if len(hit) == 0:
            return None
        return hit.iloc[0]

    def cards(self):
        """Every reason card as (gene, text). Empty when none are present."""
        if self.cards_dir is None or self.cards_dir.exists() == False:
            return []
        out = []
        for p in sorted(self.cards_dir.glob("*.md")):
            out.append((p.stem, p.read_text()))
        return out


def numbers_in(text, label):
    """Pull the number that follows a label in a reason card.

    Returns None when the label is not present, so a missing line is
    never confused with a wrong value.
    """
    m = re.search(re.escape(label) + r":?\s*(-?\d+\.?\d*)", text)
    if m is None:
        return None
    return float(m.group(1))
