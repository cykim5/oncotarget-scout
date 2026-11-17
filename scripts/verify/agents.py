"""The specialists.

Each one has a narrow remit, its own tools, and no knowledge of what
the others are doing. That separation is what gives the cross check
agent something independent to compare.

Each agent works in two phases. collect() calls its tools and returns
the raw records, which is cheap and local. judge() turns one record
into a Finding by asking the backend, which is the expensive part when
the backend is a model. Keeping them apart is what lets the
orchestrator parallelize at whichever level it wants.
"""
from . import tools
from .findings import Finding

CARD_FIELDS = {
    "Tumor specificity layer value": "tumor_specificity",
    "Dependency layer value": "dependency",
    "Safety window layer value": "safety_window",
    "Malignant cell share layer value": "malignant_share",
    "Selective essentiality": "selective_essentiality",
    "Mean Chronos in mutant lung lines": "mean_chronos_mutant_lung",
}


class Agent:
    name = "agent"

    def __init__(self, backend):
        self.backend = backend

    def collect(self, ctx):
        """Return a list of (kind, record) for judge() to decide on."""
        raise NotImplementedError

    def judge(self, ctx, kind, record):
        """Turn one record into a Finding, or None when it is fine."""
        raise NotImplementedError

    def run(self, ctx):
        out = []
        for kind, record in self.collect(ctx):
            finding = self.judge(ctx, kind, record)
            if finding is not None:
                out.append(finding)
        return out


class LayerBalanceAgent(Agent):
    """Is the top of the ranking evidenced like the rest of the list?

    This agent does not know what any layer means. It only asks whether
    the genes at the top are systematically missing something the
    population has, which would mean the score rewards absence.
    """
    name = "layer-balance"

    def __init__(self, backend, k=20):
        super().__init__(backend)
        self.k = k

    def collect(self, ctx):
        return [("depletion", r) for r in tools.layer_depletion(ctx, self.k)]

    def judge(self, ctx, kind, record):
        is_problem, severity, sentence = self.backend.judge_depletion(record)
        if is_problem == False:
            return None
        return Finding(
            agent=self.name, severity=severity, claim=sentence,
            evidence=(f"{record['present_in_top']} of the top "
                      f"{record['k']} carry {record['layer']}, against "
                      f"{record['population_rate'] * 100:.1f} percent of "
                      f"the ranked population, p={record['p_depleted']:.2e}"),
            where=ctx.ranking_file, tools_used=["layer_depletion"])


class ControlAgent(Agent):
    """Does the ranking still pass the checks fixed before it ran?"""
    name = "controls"

    def collect(self, ctx):
        out = [("answer_key", r) for r in tools.answer_key_position(ctx)]
        record = tools.shortlist_dependency_support(ctx)
        if record["checked"] > 0:
            out.append(("shortlist", record))
        return out

    def judge(self, ctx, kind, record):
        if kind == "answer_key":
            is_problem, severity, sentence = self.backend.judge_answer_key(record)
            if is_problem == False:
                return None
            rank = record["rank"]
            return Finding(
                agent=self.name, severity=severity, claim=sentence,
                evidence=(f"rank {rank} of {record['n']}" if rank is not None
                          else "not present in the ranked list"),
                where=ctx.ranking_file, tools_used=["answer_key_position"])
        is_problem, severity, sentence = self.backend.judge_shortlist(record)
        if is_problem == False:
            return None
        return Finding(
            agent=self.name, severity=severity, claim=sentence,
            evidence=f"genes without the layer: {record['unsupported']}",
            where="shortlist", tools_used=["shortlist_dependency_support"])


class CardFidelityAgent(Agent):
    """Do the written reason cards match the table they came from?"""
    name = "card-fidelity"

    def collect(self, ctx):
        record = tools.card_numbers_match(ctx, CARD_FIELDS)
        return [("card", m) for m in record["mismatches"]]

    def judge(self, ctx, kind, record):
        is_problem, severity, sentence = self.backend.judge_card(record)
        if is_problem == False:
            return None
        return Finding(
            agent=self.name, severity=severity, claim=sentence,
            evidence=record["note"],
            where=f"reason card {record['gene']}",
            tools_used=["card_numbers_match"])


class CrossCheckAgent:
    """Reads the other agents' findings and looks for one root cause.

    This agent cannot exist without decomposition. A single pass that
    produces one flat list has nothing independent to cross check, so
    this is the part of the architecture that has to earn its keep.
    """
    name = "cross-check"

    def __init__(self, backend):
        self.backend = backend

    def run(self, ctx, findings):
        depleted = [f for f in findings if f.agent == "layer-balance"]
        controls = [f for f in findings if f.agent == "controls"]
        if len(depleted) == 0 or len(controls) == 0:
            return []
        layers = sorted({f.evidence.split("carry ")[1].split(",")[0]
                         for f in depleted if "carry " in f.evidence})
        genes = ", ".join(sorted({f.claim.split()[0] for f in controls}))
        return [Finding(
            agent=self.name, severity="critical",
            claim=("two independent checks fail together, which points at one "
                   "cause in how the score is built rather than at the data"),
            evidence=(f"the top of the list is depleted of {', '.join(layers)}, "
                      f"and the known targets {genes} rank far below where "
                      f"established biology puts them; a score that lets a "
                      f"gene rank without a layer will rank genes that lack it "
                      f"above genes that have it"),
            where=ctx.ranking_file,
            tools_used=["layer_depletion", "answer_key_position"])]


def specialists(backend, k=20):
    return [LayerBalanceAgent(backend, k), ControlAgent(backend),
            CardFidelityAgent(backend)]
