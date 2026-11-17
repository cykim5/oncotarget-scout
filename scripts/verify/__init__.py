"""A verification layer for the pipeline's own output.

The pipeline produces a ranking. This package checks it, by splitting
the job into specialists that run in parallel and then cross checking
what they found against each other.

It exists because a real defect in this project survived three phases
of review: the score renormalized its weights over whichever layers a
gene had, so genes with missing data outranked genes with complete
evidence. The question this package answers is whether an automated
layer would have caught it, and whether splitting the work up catches
more than one pass over everything.
"""
from .backends import LLMBackend, RuleBackend
from .context import Context
from .findings import Finding
from .orchestrator import run_decomposed, run_single_pass, summarize

__all__ = ["Context", "Finding", "RuleBackend", "LLMBackend",
           "run_decomposed", "run_single_pass", "summarize"]
