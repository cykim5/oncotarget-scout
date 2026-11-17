"""What a verification agent returns."""
from dataclasses import dataclass, field, asdict
from typing import List

SEVERITIES = ["info", "warning", "critical"]


@dataclass
class Finding:
    agent: str
    severity: str          # one of SEVERITIES
    claim: str             # what is wrong, in one sentence
    evidence: str          # the numbers behind the claim
    where: str             # file or table the claim is about
    tools_used: List[str] = field(default_factory=list)

    def __post_init__(self):
        if self.severity not in SEVERITIES:
            raise ValueError(f"severity must be one of {SEVERITIES}")

    def key(self):
        """Identity of a finding, so two runs can be compared."""
        return (self.agent, self.claim)

    def to_dict(self):
        return asdict(self)


def worst(findings):
    """The highest severity in a list, or 'info' when it is empty."""
    rank = {"info": 0, "warning": 1, "critical": 2}
    if len(findings) == 0:
        return "info"
    return max(findings, key=lambda f: rank[f.severity]).severity
