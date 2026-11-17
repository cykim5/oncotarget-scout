"""How a tool result becomes a judgment.

RuleBackend applies thresholds that are written down here and nowhere
else. LLMBackend sends the same tool result to a model and asks for the
same decision. Both return (is_problem, severity, sentence), so the
agents and the orchestrator do not care which one is in use.

The measured experiment in this repository used RuleBackend. The model
path is implemented and needs ANTHROPIC_API_KEY to run; until it has
been run, no claim in this repository rests on it.
"""
import json
import os
import urllib.request

P_DEPLETED = 0.01        # a layer counts as depleted below this
ANSWER_KEY_PCTILE = 95.0  # a known target should sit at least this high


class RuleBackend:
    name = "rules"

    def judge_depletion(self, record):
        if record["p_depleted"] < P_DEPLETED:
            return (True, "critical",
                    f"the top {record['k']} is missing the "
                    f"{record['layer']} layer far more often than the "
                    f"ranked population is")
        return (False, "info", f"{record['layer']} is not depleted at the top")

    def judge_answer_key(self, record):
        if record["rank"] is None:
            return (True, "warning",
                    f"{record['gene']} is a known target and is not ranked")
        if record["percentile"] < ANSWER_KEY_PCTILE:
            return (True, "critical",
                    f"{record['gene']} is an approved target in this disease "
                    f"and sits at the {record['percentile']:.1f}th percentile")
        return (False, "info",
                f"{record['gene']} is at the {record['percentile']:.1f}th "
                f"percentile")

    def judge_card(self, record):
        return (True, "warning",
                f"{record['gene']}: the card says {record['field']} is "
                f"{record['card']} and the table says {record['table']}")

    def judge_shortlist(self, record):
        if len(record["unsupported"]) > 0:
            return (True, "critical",
                    f"{len(record['unsupported'])} shortlisted genes have no "
                    f"dependency evidence")
        return (False, "info", "every shortlisted gene has dependency evidence")


class LLMBackend:
    """The same decisions, made by a model instead of a threshold.

    Implemented and not yet run. It exists so the comparison between a
    rule and a model is a configuration change rather than a rewrite.
    """
    name = "llm"

    def __init__(self, model="claude-opus-5", api_key=None):
        self.model = model
        self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
        if self.api_key is None:
            raise RuntimeError(
                "LLMBackend needs ANTHROPIC_API_KEY. The measured results "
                "in docs/verification-layer.md used RuleBackend.")

    def _ask(self, question, record):
        prompt = (
            "You are checking the output of a scientific ranking pipeline.\n"
            f"{question}\n\n"
            f"Tool result as JSON:\n{json.dumps(record, indent=2)}\n\n"
            "Answer with JSON only: "
            '{"is_problem": true or false, '
            '"severity": "info" or "warning" or "critical", '
            '"sentence": "one sentence a scientist would write"}')
        body = json.dumps({
            "model": self.model,
            "max_tokens": 400,
            "messages": [{"role": "user", "content": prompt}],
        }).encode()
        req = urllib.request.Request(
            "https://api.anthropic.com/v1/messages", data=body,
            headers={"content-type": "application/json",
                     "x-api-key": self.api_key,
                     "anthropic-version": "2023-06-01"})
        with urllib.request.urlopen(req, timeout=60) as r:
            payload = json.loads(r.read())
        text = payload["content"][0]["text"].strip()
        if text.startswith("```"):
            text = text.split("\n", 1)[1].rsplit("```", 1)[0]
        out = json.loads(text)
        return (bool(out["is_problem"]), out["severity"], out["sentence"])

    def judge_depletion(self, record):
        return self._ask("Is the top of this ranking missing a layer of "
                         "evidence that the rest of the list has?", record)

    def judge_answer_key(self, record):
        return self._ask("This gene is an approved drug target in this "
                         "disease. Is its position a problem?", record)

    def judge_card(self, record):
        return self._ask("A written summary disagrees with the source "
                         "table. How serious is it?", record)

    def judge_shortlist(self, record):
        return self._ask("Does every shortlisted gene carry the evidence "
                         "the shortlist claims to require?", record)
