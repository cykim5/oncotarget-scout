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

    Nothing in this repository has been measured with this. Every
    number in docs/verification-layer.md came from RuleBackend. This is
    here so that running it is a flag and an API key, not a rewrite.

    Two providers, chosen by whichever key is set:

      ANTHROPIC_API_KEY  the Claude Messages API
      NEBIUS_API_KEY     Nebius, the OpenAI compatible endpoint that
                         Gene-ius already uses, so leftover credits
                         from that project work here

    Set VERIFY_MODEL to pick a model. The default is Claude Opus 5.5.
    claude-haiku-5-5 costs roughly a cent for a whole run and is the
    more interesting test anyway: if the cheapest model matches
    thresholds written by hand, that is the result worth reporting.
    """
    name = "llm"

    ANTHROPIC = "https://api.anthropic.com/v1/messages"
    NEBIUS = "https://api.studio.nebius.ai/v1/chat/completions"

    SCHEMA = {
        "type": "object",
        "properties": {
            "is_problem": {"type": "boolean"},
            "severity": {"type": "string",
                         "enum": ["info", "warning", "critical"]},
            "sentence": {"type": "string"},
        },
        "required": ["is_problem", "severity", "sentence"],
        "additionalProperties": False,
    }

    def __init__(self, model=None, api_key=None, provider=None):
        self.provider = provider or self._detect_provider()
        if self.provider == "anthropic":
            self.api_key = api_key or os.environ.get("ANTHROPIC_API_KEY")
            default_model = "claude-opus-5-5"
        else:
            self.api_key = api_key or os.environ.get("NEBIUS_API_KEY")
            default_model = os.environ.get(
                "NEBIUS_MODEL", "nvidia/Llama-3_3-Nemotron-Super-49B-v1")
        self.model = model or os.environ.get("VERIFY_MODEL") or default_model
        self.name = f"llm:{self.model}"
        if self.api_key is None:
            raise RuntimeError(
                "No API key found. Set ANTHROPIC_API_KEY or NEBIUS_API_KEY. "
                "Every measured number in docs/verification-layer.md used "
                "RuleBackend and needs no key.")

    @staticmethod
    def _detect_provider():
        if os.environ.get("ANTHROPIC_API_KEY") is not None:
            return "anthropic"
        if os.environ.get("NEBIUS_API_KEY") is not None:
            return "nebius"
        return "anthropic"

    @staticmethod
    def _prompt(question, record):
        return ("You are checking the output of a scientific ranking "
                "pipeline.\n"
                f"{question}\n\n"
                f"Tool result as JSON:\n{json.dumps(record, indent=2)}\n\n"
                "Reply with one sentence a scientist would write, and say "
                "whether this is a problem and how serious it is.")

    def _body(self, prompt):
        if self.provider == "anthropic":
            return self.ANTHROPIC, {
                "content-type": "application/json",
                "x-api-key": self.api_key,
                "anthropic-version": "2023-06-01",
            }, {
                "model": self.model,
                "max_tokens": 2048,
                "output_config": {
                    "effort": "low",
                    "format": {"type": "json_schema", "schema": self.SCHEMA},
                },
                "messages": [{"role": "user", "content": prompt}],
            }
        return self.NEBIUS, {
            "content-type": "application/json",
            "authorization": f"Bearer {self.api_key}",
        }, {
            "model": self.model,
            "max_tokens": 2048,
            "response_format": {
                "type": "json_schema",
                "json_schema": {"name": "judgment", "strict": True,
                                "schema": self.SCHEMA},
            },
            "messages": [{"role": "user", "content": prompt}],
        }

    @staticmethod
    def _text(provider, payload):
        if provider == "anthropic":
            if payload.get("stop_reason") == "refusal":
                raise RuntimeError("the model declined this request")
            for block in payload.get("content", []):
                if block.get("type") == "text":
                    return block["text"]
            raise RuntimeError("no text block in the response")
        return payload["choices"][0]["message"]["content"]

    @staticmethod
    def _parse(text):
        text = text.strip()
        if text.startswith("```"):
            text = text.split("\n", 1)[1].rsplit("```", 1)[0]
        out = json.loads(text)
        return (bool(out["is_problem"]), out["severity"], out["sentence"])

    def _ask(self, question, record):
        url, headers, body = self._body(self._prompt(question, record))
        req = urllib.request.Request(
            url, data=json.dumps(body).encode(), headers=headers)
        with urllib.request.urlopen(req, timeout=120) as r:
            payload = json.loads(r.read())
        return self._parse(self._text(self.provider, payload))

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


def make_backend(name=None):
    """Pick a backend by name, or from VERIFY_BACKEND, defaulting to rules.

    This is the only place that chooses, so adding a key later is a
    flag and not an edit.
    """
    name = name or os.environ.get("VERIFY_BACKEND", "rules")
    if name == "rules":
        return RuleBackend()
    if name in ("llm", "anthropic", "nebius"):
        provider = None if name == "llm" else name
        return LLMBackend(provider=provider)
    raise ValueError(f"unknown backend {name!r}; use rules, llm, "
                     f"anthropic or nebius")
