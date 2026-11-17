"""A backend wrapper that costs time, so the architecture can be priced.

With rule thresholds every check finishes in microseconds, and running
three of them in parallel costs more in thread overhead than it saves.
A model backed check costs seconds. This wrapper adds a fixed delay to
every judgment so the crossover point can be measured instead of
asserted.
"""
import time


class SlowBackend:
    def __init__(self, inner, seconds):
        self.inner = inner
        self.seconds = seconds
        self.name = f"{inner.name}+{seconds * 1000:.0f}ms"

    def _slow(self, method, record):
        time.sleep(self.seconds)
        return method(record)

    def judge_depletion(self, record):
        return self._slow(self.inner.judge_depletion, record)

    def judge_answer_key(self, record):
        return self._slow(self.inner.judge_answer_key, record)

    def judge_card(self, record):
        return self._slow(self.inner.judge_card, record)

    def judge_shortlist(self, record):
        return self._slow(self.inner.judge_shortlist, record)


class CountingBackend:
    """Counts how many judgments each run makes, to explain the speedup."""

    def __init__(self, inner):
        self.inner = inner
        self.name = inner.name
        self.calls = 0

    def _count(self, method, record):
        self.calls = self.calls + 1
        return method(record)

    def judge_depletion(self, record):
        return self._count(self.inner.judge_depletion, record)

    def judge_answer_key(self, record):
        return self._count(self.inner.judge_answer_key, record)

    def judge_card(self, record):
        return self._count(self.inner.judge_card, record)

    def judge_shortlist(self, record):
        return self._count(self.inner.judge_shortlist, record)
