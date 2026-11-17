"""Three ways to run the same checks, so they can be compared.

run_single_pass     one worker, every check in sequence. The simpler
                    approach the others have to beat.
run_decomposed      one thread per specialist. Each specialist still
                    judges its own records one after another, so the
                    speedup is capped by the busiest specialist.
run_flat_parallel   every judgment from every specialist goes into one
                    pool. The imbalance between specialists stops
                    mattering, and the cap becomes the pool size.

All three get identical tools and an identical backend, so a difference
between them is the shape and nothing else.

A note on what parallel means here. These are threads, not processes,
so while Python holds the GIL only one of them computes at a time.
Thread parallelism only buys wall clock when the work releases the GIL,
which a network call to a model does and a scipy threshold does not.
That is why the measured speedup at rule speed is below 1.0 and why the
latency sweep in docs/verification-layer.md is the honest way to read
these numbers.
"""
import time
from concurrent.futures import ThreadPoolExecutor

from .agents import CrossCheckAgent, specialists


def _finish(mode, backend, findings, start, agents, workers):
    return {"mode": mode, "backend": backend.name, "findings": findings,
            "seconds": time.perf_counter() - start, "agents": agents,
            "workers": workers}


def run_single_pass(ctx, backend, k=20):
    agents = specialists(backend, k)
    start = time.perf_counter()
    findings = []
    for agent in agents:
        findings = findings + agent.run(ctx)
    return _finish("single-pass", backend, findings, start,
                   ["single-worker"], 1)


def run_decomposed(ctx, backend, k=20, max_workers=4):
    agents = specialists(backend, k)
    start = time.perf_counter()
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        results = list(pool.map(lambda a: a.run(ctx), agents))
    findings = [f for part in results for f in part]
    findings = findings + CrossCheckAgent(backend).run(ctx, findings)
    return _finish("decomposed", backend, findings, start,
                   [a.name for a in agents] + ["cross-check"], max_workers)


def run_flat_parallel(ctx, backend, k=20, max_workers=12):
    """Parallelize at the judgment level instead of the agent level.

    collect() is cheap and local, so it runs per specialist. Every
    record those specialists produced is then judged in one pool, which
    is where the cost is when the backend is a model.
    """
    agents = specialists(backend, k)
    start = time.perf_counter()
    work = []
    for agent in agents:
        for kind, record in agent.collect(ctx):
            work.append((agent, kind, record))
    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        judged = list(pool.map(lambda w: w[0].judge(ctx, w[1], w[2]), work))
    findings = [f for f in judged if f is not None]
    findings = findings + CrossCheckAgent(backend).run(ctx, findings)
    return _finish("flat-parallel", backend, findings, start,
                   [a.name for a in agents] + ["cross-check"], max_workers)


def summarize(result):
    counts = {}
    for f in result["findings"]:
        counts[f.severity] = counts.get(f.severity, 0) + 1
    return {"mode": result["mode"], "backend": result["backend"],
            "n_findings": len(result["findings"]), "by_severity": counts,
            "seconds": round(result["seconds"], 3)}
