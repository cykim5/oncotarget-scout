# Checking the pipeline's own output, and whether splitting the job up is worth it

## Why this exists

A real defect in this project survived three phases of my own review.
The score renormalized the layer weights over whichever layers a gene
had, so a gene with two layers could outrank a gene with five. Every
one of the top 20 genes was missing the dependency layer. The decision
shortlist was unaffected, because it already required dependency
support, but the headline ranking was wrong for weeks and I read it
many times without seeing it.

That is the honest motivation. Not that automated verification is
fashionable, but that a human reading his own output missed something a
cheap invariant catches in five milliseconds.

Both rankings are committed, the broken one and the fixed one, so the
same checks can be run against each and the difference measured instead
of asserted.

## The shape

Three specialists, each with a narrow remit, its own tools, and no
knowledge of what the others are doing.

| Agent | What it asks | Tool |
| --- | --- | --- |
| `layer-balance` | Is the top of the list systematically missing a layer the rest of the list has? | `layer_depletion` |
| `controls` | Do the checks fixed before the ranking ran still pass? | `answer_key_position`, `shortlist_dependency_support` |
| `card-fidelity` | Do the written reason cards match the table they came from? | `card_numbers_match` |

Then a fourth, `cross-check`, reads what the first three produced and
looks for a single cause behind separate failures. It is the one part
of the design that cannot exist without decomposition, because a single
pass producing one flat list has nothing independent to compare.

The tools return numbers and never a judgment. The judgment is a
backend: `RuleBackend` applies thresholds, `LLMBackend` sends the same
tool result to a model and asks for the same decision. Separating them
is what makes the two comparable.

### The invariant that matters

`layer_depletion` is deliberately general. It does not know what any
layer means. For each layer it compares how often that layer is present
in the top k against how often it is present across the whole ranked
list, and gives the one sided binomial p value for the top being
depleted. If the top of a ranking is missing evidence the population
has, the score is rewarding absence.

Writing a check that looked for non-coding RNAs, or for the specific
bug, would have proved nothing. This one finds it without being told
what it is looking for.

## Result 1: detection

Run against the ranking as it stood before the fix:

| Agent | Finding |
| --- | --- |
| `layer-balance` | 0 of the top 20 carry dependency, against 59.0 percent of the ranked population, p = 1.8e-08 |
| `layer-balance` | 1 of the top 20 carry malignant share, against 46.3 percent, p = 7.3e-05 |
| `controls` | KRAS at the 85.5th percentile, EGFR at the 81.7th, CD274 at the 65.8th |
| `cross-check` | two independent checks fail together, which points at one cause in how the score is built rather than at the data |

Six critical findings. Against the fixed ranking, one.

**That one is a false positive, and it stays in.** The `controls` agent
flags CD274 at the 72nd percentile, because the rule inherited this
project's original wording, which expected every approved target near
the top. The project itself later established that CD274 ranking low is
the correct result, since it predicts who responds to immunotherapy
rather than naming a gene this subgroup depends on. Silencing it by
naming CD274 in the rule would be tuning the checker to the answer.

It is also the most interesting thing the experiment produced. The
checker reproduced a disagreement the project had already resolved in
prose but never encoded. Writing down what should be true turns out to
be harder than running the check.

## Result 2: what decomposition bought

Exactly one finding. The cross-check synthesis, and nothing else. The
single pass and the decomposed run found the same five underlying
problems, because they were given the same tools and the same backend.

That is a smaller return than the architecture suggests, and it is the
honest number.

## Result 3: what is actually running in parallel, and what it buys

Worth being precise, because "parallel execution" can mean almost
anything. Here it is one `ThreadPoolExecutor`. Nothing else in the
project is concurrent.

Two places it can go, and they do not perform the same.

**Agent level.** One thread per specialist, three threads. Each
specialist still judges its own records one after another, so the
clock is set by the busiest specialist.

**Judgment level.** `collect()` is cheap and local, so it stays per
specialist. Every record those specialists produced then goes into one
pool. The imbalance between specialists stops mattering.

| Agent | Judgments |
| --- | --- |
| `layer-balance` | 5 |
| `controls` | 4 |
| `card-fidelity` | 0 |

Nine judgments in sequence. The busiest specialist has five, so agent
level parallelism cannot beat 9/5 = 1.80x however many threads it is
given. Judgment level is bounded by nine.

| Cost per check | Single pass | By agent | By judgment | Best speedup |
| --- | --- | --- | --- | --- |
| 0 ms, rule thresholds | 5.6 ms | 5.9 ms | 5.9 ms | **0.94x** |
| 1 ms | 17.4 ms | 7.1 ms | 7.2 ms | 2.42x |
| 10 ms | 121.5 ms | 65.7 ms | 22.5 ms | 5.40x |
| 50 ms | 506.6 ms | 270.3 ms | 71.4 ms | 7.10x |
| 100 ms | 961.1 ms | 520.7 ms | 126.9 ms | 7.57x |
| 500 ms | 4,553 ms | 2,524 ms | 527.9 ms | **8.62x** |

Three things fall out of that table.

**At rule speed, coordination loses.** 0.94x. Thread setup costs more
than the work. The simpler approach wins and the honest thing is to
say so rather than report the 8.62x and stop.

**Agent level decomposition converges on 1.80x, exactly its ceiling.**
Adding a fourth specialist would not move it. The limit is the shape of
the work, not the number of workers.

**Where you cut matters more than whether you cut.** Same checks, same
backend, same threads: 1.80x one way and 8.62x the other. Agent level
decomposition is the intuitive design and the one most agent
frameworks push you toward, and here it leaves a factor of five on the
floor.

### The caveat that matters

These are threads, not processes, so while Python holds the GIL only
one of them computes at a time. Thread parallelism only buys wall clock
when the work releases the GIL. A network call to a model does. A scipy
threshold does not.

So the 0 ms row is real measured CPU work and the rows below it are a
simulation: the delay is `time.sleep`, which releases the GIL exactly
as a network call would. That makes it a fair model of a model backed
check and it is still a model. The speedups above are what the
architecture should deliver against an API, not what it has delivered
against one.

## What this does and does not show

It shows that a general invariant, costing five milliseconds, catches a
defect that three phases of human review missed, and stays quiet on the
corrected ranking apart from one false positive whose cause is
understood.

It does not show that coordination beats a simpler approach in general.
At this cost profile it does not: the parallelism is overhead and the
decomposition bought one extra finding. The conditions under which it
would pay are now measured rather than assumed, which is a more useful
thing to know than a speedup number taken on its own.

Three limits, stated plainly:

- **The measured runs used `RuleBackend`.** `LLMBackend` is implemented
  against the Messages API and has not been run, so nothing here rests
  on it. A model backed check costs hundreds of milliseconds at least,
  which is the bottom of the table above, where judgment level
  parallelism wins clearly on time. Whether a model finds *more* than a
  threshold does is the open question, and it needs a key to answer.
- **The speedups below the first row are simulated**, with `time.sleep`
  standing in for a network call. The GIL note above says why that is a
  fair model and still a model.
- **One defect is one defect.** The experiment has a sample size of one
  historical bug. It tests whether these checks catch that, not whether
  they catch the next one.
- **The checks encode what I already believe should be true.** That is
  the real limit, and the CD274 false positive is the proof of it.

## Running it

```bash
python3 scripts/run_verification_experiment.py   # the measured results above
python3 -m pytest tests/test_verify.py -q        # including the regression test
```

Nothing above needs an API key. **There is no model anywhere in the
measured results.** The specialists are Python functions and the
judgments are thresholds written in `scripts/verify/backends.py`. The
architecture is agent shaped, meaning the work is decomposed, each part
has its own tools, the parts run in parallel and one part reads the
others, but nothing in it reasons.

### Adding a model later

One flag, no edit:

```bash
export ANTHROPIC_API_KEY=...        # or NEBIUS_API_KEY
python3 scripts/run_verification_experiment.py --backend llm
```

It picks the provider from whichever key is set, writes to
`results/verification-experiment-llm.json` so the committed numbers are
not overwritten, and skips the latency sweep, which simulates cost with
`sleep` and is meaningless once every call costs real time.

`VERIFY_MODEL` picks the model. The default is Claude Opus 5.5.
`claude-haiku-5-5` costs about a cent for a whole run and is the more
interesting test: if the cheapest model matches thresholds written by
hand, that is the result worth reporting, and if it is less consistent
than a threshold on the borderline CD274 case, that is worth reporting
too.
