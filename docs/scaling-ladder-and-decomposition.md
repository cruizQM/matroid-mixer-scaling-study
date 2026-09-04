# Does the cost-aware mixer actually scale, under assumptions closer to real feeders — and is any of it runnable today?

`docs/bounded-witness-mixer.md` established the bounded-witness mixer and
fixed its first, badly cost-blind default. This document picks up right
after that fix and follows the investigation to its actual conclusion —
which is not where it looked like it was heading partway through. In
order:

1. Headline result 1's own "cost decreases with network size" claim
   turns out to depend on an assumption (tie count held perfectly flat)
   that inverts under a mild, more realistic alternative.
2. An escalating ladder of increasingly realistic assumptions
   (short/long-range ties, log/linear tie-count growth) for the
   cost-aware mixer, checked against rough NISQ hardware feasibility.
3. A deliberate push to force the circuit cheap enough to be
   NISQ-relevant (50-500 CX) uncovers a real bug — a construction that
   looks safe because it has stopped doing anything.
4. An attempt to fix the search itself (a per-term adaptive cost
   coefficient) that looked like a clean win on an initial sweep, and
   was **not**, once properly re-validated after the bug fix above. Kept
   in the record as a real, instructive dead end.
5. The fix that actually works: combining the bounded-witness mixer with
   zone decomposition — validated cleanly, dominating the whole-graph
   construction everywhere it applies.

Every numeric claim below was re-measured after every code fix mentioned
in this document was already in place — an earlier draft of this
investigation compared numbers across a code change without re-running
both sides, and the "adaptive alpha wins" and "decomposition sometimes
loses" conclusions it produced were both wrong. See "A methodology
mistake, made and caught" near the end for what that looked like and how
it was caught, since it's as much a part of this repo's discipline as
any of the numbers.

## 1. Is headline result 1's flat-tie-count assumption load-bearing?

`run_scaling_study.py` holds `k_ties` fixed (`=3`) across its entire
sweep. `graphs.py`'s own docstring already notes this makes larger
networks progressively SPARSER, since the same small tie count gets
spread across more of the graph — which is exactly why headline result
1's circuit cost trends flat-to-*decreasing* with size. Whether tie count
actually stays flat as real feeders grow is a different question, and
one this repo had exactly one real anchor point for (5 ties at the
33-bus feeder) but had never checked against.

`run_scaling_study_log_ties.py` compares, at identical sizes/seeds, the
existing fixed `k_ties=3` against `k_ties_log(n) = round(1.43*ln(n))`
— calibrated so `k_ties_log(33) = 5`, matching that one real anchor
point exactly, not chosen to produce any particular result. Both use the
EXACT construction (`build_matroid_mixer`), unmodified — this section is
entirely about headline result 1's own methodology, independent of
anything in the bounded-witness/cost-aware work.

| | mean CX, n_nodes≤20 | mean CX, n_nodes≥30 | trend |
|---|---|---|---|
| fixed `k_ties=3` | 331 | 197 | **decreasing**, as headline result 1 reports |
| log-scaled `k_ties(n)` | 685 | 964 | **increasing** |

Averaged across the whole tested range (`n_nodes` 8-40), the log-scaled
condition costs **3.4x more** than the fixed condition; at `n_nodes=33`
— the exact size headline result 1's own sweep stops at, and the
calibration anchor — it's **5.8x more** (764 vs. 132 CX). The exact
construction stayed fully connected at every point in both conditions
(`k_ties` never exceeded 5, well below `tie_density_sweep.py`'s observed
failure threshold of ~6-8), so this is a genuine "more expensive, still
fully working" result, not confounded by connectivity loss.

**This doesn't make headline result 1 wrong** — it's a real, faithful
consequence of the assumption it states plainly. It does mean "cost
decreases with scale" is a property of that specific assumption, not of
the construction in general, and it hadn't been checked against even a
mild alternative before now.

## 2. An escalating ladder of realism for the cost-aware mixer

Given that headline result 1's own construction can look good or bad
depending on how tie count is assumed to grow, the natural next question
is whether the mixer actually meant for real-topology, real-scale use —
`truncated_mixer.py`'s bounded-witness construction — holds up better.

Four conditions, two independent axes, both calibrated to the same
single real anchor point (`k_ties(33)=5`):

- **tie placement**: short-range (`generate_feeder_graph`, headline
  result 1's family) vs. long-range (`generate_feeder_graph_long_range_ties`,
  the family `docs/circuit-validity.md`'s real-feeder section and
  the bounded-witness safety survey are built on).
- **tie-count growth**: log (`round(1.43*ln(n))`, mild) vs. linear
  (`round(0.152*n)`, aggressive — 23 ties at `n_nodes=150`).

At each size, a walked-exchange-graph sample
(`random_trees.random_walk_exchange_sample` — full enumeration is
intractable at these densities) feeds
`truncated_mixer.build_truncated_witness_mixer` at `cost_alpha=0.01`
(this repo's validated default — see section 4 for why NOT the adaptive
variant). Reproduce: `python scripts/run_cost_aware_scaling_ladder.py`.

| condition | n_nodes=10 | 30 | 60 | 100 | 150 |
|---|---|---|---|---|---|
| short_log | 264 | 1,303 | 3,191 | 2,592 | 2,745 |
| short_linear | 88 | 4,404 | 1,349 | 1,080 | 1,351 |
| long_log | 3,303 | 3,452 | 11,504 | 9,921 | 10,712 |
| long_linear | 4,182 | 3,452 | 11,980 | 10,431 | 13,273 |

(mean CX over 3 seeds; `results/cost_aware_scaling_ladder_results.csv`
has per-seed detail and depth.)

**The pattern**: cost roughly plateaus by `n_nodes=60` in every
condition rather than continuing to grow — placement (short vs. long
range) drives the cost LEVEL (short: hundreds to low thousands; long:
several thousand to low tens of thousands), growth rate (log vs. linear)
matters far less than placement does, and the mixer stays fully
connected at every single point tested. Variance is high at the smallest
size (`n_nodes=10`, coefficient of variation 50-65%+) and drops to a
much tighter band by `n_nodes≥60` — "small variance" is real, but only
past a floor size.

## 3. Is any of this executable on current hardware?

Rough estimate using published two-qubit gate fidelities (general
published ranges — best-case trapped-ion ~99.8-99.9%, i.e. ~0.1-0.2%
error per gate; typical superconducting ~99-99.5%, i.e. ~0.5-1% — not a
live lookup for a specific device), via the standard circuit-survival
approximation `fidelity ≈ (1-p)^N_CX`:

| CX count | best-case trapped-ion (p=0.001) | typical superconducting (p=0.005) |
|---|---|---|
| 100 | 90% | 61% |
| 500 | 61% | 8% |
| 1,000 | 37% | 0.7% |
| 5,000 | 0.7% | ~0 |
| 13,000 (this ladder's `long_linear` cases) | ~10⁻⁶ | ~0 |

None of section 2's conditions land anywhere near a NISQ-plausible
regime except the cheapest corners of the short-range conditions
(tens-to-low-hundreds of CX at the smallest sizes, before growth sets
in). Everything in the `long_*` conditions, and most of the
`short_*` conditions past `n_nodes=30`, needs fault tolerance, not
physical qubits, to run with any useful fidelity. This motivated
pushing cost down further, deliberately, even at the cost of safety —
section 4.

## 4. A deliberate push to 50-500 CX, a real bug, and a dead end

### The push and the collapse

The natural first lever is the search's own cost-pressure parameter.
Pushing `cost_alpha` (or, once it existed, the adaptive `gain_price`)
higher, combined with `exact_search_max_size=0` (routing every candidate
through the cost-aware path instead of letting brute-force exact search
handle some of them for free), DID reach 50-500 CX on the hardest tested
instance (`long_linear`, `n_nodes=150`) — but not for a good reason.
Direct inspection: **146 of 151 terms had collapsed to "never fire"** —
`unsafe_rate=0.0`, `mean_feasible_mass=1.0`, and `cx=20`, all because the
circuit had stopped doing almost anything, not because it preserves
feasibility well. `fully_connected=True` reported success on this
instance anyway — that flag is computed from which candidate pairs the
exchange-graph SAMPLE says *could* connect components, not from whether
the compiled circuit's terms actually fire, so it cannot see this
failure mode at all.

**Why "never fire" is an unbeatable shortcut without a fix**: a
majority-vote rule at witness size 0 is forced to pick "always fire" or
"never fire" (nothing in between, at size 0). Since a term is only ever
selected because `build_matroid_mixer`'s own union-find confirmed it
`connects_something`, most candidates' trigger sample is NOT majority-valid
— so "never fire" often has genuinely low LEAK (misclassifying only
the minority of states where it should have fired) at literally ZERO
circuit cost. Under any `leak + cost` scoring, that combination is
unbeatable once cost pressure is high enough, regardless of how the
pressure is applied.

### The fix

`_search_truncated_witness`/`_search_truncated_witness_adaptive` in
`truncated_mixer.py` now track, separately from the plain best-leak
option, the best leak found where `valid_patterns` is NON-empty
("active"), and restrict the final witness choice to active options
whenever at least one size offers one — an inert witness only wins if
truly NOTHING active was found at any size (a rare, genuinely
all-invalid case). Re-checking the same instance after the fix: 69 of
151 terms active, cost back up to ~9,500-10,800 CX, and **completely
unresponsive to further increases in cost pressure** — every setting
from `gain_price=0.001` to `3.0` gives the identical result, because
every term is already at its cheapest genuinely active option.

### The honest conclusion: there is no real 50-500 CX regime here

Checked across multiple conditions after the fix, every floor sits well
above 500 CX:

| condition (size) | floor CX, forced active | active/total terms |
|---|---|---|
| short_log, n=60 | ~2,300-2,900 | 19/19 (all) |
| short_linear, n=100 | ~4,600 | 26/36 |
| long_log, n=60 | ~27,700 | 152/206 |
| long_linear, n=150 | ~9,500-10,800 | 69/151 |

Every term that exists does so because it's genuinely needed to connect
otherwise-unreachable parts of the feasible space — a needed term can't
be made free, only made to fire on a real fraction of states, which
costs some minimum number of gates. The apparent 50-500 CX regime found
before the fix was "make most of the needed terms silently stop being
needed," not a real cost-safety tradeoff — the fixed search correctly
refuses to make that trade invisibly. Reaching a genuinely NISQ-relevant
gate count at this scale needs a structurally different fix, not further
tuning of this search's cost pressure — which motivated section 5.

Reproduce: `python scripts/run_cost_aware_scaling_ladder_aggressive.py`
(the collapse, `exact_search_max_size=0, gain_price=0.05`, tracks
`active_terms` explicitly per the fix above so the collapse is visible
in the CSV, not just in this narrative).

## 5. A dead end that looked like a win: adaptive per-term cost pressure

### The idea

A single global `cost_alpha` doesn't fit terms of different difficulty
— pushing it hard enough to help an easy term (little achievable leak
improvement, so width isn't worth paying for) also over-penalizes a hard
term (real improvement available, worth the width). The natural fix:
make the coefficient per-term, derived from how much THAT term actually
has to gain: `alpha_term = gain_price / total_achievable_gain`, where
`total_achievable_gain` is the term's own leak(0) minus its best
achievable leak at any size up to the cap. A term with a lot to gain
gets a small effective alpha (width stays affordable); a term with
little to gain gets a large one (settle for cheap).

An initial sweep, run BEFORE the never-fire fix in section 4 existed,
found this dominating a fixed `cost_alpha=0.01` — cheaper AND
safer-or-equal on every condition tested, no per-condition tuning
required. That result does not survive re-measurement.

### Why it looked like a win, and wasn't

The per-term formula turned out to be MORE susceptible to the
never-fire shortcut than a single shared coefficient, not less: `alpha_term`
inflates whenever a term's own achievable gain looks small — and the
never-fire option itself is exactly what makes a term's achievable gain
look artificially small (if "never fire" already has low leak, there's
"little room left to improve", which the formula reads as "this term
doesn't need width" rather than "this term found a free shortcut"). The
initial sweep's apparent win was measuring this artifact, at a milder
degree than section 4's explicit 146/151 collapse, but the same
mechanism.

Re-measured with the never-fire fix in place on both sides —
`python scripts/run_fixed_alpha_ladder.py` vs.
`python scripts/run_cost_aware_scaling_ladder.py` — the adaptive variant
is not a clean win. It wins narrowly on `short_log`, and loses,
sometimes badly, everywhere else:

| condition | n_nodes | fixed `cost_alpha=0.01` | adaptive `gain_price=0.01` | who wins |
|---|---|---|---|---|
| short_log | 30 | 4,404 | 1,303 | adaptive, 3.4x |
| short_linear | 100 | 1,080 | 7,889 | **fixed, 7.3x** |
| long_log | 30 | 3,452 | 17,225 | **fixed, 5.0x** |
| long_linear | 60 | 11,980 | 23,442 | **fixed, 2.0x** |

Safety is a wash-to-slightly-worse for adaptive throughout, not better,
in this corrected comparison. `build_truncated_witness_mixer`'s default
was reverted to `adaptive=False, cost_alpha=0.01` — the simpler option,
now that the more complex one isn't shown to earn its complexity.
`adaptive=True` remains available in the code (and
`_search_truncated_witness_adaptive`'s docstring documents this history
directly) for further investigation, not because it's recommended.

**Why this is in the document at all, not just reverted quietly**: this
repo's own discipline (`docs/circuit-validity.md`, `docs/bounded-witness-mixer.md`)
is to report what was tried and found wrong, not just what worked. A
plausible-sounding fix that turned out to be riding the same bug it was
supposed to be independent of is exactly the kind of thing that
discipline exists to catch and record.

## 6. The fix that actually works: decomposition

### The idea

Section 4's floors (2,300-27,700 CX) are driven by TERM COUNT (19-206
active terms), not individual gate width — most active terms already
have small witnesses. Zone decomposition (`zone_decomposition.py`,
`docs/circuit-validity.md`'s headline result 2 fix, originally built for
the EXACT construction) directly attacks term count: most of a
whole-graph construction's terms only exist to connect distant parts of
the network, something a small local zone doesn't need nearly as many
terms to do internally.

### Method

At each ladder instance, partition into zones
(`zone_decomposition.partition_zones_by_size`, `target_zone_size=8`,
matching `run_decomposition_scaling_study.py`'s own convention), then
build a `cost_alpha=0.01` bounded-witness mixer on EACH zone subgraph
and the assembly graph independently — exact enumeration when its own
cost (`C(n_edges, n_nodes-1)`) stays under a cap, falling back to a
walked sample only when it wouldn't (the assembly graph, one supernode
per zone with every cross-zone tie collapsed into a boundary edge, can
end up genuinely dense even at small `n_nodes` — one real instance had
19 zone-supernodes and 31 boundary edges, and brute-force enumeration
there alone took 81 seconds; a walked sample used unconditionally
everywhere, tried first, gave a meaningfully WORSE result on the same
instance than exact enumeration where exact was actually affordable, so
exact-when-cheap is preferred, not sampling by default). Total CX/depth
is summed across every zone + assembly subproblem; safety is checked
per-subproblem (the same scope `run_decomposition_scaling_study.py`
itself used for the exact construction's leak-freeness check, not a full
joint end-to-end simulation).

Reproduce: `python scripts/run_decomposed_cost_aware_ladder.py`.

### Result: it dominates, everywhere it applies

Comparing decomposed (`cost_alpha=0.01` per subproblem) against the
whole-graph construction from section 2 (same `cost_alpha=0.01`,
re-measured — see the methodology note below):

| condition | n=30 | n=60 | n=100 | n=150 |
|---|---|---|---|---|
| short_log | 5.6x | 14.4x | 18.7x | **26.4x** |
| short_linear | 5.6x | 10.0x | 9.8x | 5.3x |
| long_log | 16.7x | 17.4x | 4.2x | 8.2x |
| long_linear | 16.7x | 13.1x | 2.8x | 2.0x |

(cost reduction from decomposing, mean CX over 3 seeds;
`results/decomposed_cost_aware_ladder_summary.csv` /
`results/cost_aware_scaling_ladder_summary.csv`.) **Decomposition won
every single seed at every size tested where it actually applies.** At
`n_nodes=10` both methods tie exactly — the graph is too small to split
into more than one zone (`target_zone_size=8` → `round(10/8)=1`), so
decomposition trivially reduces to the whole-graph construction there.
`results/best_of_both_ladder_summary.csv` confirms this directly: a
"build both, keep whichever is cheaper" selector chose decomposition on
3/3 seeds at every size ≥30, and never once chose whole-graph where
decomposition had the chance to lose.

**Practical consequence**: the "build both and compare" selector is more
machinery than this result needs. A much simpler rule — decompose
whenever the graph splits into more than one zone — captures essentially
all of the benefit, since whole-graph never won a single seed once given
the chance to lose. This is the actual recommended construction:
`truncated_mixer.build_truncated_witness_mixer` on each zone + assembly
subproblem from `zone_decomposition.py`, with `cost_alpha=0.01` (or
`0.0` — section 5's finding that adaptive doesn't help says nothing
about whether ANY cost pressure at all helps within a zone; that
narrower question wasn't re-tested here), not the whole-graph
construction directly, for any network too large to be a single zone.

## A methodology mistake, made and caught

Worth stating plainly, since catching it is part of what makes the
final numbers here trustworthy: an earlier pass through this
investigation compared the whole-graph and decomposed constructions, and
separately the adaptive and fixed search variants, using result files
generated at DIFFERENT points in the code's history — some from before
the never-fire fix (section 4), some from after. Both comparisons looked
clean and were wrong: decomposition appeared to sometimes LOSE to
whole-graph on easy instances (it doesn't — that whole-graph baseline
was quietly benefiting from undetected term collapse), and adaptive
search appeared to dominate fixed `cost_alpha` (it doesn't — same root
cause). Both were caught by noticing a suspiciously large jump in a
recomputed number (`short_linear n=60`'s whole-graph cost going from 187
to 3,164-5,448 CX on the identical instance) rather than by any
systematic check — a large the effect that happened to be large enough
to be conspicuous, not something the process caught automatically. Every
comparison in this document was re-run with the never-fire fix already
in place on every side being compared, specifically because of this.

## Honest scope of this document

- Section 1's log/linear tie-count growth models are both calibrated to
  ONE real data point (the 33-bus feeder's 5 ties) — neither is a
  validated growth law for how real feeder tie counts actually scale
  with network size, just two differently-shaped curves through the one
  anchor available.
- Section 6's safety check is per-subproblem, not a full joint
  simulation of the assembled whole-network mixer — the same scope
  limitation `docs/circuit-validity.md`'s own decomposition work has for
  the exact construction.
- Section 6's decomposition result used `cost_alpha=0.01` specifically
  (not adaptive, and not re-swept across other fixed values) on each
  subproblem — whether a different fixed coefficient, or genuinely no
  cost pressure at all, does even better WITHIN a zone (zones are much
  smaller, so the never-fire dynamics found in section 4 may or may not
  reappear at that scale) is not tested here.
- The NISQ feasibility numbers (section 3) are order-of-magnitude
  illustrations from published two-qubit gate fidelity ranges, not a
  claim about any specific current device's real specifications.
