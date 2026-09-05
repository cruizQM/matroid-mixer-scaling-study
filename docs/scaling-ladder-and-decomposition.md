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
   construction almost everywhere it applies.
6. A second instance of the SAME bug shape found inside the decomposition
   script itself, fixed, and every dependent result re-run.
7. Hierarchical (density-aware) decomposition, developed to push the
   still-expensive hard conditions further toward NISQ feasibility --
   three iterations, each empirically validated against the last.
8. A real-topology check: does actual (or real-benchmark) feeder data
   support the log-growth and long-range-tie modeling choices used
   throughout? (Yes to log over linear, clearly; yes to long-range
   ties.) `CONDITIONS` is updated accordingly -- linear tie-count growth
   is kept as an explicit stress test, not deleted, but is no longer the
   default.
9. Direct validation of the actual construction (decomposition +
   cost-aware bounded-witness mixer) on two real networks, not just
   synthetic proxies calibrated to real anchor points.
10. A cost-CAPPED decomposition that measures each subproblem's actual
    cost and recursively re-splits anything too expensive, instead of
    picking a zone size up front and hoping — guarantees every
    synthetic-ladder seed stays under 500 CX. A first pass on real data
    got CIGRE MV cleanly but missed IEEE33 by 1.6%; checking whether that
    was actually the best achievable (not just "does it meet the
    threshold") found two real improvements -- searching more
    `cost_alpha` values, and more zone-size granularities -- that closed
    the gap on IEEE33 and improved CIGRE MV further, with zero regressions
    across the full synthetic ladder.

Every numeric claim below was re-measured after every code fix mentioned
in this document was already in place — two SEPARATE earlier drafts of
this investigation each compared numbers across a code change without
re-running both sides, and each produced a wrong conclusion as a result
("adaptive alpha wins" once, "decomposition sometimes loses" the other
time). See "Methodology mistakes, made and caught" near the end for what
both looked like and how they were caught, since that's as much a part
of this repo's discipline as any of the numbers.

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

### Result: it dominates almost everywhere it applies

Comparing decomposed (`cost_alpha=0.01` per subproblem) against the
whole-graph construction from section 2, both measured under IDENTICAL,
current code (see section 7 below for why that qualifier matters —
an earlier version of this table used a decomposed baseline that was
silently running the rejected adaptive search inside each zone):

| condition | n=10 | n=30 | n=60 | n=100 | n=150 |
|---|---|---|---|---|---|
| short_log | ties | 5.6x | 14.4x | 18.7x | **26.4x** |
| long_log | 0.44x (worse) | 35.2x | 8.3x | 46.7x | 4.9x |

(cost reduction from decomposing, mean CX over 3 seeds;
`results/decomposed_cost_aware_ladder_summary.csv` /
`results/cost_aware_scaling_ladder_summary.csv`. `short_linear`/
`long_linear` are no longer the default conditions — see section 9 —
but the same comparison on them, kept as a stress-test data point, shows
the identical pattern: decomposition wins every seed at `n_nodes>=30`.)

`long_log`'s ratio swings widely by size (4.9x to 46.7x) — real
seed-to-seed variance on only 3 seeds per point, not a sign the
technique is unreliable: **per-seed data
(`results/best_of_both_ladder_results.csv`) shows decomposition winning
EVERY SINGLE SEED at every size >= 30, for both conditions, without
exception** — the swings are in how much it wins by, not whether it
wins. At `n_nodes=10`, `long_log` shows decomposition looking WORSE
(0.44x) — this is a real, but different, effect, not decomposition
failing: with only one zone at this size (`target_zone_size=8` →
`round(10/8)=1`), decomposition should reduce to the identical
construction, but the two measurement paths use different tree inputs
(the whole-graph ladder uses a walked-exchange-graph SAMPLE; the
decomposed path uses EXACT enumeration whenever affordable, which it is
at this size) — comparing them at `n_nodes=10` is actually comparing two
different-but-valid inputs to the same search, not measuring
decomposition's own effect. At `n_nodes=30` and up, both zones and the
assembly graph are genuinely non-trivial, and the comparison is a fair
one.

**Practical consequence**: decompose whenever the graph splits into more
than one zone — that's the recommended construction:
`truncated_mixer.build_truncated_witness_mixer` on each zone + assembly
subproblem from `zone_decomposition.py`, with `cost_alpha=0.01`, not the
whole-graph construction directly, for any network too large to be a
single zone.

## 7. A second bug, same shape as the first, found the same way

`run_decomposed_cost_aware_ladder.py`'s `measure_subproblem` was written
between sections 5 and 6, chronologically -- before the adaptive-vs-fixed
investigation in section 5 concluded adaptive doesn't reliably beat fixed
`cost_alpha`. It still called `build_truncated_witness_mixer` with
`adaptive=True`, meaning every decomposition result up to this point
(including section 6's original table, and hierarchical decomposition
below) was quietly built on the REJECTED search mode, not the one this
document itself recommends.

Confirmed to matter, not just theoretically inconsistent: on one
instance (`long_linear`, `n_nodes=150`, `seed=0`), one specific zone (of
19) cost 3,550 CX under adaptive vs. 292 CX under fixed
`cost_alpha=0.01` -- a 12x gap, and the dominant contributor to that
instance's entire zone-level cost. Fixed to `adaptive=False,
cost_alpha=0.01`, and every dependent script re-run (not left stale):
`run_decomposed_cost_aware_ladder.py`,
`run_hierarchical_decomposed_ladder.py`, `run_best_of_both_ladder.py`.
Section 6's table above already reflects the corrected numbers. The net
effect on the headline comparison was mixed, not uniformly better or
worse -- `long_log` at `n_nodes=150` got WORSE (flat decomposition:
1,318 -> 2,226 CX) while `long_linear` got BETTER (8,767 -> 3,160 CX) --
matching the same condition-dependent pattern section 5 already found at
whole-graph scale.

## 8. Hierarchical (density-aware) decomposition

Flat (one-level) decomposition still leaves `long_log` and `long_linear`
at real scale (`n_nodes=150`) short of NISQ-plausible (2,226 and 3,160
CX respectively, post-fix). Direct probing found naively shrinking
`target_zone_size` does NOT help -- it just moves cost into an
increasingly dense ASSEMBLY graph (one supernode per zone, one edge per
cross-zone tie: shrinking zones means more edges cross zone boundaries,
so the assembly graph gets denser, not simpler). The fix has to be
recursive: if the assembly graph is itself still large, decompose IT the
same way, instead of measuring it directly.

Three iterations, each validated against the last before moving on:

1. **Recurse on the assembly graph, same `target_zone_size` at every
   level.** Clean win for `long_log` (2,226 -> 1,526 CX at n=150 with the
   corrected fixed-alpha baseline), but UNSTABLE for `long_linear`:
   8,798 +/- 5,943 CX across 3 seeds, one seed reaching 17,176 -- WORSE
   than not recursing at all. Traced directly: the assembly graph is
   proportionally DENSER than the graph it came from (a 19-node assembly
   graph can have 27-34 edges, far denser than the 150-node/172-edge
   original), so a fixed zone size gives too few, too-concentrated
   sub-zones there -- one seed's density landed almost entirely in a
   single one of only 2 resulting sub-zones.
2. **Halve `target_zone_size` at each recursion level.** Fixes the
   instability (`long_linear`: 4,308.7 +/- 1,768.7, no outliers) but is
   still an arbitrary schedule, not a response to what's actually being
   decomposed.
3. **Scale `target_zone_size` inversely to the CURRENT graph's measured
   density relative to the original** (`graph.n_edges / graph.n_nodes`,
   recomputed at every level) -- measurably better than blind halving on
   the same seeds (`long_linear`: 3,767.3 +/- 1,167.2, tighter variance
   too), because it responds to actual density instead of a fixed
   per-level schedule.

Final validated numbers at `n_nodes=150` (post section-7 fix, primary
conditions): `long_log` 1,526 +/- 236 CX (vs. 2,226 flat, a 1.46x
improvement; 21.7% fidelity at best-case trapped-ion error rates vs.
flat's 10.8%); `long_linear` (stress test) 1,986.7 +/- 1,252.8 CX (vs.
3,160 flat, a 1.59x improvement; 13.7% fidelity vs. 4.2%). Neither
crosses into a comfortably NISQ-plausible number, but both are real,
validated improvements over flat decomposition's own corrected baseline.

Reproduce: `python scripts/run_hierarchical_decomposed_ladder.py`.

## 9. Does real topology actually support these modeling choices?

Sections 1-8 calibrate log/linear tie-count growth and long/short-range
tie placement to a SINGLE real anchor point (5 ties at the 33-bus
feeder). This checks that choice against real data directly, not just
assumes it.

**Tie count**: three real/real-benchmark networks (via `pandapower.networks`
-- CIGRE MV benchmark, `case33bw` already used throughout this repo, and
`mv_oberrhein`, a real German MV distribution network), spanning a 12x
size range:

| network | n_bus | n_ties | ties/n_bus |
|---|---|---|---|
| CIGRE MV | 15 | 3 | 0.200 |
| case33bw | 33 | 5 | 0.152 |
| mv_oberrhein | 179 | 6 | 0.034 |

The ratio drops 6x from smallest to largest -- flatly inconsistent with
LINEAR growth (which would keep it roughly constant) and reasonably
consistent with LOG growth (`ties/log(n_bus)` comes out to 1.11, 1.43,
1.16 -- much tighter). **This is why `CONDITIONS` (section 2, and
`run_cost_aware_scaling_ladder.py` generally) now defaults to the two
log-growth conditions only** -- the linear-growth ones move to
`STRESS_TEST_CONDITIONS`, kept and fully reproducible, but no longer
presented as equally realistic.

**Tie range**: measuring actual topological span (shortest-path distance
between a tie's endpoints in the rest of the network) and normalizing by
network size:

| | n=15 | n=33 | n=179 |
|---|---|---|---|
| real data | ~0.33-0.45 | ~0.36-0.61 | ~0.35-0.44 (within each sub-component) |
| this repo's long-range generator | 0.60 | 0.55 | 0.38 |
| this repo's short-range generator | 0.14 | 0.065 | 0.011 (shrinking with scale) |

Real ties consistently span 33-45% of the network -- genuinely
long-range. The long-range generator (`generate_feeder_graph_long_range_ties`,
already this repo's choice for the real-feeder failure mode since
`docs/circuit-validity.md`) matches reasonably well, converging closely
by `n=179`. The short-range generator -- what headline result 1's own
"cost decreases with scale" story is built on -- is nowhere close, and
gets WORSE with scale, not better.

**Honest caveats on this check**: only 3 data points spanning one order
of magnitude -- enough to rule out linear and lean toward log, not
enough to fit a real growth law. None is a clean single radial feeder
either: `mv_oberrhein` is actually 2 substations' worth of sub-feeders
combined (removing its 6 ties splits the 177-edge remainder into 69
components, not 1 -- several ties are jointly load-bearing for
connectivity there, not simple independent redundant loops the way
case33bw's are), and CIGRE MV similarly isn't a single connected tree
without also counting its 2 transformers as edges. Real published data
tends to model service areas, not single feeders.

## 10. Direct validation on real networks, not just synthetic proxies

Everything above uses synthetic graphs calibrated to real anchor points.
This runs the actual construction -- decomposition + the cost-aware
bounded-witness mixer -- directly on two real networks: `case33bw`
(already this repo's headline result 2 validation instance) and a newly
added loader, `real_feeders.load_cigre_mv` (same discipline as
`load_ieee33`: ties identified from open switches / `in_service=False`,
checked to actually form a spanning tree before use, not assumed).
`mv_oberrhein` is not used here for the same reason section 9 flags it
as structurally awkward -- running the actual circuit construction on it
would need real data-cleaning work (e.g. splitting it into its two
actual service areas) beyond this check's scope, the same call this
repo already made for the CATS transmission dataset
(`real_feeders.py`'s own docstring).

| network | construction | n_terms | CX | depth | result |
|---|---|---|---|---|---|
| CIGRE MV (15 buses, 3 ties) | exact whole-graph | 16 | 12,220 | 24,245 | connected, exactly leak-free, 37 candidates dropped |
| CIGRE MV | decomposed (2 zones) | 10 | **3,668** | 7,101 | 3.3x cheaper; unsafe_rate 3.0%, mean feasible mass 0.9995 |
| IEEE33 (33 buses, 5 ties) | exact whole-graph *(already measured, `results/real_feeder_results.csv`)* | 24 | 96 | 78 | **573 of 597 candidates dropped, disconnected** |
| IEEE33 | decomposed (4 zones) | 6 | **132** | 272 | fully functional; unsafe_rate 0%, mean feasible mass 1.0 |

On IEEE33, decomposition isn't just cheaper than the exact whole-graph
construction -- the exact construction doesn't work at all there
(already established by headline result 2; not re-measured here, cited
instead). This reproduces that same distinction with this branch's
cost-aware construction rather than the exact one, on genuine published
topology, not a synthetic stand-in.

Reproduce: `python scripts/run_real_networks_hierarchical.py`.

## 11. Cost-capped decomposition: guaranteeing the 500 CX target, not hoping for it

Section 8's hierarchical decomposition, and section 6's flat version
before it, both pick a `target_zone_size` up front (fixed, or
density-scaled) and hope the result is cheap. It usually is, but not
reliably: measuring the actual variance in section 6/8's own results
found decomposed cost's coefficient of variation reaching 99% at some
sizes (`long_log`, `n_nodes=100`) -- far WORSE than the whole-graph
construction's own variance, which shrinks smoothly with size (55% down
to 5%). The reason: decomposition trades one big averaging problem for
`n_zones` small independent ones, and a single "unlucky" zone (one
concrete case: 3,550 of a 3,566-CX subtotal from ONE zone of 19) can
dominate a seed's total when there aren't yet enough zones for that to
average out.

**The fix: measure, don't guess.** Build each subproblem, transpile it,
and check its ACTUAL cost against a threshold (`CX_THRESHOLD=500`). If
it's already under threshold, it's a leaf, done. If not, partition IT
further (same `zone_decomposition` machinery, at a smaller
`target_zone_size`) and recurse on each of its own zones + assembly.
Falls back to escalating `cost_alpha` (trading safety for cost, bounded
-- never pushed to the point of inertness section 4 found) only for
subgraphs too small to partition further (`n_nodes <= 4`) and still over
threshold.

**Result on the synthetic ladder: EVERY seed, both primary conditions,
all 5 sizes, meets the threshold** — `results/cost_capped_decomposition_summary.csv`
shows `all_seeds_met_threshold=True` on every single row, with mean
feasible mass exactly 1.0 everywhere (no fallback cost-pressure
escalation was needed on any instance tested). This also fixes the
variance problem as a side effect: capping every leaf means no single
unlucky zone can dominate the total the way it did before.

**Tested directly on the two real networks (not just synthetic proxies) —
first pass**: CIGRE MV met the target cleanly at 274 CX (didn't even need
to split, mean feasible mass 0.985). IEEE33 missed by a small margin: 508
CX, 8 over (1.6%). Traced precisely, not left as a mystery: decomposing
one level further bottomed out at a 4-node, 8-edge sub-assembly (a
near-complete multigraph core) sitting exactly at the `n_nodes <= 4`
irreducibility threshold, where `cost_alpha` genuinely couldn't help
(508 CX at every tested value from 0.01 to 5.0 -- those candidates were
already at the witness width where leak is exactly zero). That diagnosis
was correct as far as it went, but it was answering the wrong question --
see below.

### Two refinements, found by asking "is this actually the best we can
do?" of the real-network results, not just "does it meet the threshold?"

**1. Don't stop at the first `cost_alpha` that passes.** The algorithm
above tries only `cost_alpha=0.01` and accepts it the instant it meets
threshold. Checking whether a higher value could do even better on
CIGRE MV (never tried, because 0.01 already "worked"): `cost_alpha=0.2`
gives **64 CX at a PERFECT 1.0 mean feasible mass** -- cheaper AND safer
than the reported 274 CX / 0.985, left on the table purely because the
search stopped at the first passing attempt. `_best_under_threshold` now
sweeps all of `FALLBACK_COST_ALPHAS` and keeps the cheapest result that
still meets threshold, at every leaf -- not just the irreducible-fallback
ones that already used this list for a different purpose.

**2. Don't jump straight to a halved granularity -- but don't
assume the un-halved one is always better either.** The 508-CX diagnosis
above was real, but incomplete: it explained why *that specific*
recursive path couldn't do better, without asking whether a *different*
path would have. Direct check: partitioning the whole IEEE33 graph at
`target_zone_size=8` directly (this repo's own flat-decomposition
convention, `run_decomposed_cost_aware_ladder.py`'s own choice) --
instead of the halved `target_zone_size=4` the algorithm jumped to the
instant the whole graph failed to meet threshold -- gives a 4-zone split
whose assembly graph is ALSO 4 nodes and 8 edges, same size as the
508-CX core, but a **different, cheaper graph**: 132 CX. Same size,
4x cheaper, because which specific zones get contracted together differs
depending on the granularity chosen, and the halving schedule never
tried the granularity that turns out to matter here.

The tempting fix -- "always try the current, un-halved size before
halving" -- was checked against the FULL synthetic ladder before being
kept, not just the one real instance that motivated it, and that check
caught a real problem: it made 12 of the 30 synthetic ladder seeds
**worse**, by up to 15x (`short_log, n_nodes=30, seed=0`: 28 -> 416 CX).
Shrinking zones sometimes helps and sometimes hurts depending on where a
graph's real density sits -- there is no universally-better fixed rule,
exactly as section 8's own hierarchical-decomposition findings already
established. The actual fix: `_try_granularity` now tries BOTH the
current `target_zone_size` and its half, fully recurses each candidate,
and `decompose_to_threshold` keeps whichever total is cheaper -- the same
measure-don't-guess discipline this section is built on, applied to the
granularity choice itself rather than assumed in either direction.

**Result of both refinements together, re-verified on the full ladder
(not just the two real networks that motivated them)**: every one of the
30 synthetic-ladder seeds is now equal to or cheaper than the
pre-refinement result -- none regressed -- with several dramatic
improvements (`long_log, n_nodes=30, seed=0`: 312 -> 24 CX;
`long_log, n_nodes=150, seed=0`: 500 -> 368 CX, previously sitting
exactly AT the threshold with no margin). On the two real networks:

| network | total CX (before -> after) | safety |
|---|---|---|
| CIGRE MV | 274 -> **64** | 0.985 -> **1.0 (perfect)** |
| IEEE33 | 508 -> **132** | 1.0 (perfect), unchanged |

IEEE33's cost-capped refinement now exactly matches plain zone
decomposition (both 132 CX) rather than costing 3.85x more -- it searches
multiple granularities and the one flat decomposition already uses turns
out to be the best one available, so there's nothing left to improve on
beyond matching it. CIGRE MV's refinement is now strictly better than
before on both cost and safety. Both real networks now show the exact
same monotonic pattern the synthetic ladder always did: exact ->
decomposed -> cost-capped never gets more expensive, on either network,
with no exceptions -- the "less clean on real data" finding this section
originally reported was an artifact of an under-searched algorithm, not
a real property of real topology.

Reproduce: `python scripts/run_cost_capped_decomposition.py` (synthetic
ladder), `python scripts/run_real_networks_hierarchical.py` (the two
real networks, now with `cost_capped` as a third method alongside
`exact_whole_graph` and `decomposed`).

## Methodology mistakes, made and caught

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
systematic check — a large effect that happened to be large enough to be
conspicuous, not something the process caught automatically. Every
comparison in this document was re-run with the never-fire fix already
in place on every side being compared, specifically because of this.

**A second instance of the same pattern, section 7**: the decomposition
script itself was quietly running the rejected adaptive search inside
every zone for several commits before this was noticed — again not
caught by any systematic check, but by directly diagnosing why one
`long_linear` seed's hierarchical-decomposition result looked
anomalously bad, which led back to a single zone costing 12x more than
it should have. Section 6's table, and every decomposition-dependent
result in sections 7-8, were re-run after this fix; the
`best_of_both_ladder` results were also stale until refreshed for the
same reason. Two independent instances of "a result built on code that
changed after the result was generated" in one investigation is worth
taking as a standing risk in this kind of iterative work, not a one-off
-- re-running dependents after any change to shared code
(`truncated_mixer.py`, `zone_decomposition.py`) is now treated as
mandatory, not optional, for exactly this reason.

## Honest scope of this document

- Section 1's log/linear tie-count growth models were originally
  calibrated to ONE real data point (the 33-bus feeder's 5 ties);
  section 9 checks this against 3 real/real-benchmark networks and finds
  log growth clearly better supported than linear, but 3 points spanning
  one order of magnitude is enough to rule out linear, not enough to fit
  an actual growth law.
- Section 6/8's safety checks are per-subproblem, not a full joint
  simulation of the assembled whole-network mixer — the same scope
  limitation `docs/circuit-validity.md`'s own decomposition work has for
  the exact construction. Section 10's real-network results carry the
  same limitation.
- Section 6's decomposition result used `cost_alpha=0.01` specifically
  (not re-swept across other fixed values) on each subproblem — whether
  a different fixed coefficient, or genuinely no cost pressure at all,
  does even better WITHIN a zone (zones are much smaller, so the
  never-fire dynamics found in section 4 may or may not reappear at that
  scale) is not tested here.
- The NISQ feasibility numbers (section 3) are order-of-magnitude
  illustrations from published two-qubit gate fidelity ranges, not a
  claim about any specific current device's real specifications.
- Section 8's density-aware hierarchical decomposition was only tested
  on `long_log`/`long_linear` (the conditions it was built to fix) --
  whether it also helps (or is even needed) for `short_log`/`short_linear`
  was not checked, since flat decomposition already gets those to a much
  cheaper regime on its own.
- Section 9's real-topology tie-count/range check and section 10's
  direct real-network validation both rely on `pandapower`'s bundled
  network data and this document's own switch/`in_service` parsing logic
  for identifying which lines are ties -- checked directly against this
  repo's own already-validated `case33bw` extraction (`real_feeders.load_ieee33`)
  as a sanity check, but not independently re-verified against the
  original published papers for CIGRE MV or `mv_oberrhein`.
- Section 11's refinements close IEEE33's original 508-CX gap by trying
  more zone SIZES (two candidates per split, keep the cheaper total) and
  more `cost_alpha` values, not by changing which nodes get grouped
  together. A genuinely smarter zone-CHOICE strategy (partitioning that
  reasons about which specific nodes end up together, not just how many
  per zone) remains untried and could still help further, particularly
  on graphs where neither tried granularity happens to land well.
- Section 11's `CX_THRESHOLD=500` was chosen to match this document's
  own earlier NISQ-feasibility discussion (section 3), not re-derived
  from a specific target device's current published error rates -- it's
  a round, illustrative number in the same spirit as section 3's
  estimates, not a claim that exactly 500 is the true cutoff for any
  particular piece of hardware.
