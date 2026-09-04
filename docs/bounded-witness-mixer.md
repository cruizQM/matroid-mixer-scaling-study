# Bounded-witness mixer: a second failure axis, and an alternative to decomposition

`docs/circuit-validity.md` already establishes that the exact whole-graph
construction's witness-search bound (`size ≤ 2`, found on this repo's own
synthetic sweep) "does NOT hold in general," tracing the failure on real
data to tie **range** — real tie switches connect distant parts of a
feeder, giving long fundamental cycles — and fixing it with zone
decomposition. This document covers two more pieces, developed in a
related internal research effort and ported back here on this branch:

1. **A second, independent failure axis**: tie **density** alone, holding
   range and network size fixed, breaks the same witness-search bound —
   on this repo's own nearest-tie graph family, not just the long-range
   one.
2. **An alternative to decomposition**: for cases where a bounded witness
   is preferred over splitting the problem into zones, `truncated_mixer.py`
   accepts a fixed-size witness cap and the resulting leakage, rather than
   dropping the candidate exchange outright — provided the leakage stays
   small enough, and rare enough across realistic starting points, to be
   worth the trade.
3. **A circuit-cost investigation and fix**: the first version of (2) was
   measured (not assumed) to produce circuits costing up to 38x more than
   necessary, for a construction whose whole point was practicality — a
   real gap between "this preserves feasibility reasonably well" and
   "this is cheap enough to be worth building," found only once a cost
   measurement was actually pointed at it. A fix (`cost_alpha`, a
   cost-aware search objective) closes most of that gap and is now this
   construction's default. This document covers all three.

## Finding 1: tie density breaks the exact bound on its own

`mixer.py`'s own module docstring already states the size-≤-2 bound is
"plausibly because a sparse graph whose extra edges are geometrically
short-range has few alternative ways to reconnect around any given edge" —
a claim about the nearest-tie generator's fixed, *small* `k_ties`, not
about sparsity in general. `tie_density_sweep.py` (this branch) tests that
directly: fix `n_nodes=8` and the nearest-tie selection rule (unchanged
from `generate_feeder_graph`), and grow `k_ties` from 2 to 10.

| k_ties | n_qubits | \|B\| (spanning trees, mean) | max witness (mean/max) | dropped candidates (mean/max, of *n_candidates*) | any seed disconnected? |
|---|---|---|---|---|---|
| 2 | 9 | 9 | 0.75 / 2 | 0 / 0 (of 36) | no |
| 4 | 11 | 58 | 3.0 / 4 | 0 / 0 (of 55) | no |
| 6 | 13 | 287 | 3.25 / 4 | 39.5 / 71 (of 78) | **yes** |
| 8 | 15 | 1173 | 2.0 / 4 | 73.5 / 101 (of 105) | **yes** |
| 10 | 17 | 3348 | 0.75 / 1 | 129.75 / 136 (of 136) | **yes** |

(`n_nodes=8` fixed throughout, 4 seeds per `k_ties`, brute-force witness
search only — `prefer_structural=False`, `MAX_WITNESS_SEARCH_SIZE=4`, this
repo's own defaults, unchanged, for a fair comparison against the
short-cycle bound `mixer.py`'s docstring already claims for this
generator. Full per-seed data: `results/tie_density_sweep_results.csv`;
per-`k_ties` summary: `results/tie_density_sweep_summary.csv`.)

**The failure doesn't show up primarily as growing witness size — it
shows up as growing incompleteness.** Only `k_ties=4` and `k_ties=6` push
`max_witness_size` up against the search cap (3-4); by `k_ties=10`, mean
max witness is back down to 0.75. That's not the bound holding — it's the
search *giving up before finding a large witness at all*: `dropped_candidates`
climbs from 0% of candidates (`k_ties≤4`) to 51% (`k_ties=6`) to 70%
(`k_ties=8`) to 95% (`k_ties=10`), and by `k_ties=6` at least one of the 4
seeds already leaves the mixer **disconnected** (`fully_connected=False`)
— every seed is disconnected by `k_ties=8`. The true minimal witnesses at
this density plausibly exceed the size-4 cap entirely (the same situation
`docs/circuit-validity.md`'s real-feeder section found on long-range ties,
just reached here via density instead of range) — brute force at a larger
cap could in principle find them, but at combinatorially growing cost
(`elapsed_s` already reaches 12-15s per single `k_ties=10` instance at
only `n_qubits=17`, using the search algorithm's own vectorized
implementation).

**A parallel honest note, matching this repo's own pattern in
`docs/circuit-validity.md`'s real-feeder section**: don't read the
`max_witness_size` column alone as "the bound holds better at high
density" — a shrinking observed-witness figure next to a collapsing
`dropped_candidates`/`fully_connected` figure is the search failing
faster and more completely, not the construction getting cheaper.

This is a **different** failure from the one `docs/circuit-validity.md`
documents: nothing about tie *range* changed here (still nearest-neighbor
by construction), and the graph is still sparse by the repo's own
definition (`|E| = n_nodes - 1 + k_ties`, `k_ties` held fixed within each
column). What changed is only how many ties there are. Zone decomposition,
as built (`zone_decomposition.py`), targets the range failure specifically
— it shrinks each zone so any tie *inside* it looks short-range again, but
says nothing about a zone that is itself tie-dense. Whether decomposition
also happens to help against density (by shrinking each zone's tie count
along with everything else) is a fair question this document doesn't
settle; what it does establish is that density is a distinct axis worth
checking on its own, not automatically covered by "already fixed the
range problem."

Reproduce: `python scripts/tie_density_sweep.py` →
`results/tie_density_sweep_results.csv` / `_summary.csv` (plus
`_plot.png` if `matplotlib` is installed — the numbers above don't depend
on it).

## Finding 2: a bounded-witness alternative to decomposition

`truncated_mixer.py` (`build_truncated_witness_mixer`) is built directly
on top of `mixer.py`'s own machinery — same candidate-pair selection loop,
same `ExchangeTerm`/`MixerConstruction` types, same
`_minimize_patterns`/`exchange_term_circuit` circuit synthesis — with one
change: when `find_witness_set`'s brute-force search finds no exact
witness within a small cap, instead of dropping the candidate (what
`build_matroid_mixer` does), it falls back to a random-restart search over
the graph-derived cycle-union candidate pool for a witness of a fixed,
capped size, using a **majority-vote** validity rule instead of an
exactly-verified one. The resulting term still compiles through the exact
same circuit synthesis path — `verify_term_no_leakage` still checks
whether the compiled circuit implements the *declared* `valid_patterns`
exactly (it does; nothing about circuit synthesis changed). What's new is
that the declared rule itself is no longer guaranteed semantically
correct against the true feasible set — by design, in exchange for never
dropping a connecting candidate.

**A real bug this cross-check caught.** Porting this construction onto
`mixer.py`'s shared `exchange_term_circuit` surfaced a genuine correctness
bug in code this repo already had: a `witness_qubits == ()` shortcut that
fired the swap block unconditionally, which is only correct if a
zero-qubit witness always means "always valid" — true for the exact
construction (where `find_witness_set` only returns an empty witness when
every trigger state shares the same, already-known-nonempty, validity),
but not for `truncated_mixer.py`'s majority-vote search, which can
legitimately land on "zero witness qubits, majority says never fire."
`verify_leakage_trace.py` caught this by cross-checking the actual
compiled circuit's behavior against the declared `valid_patterns`, not
assuming the shortcut was fine because it looked fine on the
exact-construction case it was written for — fixed in
`mixer.exchange_term_circuit` (previous commit on this branch), with
`verify_correctness.py` re-run to confirm no regression on the exact
construction.

### Does the declared leakage rate actually cost feasible mass?

A term's majority-vote leakage rate is a property of the *abstract*
validity function, measured uniformly over the sample it was built from.
Whether it costs any probability mass on a *specific* trajectory (one
starting tree, evolving through the mixer's full term sequence) is a
different, generally smaller, question — leakage from one term can land
on a state a later term's rule happens to classify correctly, or never
get revisited at all. `leakage_trace.py` traces this exactly (not
sampled): `trace_danger_mass` (full statevector, small instances only)
and `trace_danger_mass_sparse` (dict-of-active-basis-states, validated
exactly equivalent to the dense method in `verify_leakage_trace.py`, and
the only one of the two that's tractable at real-feeder qubit counts,
since the reachable subspace from one starting tree stays small even
where `2**n_qubits` does not).

### A circuit-cost investigation: staying connected isn't the same as being cheap

Everything above is about *leakage* — whether the mixer stays on feasible
states. Staying CONNECTED where the exact construction fails (verified
above and in the safety survey below) is a different claim from "the
resulting circuit is cheap," and the two were never checked together
until this investigation. They should have been: this repo's entire
methodology (headline results 1 and 2) is built around transpiled CX
count and depth, not just correctness or connectivity — a construction
this document was about to call "an alternative to decomposition"
without ever measuring the one thing that phrase implies.

**The first measurement was bad.** `build_truncated_witness_mixer`'s
original default (`max_witness_size=6`, no cost awareness) produces
circuits costing **17,386–32,520 CX gates at just 15-17 qubits**
(`n_nodes=8`, the density-failure family) — for comparison, headline
result 1's worst case anywhere in this repo was 1,556 CX at 10 qubits,
and the real 33-bus decomposed construction topped out around 7,904 CX
on a 37-qubit graph. A 17-qubit circuit costing 2-20x more than the
worst case ever measured on a *37-qubit real feeder* is not a
scaling win.

**Root cause, verified directly (not assumed)**: `_search_truncated_witness`
selects witnesses purely by leakage. Adding a witness qubit can only
weakly *reduce* leakage, never increase it (it refines the majority-vote
partition) — so with no cost penalty, the search has every reason to walk
to `max_witness_size` every time, and no reason not to. Checking one
dense instance directly: all 16 terms used witness size 5-6, each with
1-9 minimized cubes (`_minimize_patterns` — the exact construction's own
Quine-McCluskey machinery, reused here unchanged — IS working, pattern
*count* is small) but most cubes still command 5-6 controls. A generic
k-controlled 2-qubit gate's synthesis cost grows steeply in `k`, and
minimization shrinks cube *count*, not cube *width* — it only collapses
two patterns into a don't-care when they differ in one bit and agree
everywhere else, far more likely for the exact construction's
matroid-structured validity functions than for a majority-vote rule over
a fairly arbitrary qubit subset. So minimization, which fixed the exact
construction's own cost problem (`docs/circuit-validity.md`), does not
fix this one — same underlying principle (witness size alone doesn't
predict circuit cost), a construction it hadn't been checked against
before.

**Does lowering the cap fix it? Partially, and unevenly.** A direct sweep
of `max_witness_size` (2-6) on the density-failure family shows cost
scales steeply and cleanly with the cap while safety barely moves:

| k_ties | cap=2 CX | cap=6 CX | cost ratio | cap=2 mean feasible mass | cap=6 mean feasible mass |
|---|---|---|---|---|---|
| 6 | 759 | 13,276 | 17.5x | 0.969 | 0.997 |
| 8 | 543 | 8,914 | 16.4x | 0.984 | 0.995 |
| 10 | 1,084 | 25,020 | 23.1x | 0.953 | 0.983 |

(4 seeds/point, `n_nodes=8`; `results/truncated_witness_cap_sweep_results.csv`.)
Checked against the actual family this document's real-scale claims are
built on (`generate_feeder_graph_long_range_ties`) and across the full
size range (`n_nodes` 10-150, `results/truncated_witness_cap_sweep_longrange_results.csv`),
the same qualitative pattern holds, but the SIZE of the win shrinks
sharply with scale — because at large `n_nodes` most terms already find a
cheap exact witness in the build sample, so only a shrinking minority
ever pay the wide-witness cost at all:

| n_nodes | cap=2 CX | cap=6 CX | cost ratio |
|---|---|---|---|
| 10 | 2,321 | 88,937 | 38.3x |
| 60 | 10,404 | 53,205 | 5.1x |
| 150 | 9,461 | 17,722 | 1.9x |

Pushing the cap lower still (0, 1) does NOT continue this trend cleanly —
`results/truncated_mixer_cap_0_1_results.csv` shows the low end is
unstable and sometimes non-monotonic: at one density instance, cap=0 gave
*perfect* safety (unsafe_rate=0, mean mass=1.0 on every seed) while cap=1
was *worse* (unsafe_rate=0.51) at higher cost; at another, cap=0's mean
unsafe_rate (0.64) beat cap=1's (0.80). Below the cap 2-3 range, the
search's leakage estimate (computed on the small `trigger` set for one
candidate pair) stops reliably predicting generalization to fresh
starting trees. A single global cap, in other words, is a blunt
instrument: pushed low it saves real money, pushed too low it gets
unreliable, and there's no single number that's right for every term.

**Making the search cost-aware, instead of picking one global cap.**
`_search_truncated_witness` now accepts `cost_alpha`: score each
candidate witness as `leak + cost_alpha * 2**size` instead of `leak`
alone (`2**size` — a cheap proxy for the controlled-gate cost width
actually drives, not a transpiled measurement, which would make the
search itself expensive). `cost_alpha=0.0` reproduces the original,
cost-blind behavior exactly. Sweeping `cost_alpha` at a fixed
`max_witness_size=6` (letting the search choose its own size per term,
rather than forcing every term to the same cap) on the same density
instances finds points that *dominate* the fixed-cap curve above —
better or comparable safety at a fraction of the cost, not just a
different point on the same tradeoff:

| k_ties | fixed cap=6 (baseline) | cost-aware, `cost_alpha=0.01` |
|---|---|---|
| 6 | cx=13,276, mass=0.995, unsafe=0.16 | cx=1,065 (12x cheaper), mass=0.988, unsafe=0.26 |
| 8 | cx=8,914, mass=0.995, unsafe=0.17 | cx=651 (14x cheaper), mass=**0.996**, unsafe=**0.11** |
| 10 | cx=25,020, mass=0.983, unsafe=0.54 | cx=1,665 (15x cheaper), mass=0.978, unsafe=**0.30** |

(`results/truncated_mixer_cost_aware_results.csv`.) At `k_ties=8`,
cost-aware search beats the most expensive fixed-cap setting on cost,
safety, AND unsafe-rate simultaneously — because it lets terms that don't
need width skip paying for it, instead of forcing every term through the
same global ceiling. Push `cost_alpha` too high (0.05 in the same sweep)
and it does occasionally misfire on harder instances (one `k_ties=10`
seed dropped to mean mass 0.90) — there is a real sweet spot, not "more
alpha is free." `cost_alpha=0.01` is now `build_truncated_witness_mixer`'s
default (`max_witness_size` can stay generous — the objective itself now
discourages using it unnecessarily, rather than a hard cap forcing every
term to the same width regardless of whether it's needed).

Reproduce: `python scripts/measure_truncated_mixer.py`,
`python scripts/truncated_witness_cap_sweep.py`,
`python scripts/truncated_witness_cap_sweep_longrange.py`,
`python scripts/truncated_mixer_search_refinement.py`.

**This is not the end of the story.** A follow-up investigation
(`docs/scaling-ladder-and-decomposition.md`) checked this fixed
`cost_alpha=0.01` default against an escalating ladder of more realistic
assumptions, tried replacing it with a per-term ADAPTIVE coefficient
that looked like a clean win on an initial sweep and was not once
properly re-validated (a real, instructive dead end, kept in the record
rather than erased), found and fixed a genuine bug along the way (a
search objective that can make most of a mixer's terms silently stop
firing at all), and landed on the fix that actually works: combining
this construction with zone decomposition, which dominates the
whole-graph version everywhere tested. Read that document for what
happened after this one.

### Real-scale safety survey

`run_bounded_witness_safety_survey.py` builds a truncated mixer on
`generate_feeder_graph_long_range_ties` — the family `docs/circuit-validity.md`
already established breaks the *exact* whole-graph construction — at
increasing `n_nodes`, using a random-walk-on-the-exchange-graph sample
(`random_trees.random_walk_exchange_sample`) to build the term set (full
enumeration is exactly what's intractable at this scale). It then draws
`N_STARTING_TREES` independent, exactly-uniform starting trees via
Wilson's algorithm (`random_trees.random_spanning_tree` — independent of
the sample used to build the mixer, so this measures generalization to
fresh starting points, not just trees the construction was built from),
and traces each one's exact final feasible mass via
`leakage_trace.final_feasible_mass` (`sparse=True`, required at this
qubit count). Also transpiles the resulting circuit (`TRANSPILE_BASIS`,
`optimization_level=1`, matching this repo's own convention) at each
point, so cost and safety are reported together as one picture of the
construction actually being recommended — this script originally
measured safety alone and left cost as a separate, easy-to-miss finding;
that gap is what the circuit-cost investigation above was about.

**These numbers use `cost_alpha=0.01`, the construction's current
default** (not the original cost-blind behavior) — this is the
recommended construction's real-scale profile, not a "before" snapshot.

| n_nodes | n_qubits | approx terms (mean, of total) | mean CX count | unsafe rate (mean, min-max over 3 seeds) | worst feasible mass observed | mean feasible mass (mean over 3 seeds) |
|---|---|---|---|---|---|---|
| 10 | 14 | 28.7 / 28.7 (100%) | 2,169 | 0.232 (0.180–0.317) | 0.786 | 0.990 |
| 30 | 34 | 92.7 / 96.7 (96%) | 3,452 | 0.048 (0.000–0.137) | 0.868 | 0.997 |
| 60 | 64 | 80.3 / 149.7 (54%) | 11,106 | 0.143 (0.073–0.243) | 0.634 | 0.986 |
| 100 | 104 | 43.0 / 133.0 (32%) | 10,231 | 0.260 (0.227–0.300) | 0.326 | 0.962 |
| 150 | 154 | 19.0 / 123.7 (15%) | 9,755 | 0.120 (0.063–0.193) | 0.756 | 0.992 |

(`generate_feeder_graph_long_range_ties`, `K_TIES=5`, 3 seeds/size, 300
independent Wilson's-algorithm starting trees traced per seed —
`results/bounded_witness_safety_survey.csv` (per seed) /
`_summary.csv` (per size). "unsafe" = any measurable feasible-mass loss
(`< 1 - 1e-6`) by the end of the term sequence; "worst feasible mass" is
the single lowest value observed across all 900 traced trajectories at
that size.)

**Compared to the pre-fix numbers (`cost_alpha=0.0`, the original
default)**: mean CX dropped 1.8x (`n_nodes=150`) to 41x (`n_nodes=10`) —
the same shrinking-with-scale pattern the circuit-cost investigation
found on the density family, for the same reason (more terms find a
cheap exact witness at large scale, so cost-awareness has less to fix).
Safety did not get worse anywhere, and mostly got measurably *better* —
`unsafe_rate` at `n_nodes=30` fell from 24.1% to 4.8%, and
`mean_feasible_mass` improved or held at every single size. There was no
tradeoff to make here; the cost-aware default was strictly better on
every metric this survey tracks, on the actual family and size range
this document's real-scale claims are built on.

**Two real, and slightly in tension, patterns remain in the safety
numbers themselves.** `unsafe_rate` trends down with scale overall but
*not* monotonically — `n_nodes=100` is still a real outlier (26%, all
three seeds agree, so this isn't one unlucky seed), tracking a plausible
mechanism: `unsafe_rate` follows how much of the construction ended up
approximate at all (`n_approximate_terms`/`n_terms`, shown above), which
itself depends on how many candidates happened to have small exact
witnesses at each random instance — not a quantity this survey controls
directly. `n_nodes=150`'s construction happened to need few approximate
terms (15% of the total); `n_nodes=100`'s needed more (32%) — that's the
more likely driver of the outlier than network size itself. Reported as
measured, not smoothed into a monotonic story it doesn't quite tell.

**What stays consistent, and is arguably the more load-bearing number**:
`mean_feasible_mass` — the average, not worst-case, final feasible mass
across all 900 trajectories at a given size — stays in a tight 0.96-1.0
band at every size tested, regardless of the `unsafe_rate` swings above.
Most "unsafe" trajectories lose only a little mass, not most of it; the
worst-case column shows real tail risk exists (down to 0.33 at one
`n_nodes=100` instance) but is the exception, not the typical outcome.

Reproduce: `python scripts/run_bounded_witness_safety_survey.py` →
`results/bounded_witness_safety_survey.csv` (per seed) /
`_summary.csv` (per size).

## Honest scope of this document

- The density sweep (Finding 1) uses `n_nodes=8` throughout, kept small
  enough for `enumerate_spanning_trees` to stay tractable across the
  whole `k_ties` range — it establishes that density is *a* failure axis,
  not how it scales with network size on top of density. Extending it to
  larger `n_nodes` at fixed density (mirroring
  `run_decomposition_scaling_study.py`'s own scaling sweep) is future
  work, not done here.
- The safety survey (Finding 2) builds each mixer from a *sample* of
  trees (the exchange-graph walk), not the full feasible set — the same
  statistical-vs-exhaustive caveat `docs/circuit-validity.md`'s
  decomposition work doesn't have to make, since `enumerate_spanning_trees`
  stays exact there. What stays exact regardless: `verify_term_no_leakage`
  confirms the compiled circuit implements the *declared* rule exactly;
  the safety survey is a separate, additional check of whether the
  declared rule itself holds up, not a substitute for that per-term
  circuit-correctness guarantee.
- Neither finding claims the bounded-witness mixer beats zone
  decomposition, or vice versa — they solve different problems (bounded
  leakage vs. exact partitioning) and this document doesn't run a QAOA
  loop to compare their downstream feasible-sampling rates end to end.
  That comparison exists in a related internal project, not included
  here (against a penalty-QAOA baseline, not against decomposition), but
  hasn't been re-run against *this* repo's decomposition construction
  specifically.
  What the circuit-cost investigation DOES now establish, that wasn't
  known before: at low `n_nodes` on the density-failure family, the
  cost-aware truncated mixer's circuit cost (543-1,065 CX at 15-17
  qubits, `cost_alpha=0.01`) is competitive with or cheaper than the
  exact construction's cost on the FEW instances where the exact
  construction still worked at all (e.g. one connected `k_ties=6` seed
  cost 5,258 CX) — so "alternative to decomposition" is now a claim with
  a cost number behind it, at least at this scale, not just a
  connectivity claim.
- `cost_alpha=0.01` was chosen from a sweep of 6 values
  ({0, 0.0005, 0.001, 0.005, 0.01, 0.05}) on 3 seeds each of 3 density
  instances (`n_nodes=8`) — it consistently sat at or near the point that
  dominated the fixed-cap curve there, and the real-scale safety survey
  (long-range family, `n_nodes` up to 150) confirms it doesn't hurt
  safety at that default when applied at scale. It was NOT independently
  re-tuned on the long-range family, or per-instance — a per-family or
  per-instance tuned `cost_alpha` might do better still; `0.01` is a
  reasonable single default, not a claimed optimum.
- The circuit-cost investigation's generalization checks (density family
  at `n_nodes=8`, long-range family at `n_nodes` 10-150) both used
  `optimization_level=1` transpilation, matching this repo's own
  convention throughout — headline result 2's own investigation found
  `optimization_level=3` only closes 1-3% of a similar gap, so this
  isn't expected to change the conclusion, but it wasn't re-checked here
  specifically for the truncated construction.
