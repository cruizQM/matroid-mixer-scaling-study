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
   worth the trade. This document measures both.

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
qubit count).

| n_nodes | n_qubits | approx terms (mean, of total) | unsafe rate (mean, min-max over 3 seeds) | worst feasible mass observed | mean feasible mass (mean over 3 seeds) |
|---|---|---|---|---|---|
| 10 | 14 | 28.7 / 28.7 (100%) | 0.433 (0.343–0.513) | 0.640 | 0.977 |
| 30 | 34 | 92.7 / 96.7 (96%) | 0.241 (0.183–0.327) | 0.375 | 0.968 |
| 60 | 64 | 80.3 / 149.7 (54%) | 0.224 (0.143–0.320) | 0.633 | 0.980 |
| 100 | 104 | 43.0 / 133.0 (32%) | 0.331 (0.293–0.367) | 0.283 | 0.958 |
| 150 | 154 | 19.0 / 123.7 (15%) | 0.157 (0.063–0.240) | 0.756 | 0.991 |

(`generate_feeder_graph_long_range_ties`, `K_TIES=5`, 3 seeds/size, 300
independent Wilson's-algorithm starting trees traced per seed —
`results/bounded_witness_safety_survey.csv` (per seed) /
`_summary.csv` (per size). "unsafe" = any measurable feasible-mass loss
(`< 1 - 1e-6`) by the end of the term sequence; "worst feasible mass" is
the single lowest value observed across all 900 traced trajectories at
that size.)

**Two real, and slightly in tension, patterns.** `unsafe_rate` trends
down with scale overall (43% → 16% from `n_nodes=10` to `150`) but *not*
monotonically — `n_nodes=100` is a real outlier (33%, all three seeds
agree, so this isn't one unlucky seed), tracking a plausible mechanism:
`unsafe_rate` follows how much of the construction ended up approximate
at all (`n_approximate_terms`/`n_terms`, shown above), which itself
depends on how many candidates happened to have small exact witnesses at
each random instance — not a quantity this survey controls directly.
`n_nodes=150`'s construction happened to need few approximate terms (15%
of the total); `n_nodes=100`'s needed more (32%) — that's the more
likely driver of the outlier than network size itself. Reported as
measured, not smoothed into a monotonic story it doesn't quite tell.

**What stays consistent, and is arguably the more load-bearing number**:
`mean_feasible_mass` — the average, not worst-case, final feasible mass
across all 900 trajectories at a given size — stays in a tight 0.96-0.99
band at every size tested, regardless of the `unsafe_rate` swings above.
Most "unsafe" trajectories lose only a little mass, not most of it;
the worst-case column shows real tail risk exists (down to 0.28 at one
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
