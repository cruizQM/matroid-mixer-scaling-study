# Circuit validity: what was checked, what broke, and how it was fixed

This document exists because the first version of this mixer construction
was wrong, and the fix materially shapes the results in `/results`. It's
included so an evaluator can see the reasoning, not just the final code.

## The initial (incorrect) assumption

The matroid basis-exchange axiom guarantees that swapping one tree edge `e`
for one non-tree edge `f` (when valid) always keeps the state a spanning
tree, and that connectivity of the whole exchange graph follows from
single-edge swaps alone. By analogy with simpler constraints, where a
standard XY-hopping gate `(X_eX_f+Y_eY_f)/2` is *unconditionally* safe (it
annihilates `|00>` and `|11>`, so it does nothing when neither or both of
`e,f` are active), the first version of this code applied that same
unconditional gate to every selected exchange pair.

## Why that's wrong here

Whether exchanging `e` for `f` actually produces a valid spanning tree
depends on whether `f` crosses the cut created by removing `e` from the
*current* tree -- a property of the tree's global structure, not just the
values of qubits `e` and `f`. An unconditional swap gate fires on *every*
state with exactly one of `e,f` active, including trees where the exchange
is invalid.

This was caught by directly checking the built circuit's unitary against
every spanning tree, not assumed to be fine because the reasoning sounded
plausible. Concretely, on a 6-node/7-edge test instance, exchange
`(edge 2, edge 3)` (graph edges `(2,3)` and `(2,4)`) had 5 trigger states
(states with exactly one of the two edges present); 1 of them mapped
outside the feasible set under the unconditional gate.

## Characterizing the fix: is the needed correction local or global?

The natural worry is that fixing this needs a full graph-connectivity
oracle (whether `f` reconnects the two components of `T - e`), which could
require circuitry whose cost grows with graph size and would undermine the
scaling claim this repo exists to support.

Direct search says no, **on the synthetic graphs this study's scaling
sweep uses**: for every selected exchange on every synthetic instance
checked, there is a *minimal witness set* -- a small set of other qubits
whose values alone determine validity -- of size 0, 1, or 2. Never larger,
and not growing with instance size in anything tested there.
`scripts/mixer.py`'s `find_witness_set` finds this witness set by bounded
brute-force search (capped at size 4; if none is found within the cap,
that candidate exchange is dropped rather than included with an incorrect
circuit -- on every *synthetic* instance tested, 0 candidates were ever
dropped for this reason).

**This bound does NOT hold on the real IEEE 33-bus topology** -- see
"Real-feeder validation: this bound fails, and why" below, added after
running the same construction against real data. Read that section before
citing the "size ≤ 2" claim above as general; it is a property of this
study's synthetic graph family, not of graphic-matroid mixers in general.

**Plausible reason this stays small on the synthetic family**: those
graphs are sparse by construction (a spanning-tree backbone plus a small,
fixed number `k_ties` of extra edges -- see `methodology.md`), and —
critically — `k_ties` extra edges are chosen as the *geometrically
nearest* non-tree pairs, which keeps each one's fundamental cycle (the
tree-path between its endpoints) short: 2-3 tree edges, empirically,
regardless of instance size (see the real-feeder section below for the
measurement). A short fundamental cycle means few alternative ways to
reconnect around any given edge, so "does this exchange cross the right
cut" reduces to checking a small, fixed number of alternative edges. This
was correctly flagged as "a plausible mechanism, not a proven theorem" in
earlier versions of this document — and it turned out not to generalize,
which is exactly why that hedge mattered.

## Real-feeder validation: this bound fails, and why

Running the identical construction (`scripts/run_real_feeder_validation.py`)
against the real, published IEEE 33-bus topology (Baran & Wu 1989, via
`pandapower.networks.case33bw`; 37 qubits, 50,751 spanning trees) produced
a materially different result from every synthetic instance tested:
**573 of 630 real candidate exchanges (91%) were dropped** — no witness
of size ≤ 4 found — and the resulting mixer is **not fully connected**
(`fully_connected=False` in `results/real_feeder_results.csv`). The 24
terms that were selected are still exactly correct (verified leakage-free
via `verify_all_terms_no_leakage`, all with witness_size=0), but the
basis-exchange graph they form does not reach every spanning tree from
every other one — a real completeness gap for QAOA, not just a
performance number.

**Root cause, measured directly, not assumed**: `scripts/
investigate_fundamental_cycles.py` computes, for a spanning tree, the
tree-path length between each tie edge's endpoints (its fundamental cycle
length — exactly the set of edges whose configuration can affect an
exchange's validity). Result:

| | IEEE 33-bus (real) | Synthetic (nearest-tie model) |
|---|---|---|
| max fundamental cycle length | **20** tree edges | 3 tree edges |
| mean | **11.8** | 2.2 |
| min | 6 | 2 |

Real published tie switches exist specifically to link **distant** parts
of a feeder for reconfiguration redundancy — the opposite bias from
"nearest non-tree pair." A witness set generally needs to encode which of
a fundamental cycle's member edges are actually present in the current
tree, so validity for an exchange tied to a 20-edge cycle can depend on up
to ~20 other qubits, far beyond the search cap of 4 (and `C(35,20)` is
astronomically infeasible to brute-force search regardless of cap).

**This means brute-force minimal-witness search is the wrong algorithm
for real feeder topology, not merely under-tuned.** No feasible increase
to `MAX_WITNESS_SEARCH_SIZE` fixes this.

![A tie edge's fundamental cycle, short vs long](../results/illustration_fundamental_cycle.png)

The mechanism, illustrated directly (`scripts/plot_illustrations.py`):
same backbone, same node positions, one nearest-neighbor tie edge versus
one long-range tie edge. A candidate exchange's validity can depend on
every tree edge on its tie edge's fundamental cycle — 2 edges on the
left, 8 on the right in this small example; 6-20 on the real graph, as
measured above.

## First attempted fix: graph-derived witnesses on the whole graph — insufficient

The first fix attempted (`fundamental_cycle_structure` in `scripts/mixer.py`)
derives a witness candidate directly from the graph's fundamental-cycle
structure (an `O(n_nodes)` tree-path computation per candidate, not a
combinatorial search) — the union of `e` and `f`'s own fundamental cycles,
verified against the real spanning-tree set rather than assumed. This is
correct (confirmed on small synthetic instances: still leak-free, still
fully connected, just with larger, non-minimal witnesses than the tuned
brute-force cap finds) — but checking its witness-size distribution
directly on the whole real IEEE 33-bus graph, *before* spending hours
building and verifying it, showed it doesn't solve the practical problem:
across the 630 real candidates, the local structural witness ranges from
5 to 34 qubits, **mean 22.3, median 23**. The reason: the 33-bus feeder's
5 tie switches don't just have individually long cycles (6-20 edges, see
above) — their cycles heavily *overlap*, so nearly the entire 36-edge
non-bridge core of the graph is entangled together. There is no
meaningfully smaller graph-derived witness available on the whole graph;
a 20-34-qubit controlled gate is impractical to synthesize, and exact
per-term leakage verification (`verify_term_no_leakage`) itself becomes
infeasible above roughly `k≈20-25` (a `2^k × 2^k` unitary). This was
caught by directly computing the witness-size distribution before
committing to the expensive run, not discovered after burning hours of
compute.

## The actual fix: zone decomposition

A standard partition/assembly strategy: partition the graph into zones
via a min tie-line-cut criterion, solve a matroid-basis mixer
independently per zone (small qubit count each), plus one small
"assembly" mixer on the graph formed by contracting each zone to a
supernode. A standard graphic-matroid contraction/deletion property
guarantees the union of every zone's spanning tree plus the assembly
problem's spanning tree is a spanning tree of the whole graph — the same
mixer construction applies unchanged at both scales, only the edge set
changes.

![A graph partitioned into zones, plus the contracted assembly problem](../results/illustration_decomposition.png)

Implemented in `scripts/zone_decomposition.py` and validated end-to-end
(construction, exact per-term leakage verification, transpiled gate
counts) in `scripts/run_zone_decomposition_validation.py`, writing
`results/zone_decomposition_results.csv`. Tested at 3, 4, 5, and 6 target
zones on the real IEEE 33-bus graph:

| target zones | max witness size (zones + assembly) | all subproblems leak-free |
|---|---|---|
| 3 | 4 | yes |
| 4 | 1 | yes |
| 5 | 2 | yes |
| 6 | 4 | yes |

Every configuration cleared the target of witness size < 10 by a wide
margin, and **every zone and assembly circuit was verified exactly
leak-free** (`verify_all_terms_no_leakage`), not just classically
connected. `target_zones=4` is the cleanest result: every zone turns out
to be internally a single tree (no internal ties at all — 0 spanning-tree
choices, 0 terms, nothing to condition on), and the one small assembly
problem (4 supernodes, 8 boundary edges, 24 spanning trees) needs only
`max_witness_size=1`, with a modest 132 CX / depth 272 transpiled
circuit.

**Two things this does NOT paper over**:

1. **A naive partitioner can silently break the decomposition's own
   precondition.** An earlier version of `zone_decomposition.py` used
   plain Kernighan-Lin bisection (a min-edge-cut heuristic) and produced
   a zone whose *internal* edges didn't connect all its own nodes (some
   zone nodes were only reachable via boundary/tie edges) in 3 of 4
   tested configurations — silently violating the decomposition
   proposition's precondition that each zone have its own spanning tree.
   Caught by checking `enumerate_spanning_trees` returned a non-empty set
   for every zone, not assumed. Fixed by `_split_connected`
   (`zone_decomposition.py`): recursively removing one *spanning-tree*
   edge from each piece, which guarantees both resulting pieces are
   internally connected by construction, tie-broken by cut size and
   balance.
2. **Small witness size does not mean a small circuit.** The
   `target_zones=3` configuration has one non-trivial zone (16 qubits, 10
   terms, `max_witness_size=4`) whose transpiled circuit is 7,904 CX
   gates and depth 15,785 — large for 16 qubits, because a term with
   several valid witness patterns compiles to one controlled gate *per
   pattern* (see "What the fixed circuit does" below), and generic
   multi-controlled-gate synthesis at `optimization_level=1` is not
   gate-count-optimal. `target_zones=5`'s assembly problem shows the same
   effect at witness size 2 (1,572 CX, depth 2,991). This is a real,
   measured cost of witness-conditioning that should be reported
   alongside the witness-size numbers, not just the witness sizes alone
   — see `results/zone_decomposition_results.csv` for every configuration's
   numbers.

## Does decomposition keep witness size bounded as the network gets larger?

Everything above validates decomposition on ONE real topology (33-bus) at
a few partition granularities — it doesn't by itself show the fix keeps
working as real-scale networks grow arbitrarily large, since only one
real instance size is available to test against. `scripts/
generate_feeder_graph_long_range_ties` (in `graphs.py`) reproduces the
specific failure mode found on real data — tie edges chosen to *maximize*
tree-path length, rather than the original synthetic model's *nearest*-
neighbor ties — at any controllable network size, so this can be tested
as a trend rather than a single data point.

`scripts/run_decomposition_scaling_study.py` sweeps `n_nodes` from 10 to
120 (12 seeds each, `k_ties=5` fixed — matching the real feeder's actual
tie count), and for each size compares two things: the whole graph's
largest fundamental-cycle-membership set (a cheap, direct proxy for the
witness size a whole-graph construction would need — this is exactly the
quantity that was 22-34 qubits on the real 33-bus graph), against the
*actual, exactly-verified* max witness size across every zone and
assembly subproblem when decomposed at a **fixed target zone size of 8**
nodes (so zone count grows with the network, but each zone's own qubit
budget does not). Results (`results/decomposition_scaling_results.csv`,
`results/decomposition_scaling_summary.csv`,
`results/decomposition_scaling_plot.png`):

| n_nodes | whole-graph naive proxy (mean) | decomposed max witness (mean) |
|---|---|---|
| 10 | 11.7 | 1.3 |
| 20 | 18.2 | 3.5 |
| 40 | 27.5 | 2.7 |
| 80 | 44.1 | 2.7 |
| 120 | 58.2 | 2.6 |

The whole-graph proxy grows roughly linearly with network size (12→58
over a 12x size increase); the decomposed max witness size **stays flat,
hovering around 2.6-3.5 from `n_nodes=20` onward, with no growth trend**
— every zone and assembly subproblem at every size and seed was exactly
verified leak-free (`all_leak_free=True` throughout). This is direct
empirical support, not just the structural argument, for "decomposition
scales to arbitrarily large real-like networks": per-subproblem qubit and
witness-size budget can be held constant by choosing a fixed zone-size
target, independent of total network size.

**Two honest wrinkles in this result, not smoothed over**:

1. At `n_nodes=10` (below the target zone size of 8), the partitioner
   forms only 1 zone — i.e. no real decomposition happens, and that
   single "zone" is the whole graph. Several of those seeds show
   `all_connected=False` (some candidates dropped) — the whole-graph
   failure mode reproducing itself even at this small size, a useful
   corroborating signal (confirms the synthetic long-range-tie model
   faithfully reproduces the real problem) but meaning the `n_nodes=10`
   row isn't a real test of decomposition, only of the undecomposed
   baseline at small scale.
2. Decomposition does not *always* reach full connectivity at this
   witness cap, even once it's actually being applied: `n_nodes=80` and
   `n_nodes=120` each had 1 of 12 seeds land on `all_connected=False` too
   — a small residual fraction (~8% of the non-trivial cases tested)
   where some candidate still exhausted the witness-size-4 search within
   a zone or the assembly problem. The overwhelming majority of
   instances reach full connectivity with small witness, but
   "decomposition always fully connects" is not accurate — it's
   "decomposition fully connects in the large majority of instances
   tested, occasionally not." See `results/decomposition_scaling_results.csv`
   for exactly which (`n_nodes`, `seed`) pairs this affects.

## What drives circuit cost, and the logic-minimization fix

Witness size staying flat does not mean transpiled circuit cost stays
flat. Investigated with the same 12-seed sample above and a direct
breakdown of what drives it, rather than just noting the observation:

- **`optimization_level=3` does not meaningfully help.** Rerunning every
  instance at both `optimization_level=1` (this repo's usual setting) and
  `3`: mean max-CX changes by roughly 1-3% between them (e.g. `n_nodes=20`:
  12,815 vs 12,783 mean, before the fix below) — a generic transpiler
  optimization pass doesn't restructure *how many separate controlled
  gates* a term compiles to, which is the actual driver (next point), so
  it was never going to fix this.
- **Mean CX does NOT show a growth trend with network size** — 3,324 →
  12,815 → 2,351 → 2,308 → 1,912 (mean max-CX at `n_nodes` 10, 20, 40, 80,
  120, before the fix below) — noisy and spiky (particularly at
  `n_nodes=20`, driven by a few seeds landing several terms at the
  witness cap with many valid patterns each) but not trending upward.
  Tracking the actual driver directly (mean max valid-patterns-per-term:
  4.0 → 5.4 → 2.4 → 2.3 → 2.3) shows the same thing. Circuit cost looks
  driven by occasional "unlucky" individual instances, not a
  size-dependent trend — but individual outliers still reach genuinely
  large absolute numbers (up to 36,588 CX in one `n_nodes=20` seed),
  which matters in practice regardless of whether it trends with size.
- **The fix: logic-minimized witness patterns, implemented and
  measured.** `_minimize_patterns` (`scripts/mixer.py`) runs
  Quine–McCluskey minimization on each term's `valid_patterns`, producing
  an EXACT cover (same truth table, verified by rerunning
  `verify_correctness.py` and `run_zone_decomposition_validation.py`
  afterward — all instances still pass exactly) using fewer, often
  narrower-controlled gates — `exchange_term_circuit` now builds one
  controlled gate per minimized *cube* (which may leave some witness
  qubits as don't-cares) instead of one per exact pattern. Rerunning the
  full `n_nodes` 10-120 sweep with this fix, same 12 seeds/size, same
  instances:

  | n_nodes | mean max-CX, before | mean max-CX, after | max max-CX, before | max max-CX, after |
  |---|---|---|---|---|
  | 10 | 3,324 | **513** (6.5x) | 19,204 | **2,404** (8.0x) |
  | 20 | 12,815 | **4,953** (2.6x) | 36,588 | **11,790** (3.1x) |
  | 40 | 2,351 | **1,597** (1.5x) | 9,216 | **4,450** (2.1x) |
  | 80 | 2,308 | **2,057** (1.1x) | 7,174 | **4,932** (1.5x) |
  | 120 | 1,912 | **1,515** (1.3x) | 10,230 | **5,464** (1.9x) |

  Improvement in every single size bucket, ranging from modest (~10-25%
  at larger sizes, where patterns per term were already close to minimal)
  to dramatic (6-8x at `n_nodes=10`, where several instances had terms
  needing many valid patterns). No case got worse. `all_leak_free=True`
  and the exact same set of `all_connected` outcomes as before
  minimization (expected: this only changes circuit synthesis, not
  witness search or term selection) — full before/after data in
  `results/decomposition_scaling_results_before_minimization.csv` vs.
  `results/decomposition_scaling_results.csv`. The real IEEE 33-bus
  decomposition result (`results/zone_decomposition_results.csv`) shows
  the same pattern: `target_zones=5`'s assembly problem dropped from
  1,572 to 552 CX (2.85x), `target_zones=6`'s from 7,216 to 3,686 CX
  (1.96x); `target_zones=3`'s hardest zone was unchanged at 7,904 CX —
  minimization only helps where patterns actually share don't-care
  structure, and not every instance has any to exploit. **Still driven by
  occasional individual instances, not network size**: even after
  minimization, mean CX doesn't show a growth trend with `n_nodes` (513 →
  4,953 → 1,597 → 2,057 → 1,515) — the `n_nodes=20` spike is the same
  handful of "unlucky" seeds as before, just proportionally smaller now.

## What the fixed circuit does

For each selected exchange `(e,f)` with witness qubits `w_1,...,w_k`
(`k` in `{0,1,2}` on the synthetic sweep, up to 4 on the zone-decomposed
real-feeder result above — see that section for why the whole real graph
needs much more, and why decomposition is what keeps `k` small there),
the mixer applies the `RXX+RYY` swap block on `e,f`, controlled on a
Quine-McCluskey-minimized cover of the witness-qubit patterns under which
the exchange is valid (`_minimize_patterns`, `scripts/mixer.py` —
determined by direct lookup against the enumerated spanning-tree set, not
assumed; see the "does decomposition keep witness size bounded" section
above for why this minimization step exists and what it measurably saves).
Each minimized cube may leave some witness qubits uncontrolled
(don't-care), so a single exchange can compile to controlled gates
narrower than its full witness size `k`, not just fewer of them. When
`k=0`, this is exactly the unconditional gate -- still used when it's
actually safe, not conditioned unnecessarily.

## Verification

`scripts/verify_correctness.py` builds the exact `2^n_qubits x 2^n_qubits`
circuit unitary (`qiskit.quantum_info.Operator`) for a handful of small
instances and checks that every spanning tree's column has *zero*
amplitude (not "small", exactly zero within floating-point tolerance) on
every non-spanning-tree basis state. This is Hadfield's constraint-
preservation condition, checked directly against the actual compiled
circuit, not argued from the construction's intent. All instances checked
pass exactly. This check is only feasible for small `n_qubits` (exact
statevector simulation); the scaling measurements in `/results` go well
beyond what this exact check can reach, and rely on the witness-search
algorithm being correct at those sizes, not on having re-verified every
individual instance's unitary.
