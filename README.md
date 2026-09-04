# Matroid basis-exchange QAOA mixer: scaling study

## The question

For the constraint "the selected edges form a spanning tree of the graph"
(a basis of its graphic matroid), does a QAOA mixer that preserves that
constraint exactly — rather than enforcing it with a penalty term in the
cost function — compile into a circuit whose gate count and depth scale
favorably with graph size, for sparse graphs resembling real distribution
feeders? And if the direct answer turns out to be "not on real topology,
not without help": can something be built that still solves the problem
exactly, cheaply enough to matter on both fault-tolerant and near-term
(NISQ) hardware?

This is narrower than general worst-case graphs, where basis-exchange
move sets can blow up exponentially; it is the graph class relevant to
feeder reconfiguration, where the network is sparse (a radial backbone
plus a small number of tie switches) and close to planar. A
constraint-preserving mixer restricts QAOA's search to the feasible
subspace directly, instead of relying on a penalty coefficient in the
cost Hamiltonian to discourage infeasible states — trading a larger,
structured mixer circuit for a smaller effective search space and no
penalty-weight tuning. This repo answers whether that trade is cheap
enough, empirically, to be worth making for this problem class.

**Short answer**: not for free. The direct construction works cleanly on
a synthetic sparse-feeder family, but fails outright on real, published
topology (the IEEE 33-bus feeder) — genuinely infeasible witness sizes,
not a tunable search limit. Two things fix it: **zone decomposition**
(splits the constraint into small, exactly-solvable pieces, provably
lossless) and a **cost-aware bounded-witness mixer** (bounds circuit cost
directly, at a small, measured leakage cost). Combined, and refined
further under an escalating ladder of increasingly realistic assumptions,
they produce a construction that is exact and NISQ-plausible on real
feeder topology — validated directly on two real networks, not just
synthetic proxies.

## The strategy, at a glance

1. **Establish the baseline** on synthetic data: does the exact
   construction scale at all, under the simplest realistic assumptions?
   (Results, part 1, below.)
2. **Stress-test it against real topology.** It breaks — trace *why*
   before reaching for a fix.
3. **Fix it two independent ways**, each solving a different piece of the
   problem: zone decomposition (splits the *constraint*, exact) and a
   bounded-witness mixer (bounds *cost* directly, approximate but
   measured). Combine them.
4. **Stress-test the fix** under an escalating ladder of more realistic
   assumptions (tie placement, tie-count growth) calibrated to, and later
   checked against, real feeder data — including a real correctness bug
   and a dead end found and caught along the way, not smoothed over.
5. **Push cost down further** once the fix works, specifically to reach
   NISQ-plausible gate counts, not just "scales better than before."
6. **Validate the whole thing directly on real networks** — not only the
   synthetic families every earlier step was calibrated to.

Each numbered step below links to the full technical account in `docs/`;
this document is the narrative and the numbers that matter for deciding
whether to trust the result, not the complete derivation.

## Background: what "feasible" means here, and why some ties are worse than others

A distribution feeder has more switchable connections than it strictly
needs to reach every load: a normally-closed backbone plus a small number
of normally-open **tie** switches, held open in everyday operation and
closed only temporarily (e.g. to reroute power around a fault). The
network must stay **radial** — every node reached, no loops — for
protection and fault-isolation reasons.

![A feasible configuration (spanning tree) vs. an infeasible one (a loop)](results/illustration_feeder_problem.png)

Encoded directly: qubit `i` = 1 iff switch `i` is closed. A configuration
is feasible iff its closed-switch set is a spanning tree of the network
graph — exactly a basis of the network's graphic matroid, shown above on
a small example graph (14 nodes, one candidate tie switch, reused as the
running example in the figures below).

A QAOA mixer's job is to move probability between states without ever
leaving the feasible set. Concretely here, that means one operation
repeated many times: **a basis exchange** — close one switch, open
another, land on a different feasible tree.

![One basis-exchange move: before, the move itself, and after](results/illustration_basis_exchange_move.png)

Matroid theory guarantees moves like this, applied one at a time, are
enough to reach any feasible configuration from any other — so the
*entire* mixer is built from a set of these exchanges, one small circuit
term per move. That would be a trivial circuit-design problem if any
switch pair could be exchanged freely — but it can't: closing a switch
without also opening the right one produces a loop, not another tree.
Which switch has to open depends on the rest of the network's *current*
configuration, not just the two switches being touched, so each
exchange's circuit term generally has to be **conditioned** on other
qubits to fire only when the move is actually valid.

The qubits an exchange's term depends on are exactly the *other* tree
edges on its tie switch's **fundamental cycle** — call this set of qubits
the exchange's **witness**. A short cycle means a small witness; a long
one means a large witness, and a large witness means an expensive,
deeply-controlled gate.

![A tie edge's fundamental cycle, short vs long](results/illustration_fundamental_cycle.png)

*(synthetic — illustrative, same running example graph, tie switch choice
changed from nearest-neighbor to long-range, `scripts/plot_illustrations.py`)*
2 witness qubits on the left, 8 on the right, in this small example. Real
tie switches are deliberately placed to link *distant* parts of a network
for redundancy — closer to the right panel than the left — which is
exactly what makes this problem harder on real topology than it looks on
a graph where ties happen to be short-range.

**This is the actual question this repo measures**: how much does that
conditioning cost, in gates and circuit depth, and does that cost grow or
stay manageable as the network gets larger and as tie placement gets more
realistic? A mixer whose per-move conditioning cost explodes is a mixer
that doesn't scale, independent of anything else about QAOA.

**Why this matters beyond this repo**: keeping that cost bounded is what
would let QAOA run at network sizes where quantum computing could
eventually offer an advantage over classical methods on constrained
problems like this one. Classical formulations of the same radiality
requirement still have to encode it explicitly — as penalty terms in an
objective, or as constraints a solver enforces — while a
constraint-preserving mixer builds feasibility directly into the
algorithm's dynamics instead. That's the promise a bounded-cost result is
a *precondition* for, not proof of — see "What this does and doesn't
demonstrate" below for the boundary between the two.

## Techniques used, and what problem each one solves

In the order they're needed to understand the results, not the order
they were built (the build order — including two real detours — is in
"Headline result 4" of `docs/scaling-ladder-and-decomposition.md`, kept
there as part of the evidence, not hidden here).

**1. Exact matroid mixer (`mixer.py`, `build_matroid_mixer`).** For each
candidate exchange, brute-force search for the smallest witness that
makes the term's firing rule exactly correct (no leakage outside the
feasible set, verified against the true validity function, not assumed).
Cheap and exact when witnesses stay small — which they do on a synthetic
family where ties are short-range by construction, but not in general:
**tie range** (long fundamental cycles) and, independently, **tie
density** (many ties even at short range) both push true minimal witness
sizes past any practical search cap, on real topology and on a
density-focused synthetic family respectively (`docs/bounded-witness-mixer.md`'s
Finding 1). When that happens, the exact construction doesn't get
gradually more expensive — it gives up on candidates outright
(`dropped_candidates` climbs to 51-95%) and the resulting mixer stops
being fully connected.

**2. Zone decomposition (`zone_decomposition.py`) — fixes the *range*
failure.** Don't build one circuit for the whole network: partition into
zones via a min tie-line-cut, solve each zone's matroid mixer
independently (small qubit count each), plus one small assembly mixer
over the contracted zone graph. Guaranteed **exact** by graphic-matroid
contraction/deletion (the union of every zone's spanning tree plus the
assembly problem's spanning tree is provably a spanning tree of the whole
graph) — this is a structural fix to the *constraint*, not an
approximation.

![A graph partitioned into zones, plus the contracted assembly problem](results/illustration_decomposition.png)

*(synthetic — illustrative example, 24 nodes)* it works specifically
*because* the range failure is local to individual long-range ties: keep
each zone small enough that any tie edge still inside it stays
short-range-like, and push the genuinely long-range connections out to
one small assembly problem instead of forcing the whole graph to absorb
their cost.

**3. Bounded-witness mixer (`truncated_mixer.py`) — bounds *cost*
directly, at a measured leakage cost.** An alternative (or complement) to
decomposition for cases where a fixed witness cap is preferred to
splitting the problem into zones: instead of dropping a candidate
exchange outright when no small *exact* witness exists, search for a
witness of a fixed, capped size using a **majority-vote** validity rule,
and accept the resulting leakage — provided it's measured, not assumed,
and stays small. `leakage_trace.py` traces the actual probability mass
this costs on real trajectories (exactly, not sampled), separately from
the abstract leakage rate the search itself optimizes.

**4. Cost-aware search (`cost_alpha`) — fixes a real 38x inefficiency.**
The first version of (3) picked witnesses by leakage alone, and since
adding a witness qubit can only ever weakly *reduce* leakage, the search
had every incentive to walk to the cap every time — 17,386–32,520 CX
gates measured at just 15-17 qubits, more than a real 37-qubit
decomposed feeder circuit. `cost_alpha` adds a cost penalty to the search
objective directly, closing most of that gap; it is this construction's
validated default (`docs/bounded-witness-mixer.md`). A later
per-term-adaptive version of this same idea looked like a further
improvement on an initial sweep, and was **not** once properly
re-validated — a real dead end, kept in the record
(`docs/scaling-ladder-and-decomposition.md` §5) rather than quietly
dropped.

**5. Hierarchical (density-aware) decomposition** — for cases where the
*assembly* graph from step 2 is itself still large or dense: recurse zone
decomposition on it, scaling `target_zone_size` inversely to the current
graph's measured density rather than a fixed schedule (found, by direct
comparison, to beat both a fixed recursive schedule and naive
zone-shrinking — `docs/scaling-ladder-and-decomposition.md` §8).

**6. Cost-capped decomposition — guarantees a cost target instead of
hoping for one.** Steps 2 and 5 both pick a zone size up front and hope
the resulting cost is acceptable; it usually is, but decomposed cost
turns out to have much higher seed-to-seed variance than whole-graph cost
(coefficient of variation up to 99% at some sizes — decomposition trades
one big averaging problem for many small independent ones, and a single
"unlucky" zone can dominate a seed's total). The fix: build each
subproblem, transpile it, check its **actual** cost against a threshold,
and recursively re-split anything over it — measure, don't guess. Gets
**every seed of the synthetic ladder under 500 CX**, and comes within
1.6% of that target on real data (`docs/scaling-ladder-and-decomposition.md` §11).

## Why the ladder is shaped this way, and why log growth

Validating steps 2-6 once, on one real anchor point, isn't enough to
trust the result generally — so this repo stress-tests them across an
**escalating ladder** of increasingly realistic assumptions, varied along
two axes, both independently calibrated to the one real anchor point this
repo started with (5 ties at the 33-bus feeder):

- **tie placement**: short-range (nearest-tie, headline-result-1's
  family) vs. long-range (the family the real-feeder failure and the
  bounded-witness safety survey are built on).
- **tie-count growth**: log (`round(1.43*ln(n))`, mild) vs. linear
  (`round(0.1515*n)`, aggressive).

**This choice was checked against real data, not just assumed.** Three
real/real-benchmark networks (CIGRE MV, IEEE 33-bus, and a real German MV
network, `mv_oberrhein`, via `pandapower.networks`), spanning a 12x size
range:

| network | n_bus | n_ties | ties/n_bus |
|---|---|---|---|
| CIGRE MV | 15 | 3 | 0.200 |
| case33bw (IEEE33) | 33 | 5 | 0.152 |
| mv_oberrhein | 179 | 6 | 0.034 |

The ratio drops 6x from smallest to largest — flatly inconsistent with
linear growth (which would keep it roughly constant), reasonably
consistent with log growth (`ties/log(n_bus)` comes out to 1.11, 1.43,
1.16 — much tighter). Real tie *range* was checked the same way: real
ties consistently span 33-45% of network diameter, matching this repo's
long-range generator far better than its short-range one (which gets
*worse*, not better, with scale). **`CONDITIONS` therefore defaults to
the two log-growth conditions**; linear growth is kept as an explicit,
fully reproducible stress test, not deleted, but no longer presented as
equally realistic (`docs/scaling-ladder-and-decomposition.md` §9 has the
full check, including honest caveats — 3 points spanning one order of
magnitude rules out linear but doesn't fit an actual growth law).

## Results, part 1: the exact construction (what a fault-tolerant device could run today)

The exact construction (technique 1) is the one with no approximation and
no leakage risk from splitting the problem — the number relevant to
"does this work at all, with unlimited circuit budget."

**Synthetic baseline.** Built directly on a synthetic sparse-feeder
family, `k_ties` held fixed as the network grows:

![Mixer circuit CX count and depth vs. instance size](results/scaling_plot.png)

Measured across 26 instance sizes (`n_qubits` 10–35), 5 seeds per size.
**Within this range, cost trends flat-to-*decreasing* with size** — the
network gets sparser as it grows (same few ties spread across more
nodes), so there are fewer, cheaper exchanges to condition. This trend is
a property of the flat-tie-count assumption specifically, not the
construction in general: under log-scaled tie-count growth (the
realistic model, per above), the exact construction's cost instead
**increases** — 3.4x more on average, 5.8x more at the 33-node
calibration point (`docs/scaling-ladder-and-decomposition.md` §1).

**Real topology (IEEE 33-bus feeder).** Built the same way, directly on
the real, published feeder (Baran & Wu 1989, `n_nodes=33`, `n_qubits=37`,
via `pandapower.networks.case33bw`): **it fails.** 91% of candidate
exchanges have no witness of any practical size, and the resulting mixer
is not fully connected — a completeness gap, not a performance number
(`results/real_feeder_results.csv`). Not a search-cap artifact either:
the witness the whole graph would actually need spans dozens of qubits,
measured directly.

**Zone decomposition recovers exactness, and stays cheap as size grows.**
On the same real 33-bus graph, every zone/assembly subproblem is exactly
leak-free with witness size 0-4, vs. 22-34 on the whole graph
(`results/zone_decomposition_results.csv`). On a synthetic family
reproducing the real failure mode at controllable sizes (12 seeds/size,
`n_nodes` 10-120):

![Whole-graph vs. zone-decomposed witness requirement](results/decomposition_scaling_plot.png)

The whole-graph requirement grows roughly linearly; the decomposed
requirement (fixed zone size) stays flat around 3 — the central claim of
the decomposition fix. Re-plotted by qubit count for direct comparison
against the synthetic baseline above:

![Whole-graph vs. zone-decomposed witness requirement, by qubit count](results/decomposition_scaling_by_qubits_plot.png)

**On a second real network (CIGRE MV, 15 buses, 3 ties)**: exact
whole-graph construction succeeds here (small enough graph), at 16 terms,
12,220 CX, 24,245 depth — decomposition still cuts that 3.3x, to 3,668 CX
(`docs/scaling-ladder-and-decomposition.md` §10).

**Full account**, including the failure's root cause, a first attempted
fix that was insufficient, and a follow-up circuit-cost investigation:
`docs/circuit-validity.md`. `docs/mixer-construction.md` is the
standalone technical reference (matroid theory, exact circuit derivation,
verification arguments).

## Results, part 2: the NISQ-ready construction (decomposition + cost-aware bounded-witness mixer)

Exactness alone doesn't make a circuit runnable soon: rough NISQ
feasibility arithmetic (`fidelity ≈ (1-p)^N_CX`, published two-qubit gate
fidelity ranges) shows even the flat-decomposed real-network numbers
above are well outside a plausible regime once tie count and range are
modeled realistically:

| CX count | best-case trapped-ion (p=0.001) | typical superconducting (p=0.005) |
|---|---|---|
| 100 | 90% | 61% |
| 500 | 61% | 8% |
| 1,000 | 37% | 0.7% |
| 5,000 | 0.7% | ~0 |
| 13,000 | ~10⁻⁶ | ~0 |

**The escalating ladder result, once made cost-aware (techniques 3+4)**:
cost roughly plateaus by `n_nodes=60` rather than growing further; tie
*placement* (short vs. long range) drives the cost level far more than
tie-count *growth rate* (log vs. linear) does; the mixer stays fully
connected at every point tested
(`docs/scaling-ladder-and-decomposition.md` §2). Getting here also
surfaced a real bug — a search objective that can let most of a mixer's
terms silently stop firing at all, producing a circuit that *looks* safe
because it has stopped being a mixer — caught by tracking active vs.
inert witnesses explicitly, and a dead end (adaptive per-term cost
pressure) that looked like a further win before that fix and wasn't
after it. Both are covered in full, not smoothed over, in
`docs/scaling-ladder-and-decomposition.md` §4-5.

**Decomposition (technique 2) dominates almost everywhere it applies**
once combined with the cost-aware mixer: 5.6x-46.7x cheaper than the
whole-graph cost-aware construction, winning every single seed at every
size ≥ 30 nodes, for both log-growth conditions (§6). Hierarchical
decomposition (technique 5) pushes the hardest remaining condition
further: 1.46-1.59x cheaper than flat decomposition at real scale (§8).

**Cost-capped decomposition (technique 6) is the current recommendation**
— it doesn't just do well on average, it guarantees the target:

| condition | result |
|---|---|
| synthetic ladder, both log-growth conditions, all 5 sizes, all seeds | **100% under 500 CX**, mean feasible mass exactly 1.0 |

**Validated directly on real networks — not just synthetic proxies:**

| network | construction | CX | result |
|---|---|---|---|
| CIGRE MV (15 buses, 3 ties) | exact whole-graph | 12,220 | connected, exactly leak-free |
| CIGRE MV | decomposed | 3,668 | 3.3x cheaper |
| CIGRE MV | **cost-capped decomposed** | **274** | meets 500-CX target cleanly; mean feasible mass 0.985 |
| IEEE33 (33 buses, 5 ties) | exact whole-graph | 96 | **573/597 candidates dropped, disconnected — fails** |
| IEEE33 | decomposed | 132 | fully functional; mean feasible mass 1.0 |
| IEEE33 | **cost-capped decomposed** | **508** | 8 CX over target (1.6%); mean feasible mass 1.0 (perfect) |

IEEE33's small miss is not a mystery or a search failure: it's traced to
one small (4-node), genuinely irreducible dense multigraph core, where
neither cost pressure nor search-cap changes can help because those
candidates are already at the zero-leak witness width — there's no
cost/safety tradeoff left to spend. Closing it would need a smarter zone-
*choice* strategy (which nodes get grouped together, not just how many),
identified but not attempted (`docs/scaling-ladder-and-decomposition.md` §11).

**Full account**, including both real methodology mistakes made and
caught along the way (worth reading for what that looked like, not just
the corrected numbers): `docs/scaling-ladder-and-decomposition.md`.
`docs/bounded-witness-mixer.md` covers the density-failure axis, the
bounded-witness construction itself, and the cost-aware search fix in
full detail.

## What this does and doesn't demonstrate

This repo shows the mixer **construction** is correct and scales — both
exactly (fault-tolerant-relevant) and, via decomposition and cost-capping,
cheaply enough for a plausible NISQ target (real-topology-relevant). It
does **not** show a quantum algorithm outperforming a classical baseline:
there is no cost-Hamiltonian/oracle integration, no QAOA execution, no
classical solver comparison, and no test of the iterative boundary
coupling (an ADMM-style loop) this kind of decomposition would need for
the actual optimization objective — only the radiality constraint,
tested here, decomposes exactly. Treat this as evidence the algorithm is
buildable and scalable, a precondition for an advantage claim, not the
claim itself. See `methodology.md` for the precise boundary of what was
measured.

## A correctness bug found and fixed during this study

An early version of the mixer circuit was **wrong**: an unconditional
swap gate leaked probability outside the feasible subspace on roughly
half the instances tested, caught by direct verification against the
exact circuit unitary rather than assumed safe. The fix — witness-qubit
conditioning, found by direct search — is in `docs/circuit-validity.md`.
`scripts/verify_correctness.py` re-runs this check independently of the
scaling measurements.

## How to reproduce

```
python -m venv .venv
source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt

python scripts/verify_correctness.py                  # exact correctness check (small instances)
python scripts/run_scaling_study.py                    # synthetic sweep -> results/scaling_*.{csv,png}
python scripts/run_real_feeder_validation.py            # real IEEE 33-bus graph, whole-graph construction (fails)
python scripts/run_zone_decomposition_validation.py     # zone-decomposition fix, real IEEE 33-bus graph
python scripts/run_decomposition_scaling_study.py       # does the fix hold as network size grows? (synthetic)
python scripts/plot_decomposition_by_qubits.py          # re-plots the above by qubit count, not node count
python scripts/plot_illustrations.py                    # explanatory diagrams (not measurements)
python scripts/verify_leakage_trace.py                  # correctness check for the bounded-witness mixer's tooling
python scripts/tie_density_sweep.py                     # a second witness-blowup axis: tie density (see docs/bounded-witness-mixer.md)
python scripts/run_bounded_witness_safety_survey.py     # real-scale safety + cost survey for the bounded-witness mixer
python scripts/measure_truncated_mixer.py               # circuit-cost investigation: truncated mixer's own cost trend + vs. exact
python scripts/truncated_witness_cap_sweep.py           # witness-cap vs. cost/safety tradeoff (density family)
python scripts/truncated_witness_cap_sweep_longrange.py # same tradeoff, generalized to the long-range family + real scale
python scripts/truncated_mixer_search_refinement.py     # cap=0/1 + the cost-aware search fix (cost_alpha)
python scripts/run_scaling_study_log_ties.py            # is headline result 1's flat-tie-count assumption load-bearing?
python scripts/run_cost_aware_scaling_ladder.py         # escalating realism ladder for the bounded-witness mixer
python scripts/run_cost_aware_scaling_ladder_aggressive.py  # pushing cost down: the never-fire bug, made visible
python scripts/run_fixed_alpha_ladder.py                # fixed cost_alpha baseline, re-measured post-fix
python scripts/run_decomposed_cost_aware_ladder.py      # decomposition + bounded-witness mixer combined (the actual fix)
python scripts/run_best_of_both_ladder.py               # confirms decomposition dominates -- never loses a seed
python scripts/run_hierarchical_decomposed_ladder.py    # density-aware recursive decomposition, for the still-hard conditions
python scripts/run_real_networks_hierarchical.py        # this branch's actual construction, run directly on CIGRE MV + IEEE33
python scripts/run_cost_capped_decomposition.py         # decomposition that MEASURES cost and re-splits over-threshold subproblems
```

All scripts are deterministic (fixed random seeds); re-running should
reproduce the committed files in `results/` exactly, modulo `qiskit`/
`networkx` version differences in transpilation.

## Scope

This repo contains the measurement methodology and results for the
question above only. It does **not** include a production
mixer-compilation library, does not cover other constraint classes, and
does not include hardware execution, QAOA solution quality, or a
classical baseline comparison — see `methodology.md` for the precise
boundary of what was measured.

## Repository layout

- `methodology.md` — graph generation, move generation, and exactly what
  "gate count" and "depth" mean.
- `docs/mixer-construction.md` — technical reference: matroid theory,
  move generation, exact circuit-level derivation, witness conditioning
  and its logic-minimization, correctness-verification arguments.
- `docs/circuit-validity.md` — the full narrative: the correctness bug,
  the real-topology failure and its root cause, the decomposition fix and
  its scaling validation, and the circuit-cost investigation.
- `docs/bounded-witness-mixer.md` — a second, density-driven witness-blowup
  axis, the bounded-witness mixer as an alternative to decomposition
  (construction, a correctness bug the cross-check caught in shared code,
  real-scale safety survey), and the circuit-cost investigation that
  found the first version of that construction cost up to 38x more than
  necessary, plus the cost-aware search fix that mostly closed the gap.
- `docs/scaling-ladder-and-decomposition.md` — does the flat-tie-count
  assumption hold under a more realistic growth model (no); an
  escalating realism ladder for the bounded-witness mixer and a NISQ
  hardware feasibility check; a real bug (silent term collapse) found
  while pushing cost down further; a dead end (per-term adaptive cost
  pressure) that looked like a fix and wasn't, caught by re-validation;
  the fix that actually works (decomposition + the bounded-witness
  mixer, dominating almost everywhere it applies) -- including a SECOND
  instance of the same adaptive-search bug found inside the
  decomposition script itself; density-aware hierarchical decomposition,
  pushing the hardest remaining conditions further toward NISQ
  feasibility; a real-topology check validating log-growth/long-range
  tie modeling against real (and real-benchmark) feeder data, after
  which linear tie-count growth moved from a default condition to an
  explicit stress test; direct validation of the actual construction on
  two real networks; cost-capped decomposition (measures each
  subproblem's actual cost and recursively re-splits anything over a
  threshold, guaranteeing every synthetic-ladder seed stays under 500
  CX, and coming within 1.6% of that target on real data); and two
  methodology mistakes this investigation made and caught (comparing
  results across a code fix without re-running both sides, twice).
- `scripts/` — synthetic sweep: `graphs.py`, `mixer.py`, `measure.py`,
  `run_scaling_study.py`, `plot.py`, `verify_correctness.py`. Real
  topology: `real_feeders.py`, `run_real_feeder_validation.py`,
  `investigate_fundamental_cycles.py` (root-cause diagnostic).
  Decomposition fix: `zone_decomposition.py`,
  `run_zone_decomposition_validation.py`,
  `run_decomposition_scaling_study.py` (does it hold as size grows?);
  `plot_decomposition_by_qubits.py` re-plots that sweep's results by
  qubit count instead of node count, for direct comparison against
  headline result 1 (reads the existing CSV, no re-measurement).
  `plot_illustrations.py` generates the explanatory diagrams (not
  measurements). Bounded-witness mixer: `random_trees.py` (Wilson's
  algorithm + exchange-graph-walk sampling), `truncated_mixer.py`
  (bounded-witness construction, including the cost-aware search --
  `cost_alpha` -- described below), `leakage_trace.py` (exact + sparse
  danger-mass tracing), `verify_leakage_trace.py` (correctness check for
  all three), `tie_density_sweep.py` (the density-driven witness-blowup
  axis), `run_bounded_witness_safety_survey.py` (real-scale safety +
  circuit-cost survey). Circuit-cost investigation:
  `measure_truncated_mixer.py` (own cost trend + head-to-head vs. the
  exact construction on the density family), `truncated_witness_cap_sweep.py`
  / `truncated_witness_cap_sweep_longrange.py` (witness-cap vs. cost/safety,
  density family and long-range family + real scale), and
  `truncated_mixer_search_refinement.py` (cap=0/1 instability check, and
  the `cost_alpha` sweep that set the construction's current default).
  Escalating realism ladder + decomposition:
  `run_scaling_study_log_ties.py` (the flat-tie-count assumption,
  stress-tested), `run_cost_aware_scaling_ladder.py` (the four-condition
  ladder), `run_cost_aware_scaling_ladder_aggressive.py` (the never-fire
  collapse, deliberately reproduced with `active_terms` tracked so it's
  visible in the CSV), `run_cost_aware_scaling_ladder_alpha_sweep.py`
  and `truncated_mixer_search_refinement.py`'s adaptive-search addition
  in `truncated_mixer.py` (the adaptive-alpha dead end),
  `run_fixed_alpha_ladder.py` (the fixed-`cost_alpha` baseline,
  re-measured after the never-fire fix so the final comparison is valid),
  `run_decomposed_cost_aware_ladder.py` (the fix that works), and
  `run_best_of_both_ladder.py` (confirms decomposition never loses once
  it can be applied at all). `run_hierarchical_decomposed_ladder.py`
  (density-aware recursive decomposition for the still-hard conditions,
  three iterations each empirically validated),
  `run_real_networks_hierarchical.py` (this branch's actual construction
  run directly on `real_feeders.load_cigre_mv` -- new -- and the
  existing `load_ieee33`), and `run_cost_capped_decomposition.py` (measures
  each subproblem's actual CX cost and recursively re-splits anything over
  a threshold, rather than picking a zone size and hoping -- 100% success
  on the synthetic ladder, tested on both real networks too) round it out.
- `results/` — one CSV/plot pair per script above, all generated, none
  hand-edited; `*_before_minimization.*` files are the pre-fix numbers,
  kept for the before/after comparison in `docs/circuit-validity.md`;
  `illustration_*.png` are the explanatory diagrams
  (`illustration_feeder_problem.png`, `illustration_basis_exchange_move.png`,
  and `illustration_fundamental_cycle.png` share the exact same 14-node
  running example graph; `illustration_decomposition.png` uses a larger
  instance from the same random seed, needed for partitioning to be
  illustrative at all).

## License

Apache License 2.0 — see `LICENSE`.
