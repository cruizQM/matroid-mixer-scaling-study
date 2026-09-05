# Matroid basis-exchange QAOA mixer: scaling study

## The question

**The pitch, one level up**: QAOA's biggest practical obstacle on
constrained problems usually isn't the cost Hamiltonian — it's keeping
the search feasible at all, typically patched with penalty terms that
waste circuit depth fighting infeasibility instead of searching it, with
no guarantee against landing on an infeasible answer anyway. A
problem-specific mixer that moves *only* between feasible states removes
that burden structurally — exactly the kind of change that could unlock
more of QAOA's actual power. The catch: such mixers are notoriously hard
to build *efficiently*, and this repo is about doing that for one
concrete, practically important constraint class.

Can a QAOA mixer that never wastes search on infeasible configurations —
no penalty terms, no loops, just valid moves — scale to real electrical
grids, cheaply enough for both fault-tolerant and near-term (NISQ)
quantum hardware?

**Yes — with the right techniques.** Built directly, the mixer works
cleanly on synthetic feeders but fails outright on real, published
topology. Two fixes close that gap: **zone decomposition** (splits the
problem into small, exactly-solvable pieces) and a **cost-aware
bounded-witness mixer** (bounds circuit cost directly, at a small,
measured leakage price). Together they give both a fault-tolerant-ready
exact construction and a NISQ-ready cheap one — validated directly on
two real networks, not just synthetic proxies (see "Results" below).

## Setup: the constraint, the moves, and why some ties cost more than others

A distribution feeder has more switchable connections than it strictly
needs: a normally-closed backbone plus a small number of normally-open
**tie** switches, closed only temporarily (e.g. to reroute power around
a fault). The network must stay **radial** — every node reached, no
loops — at all times.

Qubit `i` = 1 iff switch `i` is closed; a configuration is feasible iff
its closed-switch set is a spanning tree of the network graph — a basis
of the network's graphic matroid.

![A feasible configuration (spanning tree) vs. an infeasible one (a loop)](results/illustration_feeder_problem.png)

A mixer moves probability between feasible states only, via **basis
exchanges** — close one switch, open another, land on a different tree.
Matroid theory guarantees exchanges like this reach any feasible
configuration from any other, so the whole mixer is built from one small
circuit term per exchange.

![One basis-exchange move: before, the move itself, and after](results/illustration_basis_exchange_move.png)

The catch: which switch has to open depends on the rest of the *current*
configuration, so each term generally must be **conditioned** on other
qubits — specifically, on the tie switch's **fundamental cycle** (the
other tree edges on the loop it would create). A short cycle means a
small, cheap condition; a long one means an expensive, deeply-controlled
gate. Real tie switches are deliberately placed to link *distant* parts
of a network for redundancy — the hard case, not the easy one.

![A tie edge's fundamental cycle, short vs long](results/illustration_fundamental_cycle.png)

When a term's conditioning set (its **witness**) is deliberately capped
smaller than the full cycle (technique 2, below), it becomes only
approximately correct — it can fire on states where it shouldn't.
**Leakage** is the resulting probability that ends up infeasible; **mean
feasible mass**, measured empirically per trajectory (`leakage_trace.py`
— exact for small instances, Wilson's-algorithm-sampled at scale), is the
complement — 1.0 means nothing leaked.

**The actual question this repo measures**: how much does conditioning
cost, in gates and depth, and does that cost grow or stay bounded as the
network gets larger and tie placement gets more realistic?

## Results: fault-tolerant now, NISQ-ready with decomposition

Techniques 1-2 (exact, or cost-aware bounded-witness, applied whole-graph,
no decomposition) already form a complete, valid construction at whatever
cost the search finds — fine for fault-tolerant hardware, which doesn't
care about cost the way NISQ does. Technique 3 (zone decomposition,
refined to guarantee its own cost) exists on top of that, to make the
same construction cheap enough to matter on NISQ hardware today, against
the feasibility arithmetic below (`fidelity ≈ (1-p)^N_CX`, published
two-qubit gate error rates):

| CX count | best-case trapped-ion (p=0.001) | typical superconducting (p=0.005) |
|---|---|---|
| 100 | 90% | 61% |
| 500 | 61% | 8% |
| 1,000 | 37% | 0.7% |
| 5,000 | 0.7% | ~0 |
| 13,000 | ~10⁻⁶ | ~0 |

### On synthetic data, the two tiers split cleanly

Technique 2 is stress-tested across tie placement (short vs. long-range —
real ties span 33-45% of network diameter, matching the long-range
generator) and tie-count growth (log-scaled, `k_ties(n) ≈ 1.43 ln n` —
real and benchmark networks show the ties-per-bus ratio dropping 6x from
15 to 179 buses, consistent with log growth — `docs/scaling-ladder-and-decomposition.md`
§9). Across that ladder it plateaus by `n_nodes=60` and stays fully
connected throughout — for short-range ties it peaks at `n_nodes=30`
first, a real, reproducible effect (7.9% seed-to-seed variance, not
noise, traced in §2) — but, per the table above, well past a comfortable
NISQ regime at these sizes either way.

Technique 3 fixes that without exception on this data — at every size
≥ 30 nodes, every seed, both conditions, zone decomposition (3a) alone
already costs less than the whole-graph construction (4.2x-41.5x
cheaper), and its cost-capped refinement (3b) guarantees the rest of the
way: every seed, every condition, every size tested, lands under 500 CX
(§6-8):

![Synthetic ladder: whole-graph -> zone decomposition -> cost-capped refinement](results/construction_progression_plot.png)

*(Technique 1, the exact construction, is deliberately absent here and
below — its failure mode is dropping candidates outright rather than
getting gradually more expensive, which makes it look artificially cheap
on cost and artificially perfect on safety, e.g. 75 CX, fully
disconnected, 43% of candidates dropped, at `long_log, n_nodes=30`.
Direct check: `docs/scaling-ladder-and-decomposition.md`.)*

![Synthetic ladder: the same three stages, measured for safety instead of cost](results/synthetic_mass_progression_plot.png)

Technique 2 alone leaks real, sometimes substantial probability (down to
91% mean feasible mass on the hardest long-range condition); 3a tightens
that considerably; 3b is indistinguishable from perfect (1.0 mean
feasible mass) at every size, on both conditions — cheaper AND safer than
either stage before it, not a tradeoff between the two.

At the hardest size tested (`n_nodes=150`), directly against the
feasibility numbers above:

| condition | construction | CX | reading |
|---|---|---|---|
| short-range, log growth | whole-graph | 2,745 | borderline on trapped-ion; not usable on superconducting |
| short-range, log growth | cost-capped | **73** | comfortably NISQ-ready |
| long-range, log growth | whole-graph | 10,712 | not usable on either device |
| long-range, log growth | cost-capped | **337** | comfortably NISQ-ready |

![Where the synthetic ladder lands relative to NISQ feasibility](results/synthetic_nisq_feasibility_plot.png)

### On real data, the same clean story holds

![Real networks: exact vs. zone decomposition vs. cost-capped refinement](results/real_network_comparison_plot.png)

| network | construction | CX | reading |
|---|---|---|---|
| CIGRE MV (15 buses, 3 ties) | exact whole-graph | 12,220 | fault-tolerant-ready; far outside NISQ range |
| CIGRE MV | zone decomposition | 6,771 (3,668-9,178) | cheaper on average, but unreliable — and still outside a comfortable NISQ regime even at its best |
| CIGRE MV | cost-capped refinement | **64** | comfortably NISQ-ready, perfectly safe, identical every seed |
| IEEE33 (33 buses, 5 ties) | exact whole-graph | 96 | 573/597 candidates dropped, disconnected — not usable at any cost |
| IEEE33 | zone decomposition | 132 | already comfortably NISQ-ready, identical every seed |
| IEEE33 | cost-capped refinement | **132** | matches zone decomposition exactly — nothing left on the table |

Same pattern as synthetic data, on both networks, no exceptions — but
CIGRE MV's zone-decomposition step is a genuinely unreliable number, not
just a point estimate: it ranges 3,668-9,178 CX across 5 tested seeds,
over 2x spread, from randomness in the witness search's restart order,
not the (fully deterministic) tree enumeration. Cost-capped's refinement
is immune to this on both networks — identical across every seed tested —
which is itself a reason to prefer it beyond just the mean cost: it's not
just cheaper, it's predictable. Both real networks land comfortably in
NISQ range this way.

![Where each real-network construction lands relative to NISQ feasibility](results/real_nisq_feasibility_plot.png)

**On scale beyond these two networks**: CIGRE MV and IEEE33 (15 and 33
buses) are the only published radial test feeders this repo runs the
*full circuit construction* on — such data is scarce at this scope. The
scaling **trend**, though, isn't resting on those two points alone: it's
established on the synthetic ladder up to `n_nodes=150` above, whose
parameters are calibrated against real topology statistics spanning 15
to 179 buses (`docs/scaling-ladder-and-decomposition.md` §9). That
179-bus network hasn't had the full construction run on it directly yet
— the natural next real-data point, not a gap in the claim itself.

**Full account**: `docs/scaling-ladder-and-decomposition.md`.
`docs/bounded-witness-mixer.md` covers the density-failure axis, the
bounded-witness construction itself, and the cost-aware search in full
detail.

## How it works: three techniques

**1. Exact matroid mixer** (`mixer.py`, `build_matroid_mixer`). For each
candidate exchange, brute-force search for the smallest witness that
makes the term's firing rule exactly correct — no leakage, verified
against the true validity function. Cheap and exact when witnesses stay
small, which needs both short-range ties and low tie density; when
either fails, the construction doesn't get gradually more expensive — it
drops candidates outright (51-95% in the failure cases tested) and the
resulting mixer stops being fully connected. Full derivation and
failure-mode trace: `docs/circuit-validity.md`, `docs/mixer-construction.md`,
`docs/bounded-witness-mixer.md` Finding 1.

**2. Bounded-witness mixer** (`truncated_mixer.py`) — bounds *cost*
directly, at a measured, small leakage cost. When no small exact witness
exists, don't drop the candidate: search for a witness of a fixed, capped
size instead, using a **majority-vote** validity rule, and accept the
resulting leakage — provided it's measured, not assumed, and stays small.

![Bounded-witness mixer concept: capping witness size trades cost for a measured, non-zero leakage rate](results/illustration_bounded_witness.png)

*(concept diagram, not a measurement)* capping to 1 of 3 witness qubits
and taking a majority vote gets one group right unanimously but leaks on
the other — a wider witness reduces leakage like this but costs more
circuit, which is why this technique cannot be used without its other
half: **making the search cost-aware (`cost_alpha`) is not optional**.
Picking witnesses by leakage alone gives an uncontrolled search every
incentive to walk to the cap every time — measured at 17,386–32,520 CX at
just 15-17 qubits, more than a real 37-qubit decomposed feeder circuit.
`cost_alpha` closes most of that gap and is this construction's validated
default (`docs/bounded-witness-mixer.md`).

**3a. Zone decomposition** (`zone_decomposition.py`) — fixes the *range*
failure structurally and exactly. Partition into zones via a min
tie-line-cut, solve each zone's matroid mixer independently, plus one
small assembly mixer over the contracted zone graph. Guaranteed **exact**
by graphic-matroid contraction/deletion — a structural fix to the
*constraint*, not an approximation. When the assembly graph is itself
still large or dense, the same idea applies one level up.

![A graph partitioned into zones, plus the contracted assembly problem](results/illustration_decomposition.png)

**3b. Cost-capped refinement** — picking a zone size up front only gets
you *a* cost, not a *controlled* one: decomposed cost has much higher
seed-to-seed variance than whole-graph cost (coefficient of variation up
to 99% at some sizes — a single "unlucky" zone can dominate a seed's
total). The fix is the same discipline as technique 2's: build each
subproblem, transpile it, and check its **actual** cost against a
threshold; sweep `cost_alpha` for the cheapest passing result if it's
already under, or try more than one zone-size granularity and recurse on
whichever gives the cheaper total if it's over.

![Cost-capped decomposition: measure actual cost, recurse only where it's over threshold](results/illustration_cost_capped_decomposition.png)

## What this does and doesn't demonstrate

This repo shows the mixer **construction** is correct and scales — both
exactly (fault-tolerant-relevant) and, via decomposition and cost-capping,
cheaply enough for a plausible NISQ target (real-topology-relevant). It
does **not** show a quantum algorithm outperforming a classical baseline:
there is no cost-Hamiltonian/oracle integration, no QAOA execution, no
classical solver comparison, and no test of the iterative boundary
coupling (an ADMM-style loop) this kind of decomposition would need for
the actual optimization objective — only the radiality constraint,
tested here, decomposes exactly. Classical formulations of the same
constraint still need explicit penalty terms or solver-enforced
constraints; this mixer builds feasibility into the dynamics directly —
that's the promise a bounded-cost result is a *precondition* for, not
proof of. See `methodology.md` for the precise boundary of what was
measured.

## How to reproduce

Three steps cover everything this README claims:

```
python -m venv .venv
source .venv/bin/activate   # or .venv\Scripts\activate on Windows
pip install -r requirements.txt

python scripts/verify_correctness.py && python scripts/verify_leakage_trace.py
# correctness: the exact construction is exactly leak-free and fully
# connected wherever it applies; the bounded-witness construction's
# leakage tooling (sparse tracer, Wilson's-algorithm sampling) is verified
# against an exact reference.

python scripts/run_real_networks_hierarchical.py
# the real-network numbers the Results section reports: exact vs.
# decomposed vs. cost-capped, both CIGRE MV and IEEE33.

python scripts/plot_results_figures.py && python scripts/plot_illustrations.py
# every figure in this README, regenerated from already-committed data.
```

All scripts are deterministic (fixed seeds); re-running reproduces the
committed `results/` files exactly, modulo `qiskit`/`networkx` version
differences. This reproduces what's *shown* here — the fuller
escalating-ladder, density-axis, and hierarchical-decomposition
investigations use the additional scripts indexed in
`docs/repository-map.md`.

## Scope

Measurement methodology and results for the question above only — not a
production mixer-compilation library, and not extended to constraint
classes other than graphic-matroid radiality (`methodology.md` has the
precise boundary).

## Repository layout

- `methodology.md` — graph generation, move generation, and exactly what
  "gate count" and "depth" mean.
- `docs/mixer-construction.md` — technical reference: matroid theory,
  exact circuit derivation, correctness-verification arguments.
- `docs/circuit-validity.md` — the real-topology failure, its root cause,
  the decomposition fix, and its scaling validation.
- `docs/bounded-witness-mixer.md` — the density-driven witness-blowup
  axis, the bounded-witness mixer, and the circuit-cost investigation
  behind `cost_alpha`.
- `docs/scaling-ladder-and-decomposition.md` — the escalating realism
  ladder, NISQ feasibility, hierarchical and cost-capped decomposition,
  and direct validation on two real networks. The fullest, most detailed
  account in this repo — everything summarized in this README's Results
  section traces back to a section here.
- `scripts/` and `results/` — one script per measurement, one CSV/plot
  pair per script, all generated, none hand-edited. Full script-by-script
  index: `docs/repository-map.md`.

## License

Apache License 2.0 — see `LICENSE`.
