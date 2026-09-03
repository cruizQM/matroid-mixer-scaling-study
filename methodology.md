# Methodology

## The question this measures

For a QAOA mixer that preserves the constraint "the selected edges form a
spanning tree" (a basis of the graph's graphic matroid), does the compiled
mixer circuit's gate count and depth grow, stay flat, or shrink as the
graph gets larger — for sparse graphs resembling real distribution
feeders? See `docs/circuit-validity.md` for what "the mixer" means
precisely and how its correctness was verified; this document covers graph
generation, move generation, and the measurement definitions.

## Graph generation (`scripts/graphs.py`)

Each instance is a synthetic **feeder graph**: `n_nodes` points placed
uniformly at random in the unit square, connected by a Euclidean minimum
spanning tree (a defensible model of a feeder's primary backbone — minimum
conductor length is a real design objective), plus `k_ties` extra edges
connecting the geometrically nearest node pairs not already in the tree
(modelling normally-open tie switches used for reconfiguration/
redundancy). Edge count is therefore exactly `n_nodes - 1 + k_ties`.

**This is a synthetic model, not derived from real feeder topology data.**
It is chosen for reproducibility and structural resemblance to real
feeders (sparse, spatial, mostly-tree-plus-a-few-ties) — not because it has
been validated against utility network data.

**`k_ties` is held fixed across the whole sweep** (`k_ties=3`). This
matters for two reasons: (1) it keeps the graph's cycle space dimension
(number of independent cycles) constant regardless of `n_nodes`, which is
the property `docs/circuit-validity.md` argues keeps the witness-set size
bounded; (2) a fixed, small tie-switch count independent of feeder size is
also a reasonable model of real feeders — adding busbars doesn't
necessarily add proportionally more reconfiguration switches.

**Minimum `n_nodes` is chosen, not arbitrary.** With `k_ties` fixed, a
graph is only genuinely *sparse* once `n_nodes` is large enough that
`n_qubits = n_nodes - 1 + k_ties` is small relative to the `C(n_nodes, 2)`
possible edges. At `n_nodes=4, k_ties=3`, the graph has 6 edges on 4
nodes — that *is* the complete graph K4, not sparse at all. This was
caught by inspecting an anomalously expensive small instance during
development, not assumed away. `scripts/run_scaling_study.py` computes the
smallest `n_nodes` at which edge density (`n_qubits / C(n_nodes,2)`) drops
below 0.4, and starts the sweep there (`n_nodes=8` for `k_ties=3`).

## Constraint and move generation (`scripts/mixer.py`)

Qubit `i` = 1 iff edge `i` is closed. Feasible states are exactly the
spanning trees of the feeder graph (bases of its graphic matroid) — see
`docs/mixer-construction.md` for the full matroid formalization and
exact circuit-level derivation of the mixer built on it.
Feasible-set enumeration (`scripts/graphs.py::enumerate_spanning_trees`)
is brute force over all size-`(n_nodes-1)` edge subsets — `C(n_qubits,
k_ties)`, polynomial in `n_nodes` for the fixed, small `k_ties` this study
uses (verified against Kirchhoff's Matrix-Tree theorem during
development, not just assumed correct).

Candidate mixer moves are every pair of distinct edges `(e,f)`
(`C(n_qubits,2)`, cheap since `n_qubits` stays small for sparse graphs).
Each candidate is checked against the actual enumerated spanning-tree set
— an edge pair that doesn't connect any previously-disconnected pair of
trees is discarded. Selected candidates are added greedily via a
union-find over the spanning-tree set until every spanning tree is
connected to every other by some sequence of selected moves (all
candidates have equal weight-2 cost here, since the matroid basis-exchange
axiom guarantees single-edge exchanges always suffice — no higher-weight
moves are needed, unlike more general constraint classes).

## Circuit construction and validity ("gate count" / "depth" precisely)

Each selected exchange becomes a circuit term: an `RXX`+`RYY` swap block on
qubits `e,f`, made conditional on a small "witness" qubit set when
(as is usually the case) the exchange's validity depends on other qubits —
see `docs/circuit-validity.md` for why this conditioning is necessary and
how it was derived and verified. **This "small witness" property holds
for the synthetic graph family this sweep uses, but was found not to hold
on the whole real IEEE 33-bus graph** — see `docs/circuit-validity.md`
and `scripts/investigate_fundamental_cycles.py` for why, and for how this
was resolved (zone decomposition, `scripts/zone_decomposition.py`,
restoring witness size ≤4 on every tested configuration). The results and
claims in this document are scoped to the synthetic sweep only; see
`docs/circuit-validity.md` and `results/zone_decomposition_results.csv`
for the real-feeder numbers. The full mixer circuit is the composition
of all selected terms, in edge-index order, at a single fixed Trotter
angle (`beta=0.37`, an arbitrary nonzero value — chosen only so that no
term's rotation accidentally vanishes; not tuned).

**"Gate count" and "depth" are measured on the transpiled circuit**, not
estimated analytically:

- Backend: `qiskit.transpile`, `optimization_level=1`.
- Basis gates: `["cx", "rz", "sx", "x"]` — representative of gate-model
  hardware (single-qubit rotations + a two-qubit entangler), not a claim
  about any specific device or vendor.
- **"CX count"** = number of `cx` operations in the transpiled circuit
  (`qiskit`'s `count_ops()`); **"depth"** = `QuantumCircuit.depth()` on the
  same transpiled circuit.
- **"Qubit count" per instance** = `n_qubits = n_nodes - 1 + k_ties`
  (one qubit per feeder edge) — reported directly in every results row,
  never inferred.

## What is and isn't measured

- **Measured, per instance, for every size in the sweep, across 5 random
  seeds**: qubit count, number of spanning trees, mixer term count, the
  largest witness-set size used, whether any candidate exchange had to be
  dropped (no small witness found), whether the resulting mixer is fully
  connected, transpiled CX count, transpiled total gate count, and
  transpiled depth. Every field is in `results/scaling_results.csv`
  (per-seed) and `results/scaling_summary.csv` (mean/min/max per size).
- **Not measured**: hardware execution, noise, real device compilation
  (routing/layout for a specific coupling map), or QAOA solution quality
  (approximation ratio) — this study is scoped to mixer *construction*
  cost only, per the open question it answers.
- **Range actually tested**: `n_qubits` from 10 to 35 (`n_nodes` 8 to 33),
  5 seeds per size. This range is what brute-force spanning-tree
  enumeration keeps fast (well under a second per instance in this range);
  it is not extrapolated beyond what is in `results/`.

## Reproducing

See `README.md`'s "How to reproduce" for the full, current command list
(this document doesn't duplicate it, to avoid the two drifting out of
sync as scripts are added). All scripts are deterministic (fixed seeds);
re-running should reproduce the committed CSVs exactly, modulo `qiskit`/
`networkx` version differences in transpilation.
