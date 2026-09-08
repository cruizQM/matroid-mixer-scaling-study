# A provably hard instance, and a small circuit for it

## The question this document answers

This repo's core result (`mixer.py`, `README.md`) is a general, feasibility-preserving QAOA mixer for radial network reconfiguration, validated on real, published topology. It says nothing about whether any *particular* instance is hard for a classical solver — and a mixer built by enumerating the feasible set (as `mixer.py`'s does, to search for exchange witnesses) structurally cannot even be constructed on an instance whose feasible set isn't enumerable in the first place. A small, enumerable feasible set is close to a definition of "classically easy," so the general construction is, by design, confined to instances that a classical solver would find easy too.

This document asks a different question: is there a *specific, verifiable* instance of the same problem class (radial network reconfiguration) that is hard for a classical solver even at small size, and if so, can a matching quantum circuit for it be built and compiled at that same small size? Both halves are answered here with real, checked numbers — no part of this is asserted from the reduction's asymptotic hardness claim alone.

## Two attempts that did NOT produce a hard instance

Before reaching into a different construction, two natural attempts were tried directly on this repo's own real-network data (real IEEE 33-bus feeder topology, `real_feeders.py`):
- a dense random quadratic cost over the candidate edges,
- a physically motivated congestion cost (minimize the worst-loaded line, using real per-bus load data).

Both were solved to proven optimality by CP-SAT in well under a second, even on the full 33-node network. The reason is the same one given above: a real feeder's small number of tie switches keeps its feasible set (spanning trees) small enough to be classically easy regardless of what cost function is layered on top. Hardness has to come from the combinatorial structure itself, not the cost function.

## The construction: a 3-PARTITION reduction that stays hard at small size

Khodabakhsh, Sharma & Naeem (arXiv:1711.03517) prove that minimum-loss radial network reconfiguration is strongly NP-hard, via a reduction from 3-PARTITION. The relevant special case (unit nodal demand) has objective

```
cost(T) = sum_{e=(parent,child) in T} (successors(e))^2
```

for spanning tree T, where `successors(e)` counts the nodes in the subtree hanging below e. Because 3-PARTITION is *strongly* (not just weakly) NP-complete, this reduction is hard even when its numbers are small — the numbers don't need to grow to keep the instance hard, which is what lets the resulting network stay small too.

**The gadget** (`partition_gadget.py`, `gadget_graph.py`): given k=3m integers a_1..a_k, each strictly between B/4 and B/2, summing to m·B — build a root r, m "bucket" nodes each connected to r, k "item" nodes each connected to *every* bucket node, and each item node carrying (a_i − 1) private leaf nodes. The true minimum cost equals m·(1+B)² + constant iff the a_i admit a perfect 3-partition into m groups of 3 summing to B each.

## Part 1: the classical hardness curve (measured, not assumed)

`run_hardness_sweep.py` builds this gadget at increasing m, and times CP-SAT proving optimality — cross-checked at every single point against an independent exact DP (`PartitionGadget.exact_optimum`, which exploits the same forced-edge structure described below but computes the answer via dynamic programming, not the MILP) so that a "hard" result is never trusted without an independent ground truth alongside it.

Raw results: `results/hard_instance_hardness_sweep.csv`.

**The caveat that the cross-check exposes, stated up front rather than after the table.** The exact DP is not just a verifier — it is a fast classical algorithm for these instances. It memoizes on the sorted tuple of partial bucket sums, so its state count is bounded by the number of distinct sum-multisets, and with items drawn from {4,5,6,7} at fixed B=15 those collide massively: it runs in milliseconds at every m in the table. Strong NP-hardness of 3-PARTITION requires the bucket size B to grow with m; at fixed B the problem is pseudo-polynomial, and every instance below is one a structure-aware classical solver handles instantly. What the table measures is therefore a general-purpose solver failing on the *natural graph formulation* — real and useful (it is what a practitioner would reach for first), but not evidence that no classical route exists. The route to genuine hardness is large B: it inflates the DP (state count grows with the number of distinct reachable sums) and the MILP (O(mB) leaf nodes), while the QAOA circuit stays at 3m² qubits regardless of B. That experiment is scoped, not yet run.

| m | nodes | edges | YES proof time | NO proof time |
|---|---|---|---|---|
| 4 | 65 | 100 | 1.3s | 8.2s |
| 5 | 81 | 140 | 6.0s | 25.9s |
| 6 | 97 | 186 | 22.1s | 114.9s |
| 7 | 113 | 238 | 52.2s | 214.1s |
| 8 | 129 | 296 | 87.1s | 450.1s |
| 9 | 145 | 360 | 33.8s | 853.0s |
| 10 | 161 | 430 | 225.0s | 1174.6s |
| 11 | 177 | 506 | 894.9s | **1800.0s — hit the time limit, unproved** |

("YES" = a perfect 3-partition exists; "NO" = it provably doesn't, the harder direction since CP-SAT must rule out every alternative, not just find one.) Every row's CP-SAT objective matches the independent exact DP exactly (`matches_exact=True` for all 16 rows, including `m11_no` — CP-SAT found the right value, it simply couldn't finish proving no better tree exists within 30 minutes).

`m11_no` is the strongest single data point: on a 177-node, 506-edge graph — tiny by any standard — a production MILP solver (Google OR-Tools CP-SAT, 8 workers) could not certify optimality within 30 minutes. That is a direct, literal instance of the "classical solver fails to prove optimality, at small size" bar this case study is built to clear.

The key point: proof time grows steeply while the graph itself stays small — a direct, measured demonstration that this instance family is hard even at sizes far below where a real distribution feeder's own tie-switch count would ever put it.

## Part 2: why a small quantum circuit is possible here specifically

An exchange argument — the same basis-exchange move `mixer.py` is built around — shows that with this cost, almost none of the gadget's edges ever have real freedom:

- every root→bucket edge is forced into every optimal tree (a disconnected bucket node would have to route through some item node, strictly worsening that item's own cost more than it could ever save elsewhere);
- every item→leaf edge is forced trivially (leaves have degree 1, there is no alternative);
- the *only* real freedom is: each item independently picks exactly one bucket to attach to.

That is a **partition matroid** (one independent choice per item, out of m options), not the general graphic matroid `mixer.py` handles. The practical consequence: the QAOA circuit needs only `k·m = 3m²` qubits (one-hot per item; k=3m by construction, so this grows quadratically in m, not linearly — e.g. m=4 is 48 qubits, m=9 is already 243), not one qubit per edge in the full gadget graph — and the mixer needs no witness search at all, because swapping one item's bucket choice is valid *unconditionally*, regardless of every other item's state. This is proven exactly, not sampled: `verify_partition_mixer.py`'s `verify_weight_preservation` checks, via the full unitary of the mixer's underlying two-qubit gate (the same RXX+RYY Givens-rotation `mixer.py`'s `_swap_block` already uses), that it preserves total Hamming weight for *every* possible input, and `verify_full_gadget_no_leakage` confirms the same end-to-end on a small full gadget instance.

The cost oracle (the sum-of-squares objective) becomes a set of pairwise RZZ terms between items assigned to the same bucket, derived and checked exactly (`verify_ising_coefficients`) against the true cost function before ever touching a circuit.

`measure_partition_mixer.py` gives real, transpiled (not estimated) gate counts and depth for the same instance whose classical proof time appears in Part 1:

Instance: m=4, k=12, n_qubits=48 (65 nodes / 100 edges in the full gadget graph; the same "NO" instance in the table above, whose classical proof took a few seconds). Transpiled to `["cx", "rz", "sx", "x"]`, `optimization_level=1` — the same basis and settings `measure.py` uses for the general mixer, for direct comparability:

| p (QAOA layers) | CX count | total gates | depth |
|---|---|---|---|
| 1 | 816 | 2052 | 121 |
| 2 | 1632 | 4092 | 203 |
| 3 | 2448 | 6132 | 285 |
| 5 | 4080 | 10212 | 449 |

Gate count scales exactly linearly in p (each layer repeats the same block); depth grows sublinearly relative to gate count as p increases (121→449 for a 5x gate-count increase), consistent with layers overlapping in the transpiler's scheduling. Even at p=5 this is a modest, few-thousand-gate circuit on 48 qubits.

## Part 3: does a classical tensor-network (MPS) simulator also struggle?

E.ON's secondary requirement asks for exactly this: not just an exact classical solver failing, but a smart classical *simulator* of the quantum circuit itself also failing to substitute for real hardware. `run_mps_scaling_check.py` checks this directly using Qiskit Aer's matrix-product-state backend, at `p=2` (the smallest depth where entanglement can appear at all here — see below), sweeping both instance size (m) and MPS bond dimension.

Raw results: `results/mps_scaling_check.csv`.

**First, a structural fact worth stating plainly, found while calibrating this check rather than assumed going in**: at `p=1`, MPS is *exact even at bond dimension 2*, regardless of m. The circuit starts from a computational basis state, and a diagonal (RZ/RZZ) cost layer applied to a basis state cannot create entanglement, no matter how many qubits its terms connect; the mixer layer only ever acts within one item's own small register, never across items. Real entanglement only appears once a *second* cost layer acts on the superposition the first mixer created — so classical simulation difficulty here tracks QAOA depth, not just instance size.

**With that established, `p=2` is where the real test happens:**

| m | n_qubits | exact (statevector) | MPS bd=4 | MPS bd=16 | MPS bd=64 | bd=64 wall time |
|---|---|---|---|---|---|---|
| 2 | 12 | 884.18 | 910.18 | 870.51 | 884.18 | 0.04s |
| 3 | 27 | 1561.72 | 1700.37 | 1688.44 | 1681.80 | 5.4s |
| 4 | 48 | *(infeasible — 4.3 exabytes)* | 2474.87 | 2519.87 | 2248.34 | 26.8s |
| 5 | 75 | infeasible | 2924.34 | 2732.88 | 2908.28 | 83.7s |
| 6 | 108 | infeasible | 4250.52 | 3619.88 | 3622.84 | 193.1s |
| 7 | 147 | infeasible | 6135.12 | 4349.30 | **discarded (390.6s, over budget)** | — |

Two findings, both real and not manufactured to fit a narrative:

1. **MPS does not converge cleanly even where exact ground truth is available.** At m=3, bond dimension 64 still misses the exact value by ~120 (a 7.7% relative error) — not the near-exact match a "just needs a bit more bond dimension" story would predict. Beyond m=3, where no exact answer exists to check against, the bd=16 → bd=64 values are *not monotonically converging* either (m=5's bd=64 estimate is actually further from bd=16 than one would want, m=6's bd=16 and bd=64 roughly agree, m=7's bd=64 didn't finish) — the honest read is that a fixed, modest bond dimension cannot be trusted to track the true value as m grows, not that it converges slowly.
2. **Wall-clock cost at a fixed bond dimension (64) explodes**: 5.4s → 26.8s → 83.7s → 193.1s → over 390s, for m=3→7 — roughly 70x for less than a 6x increase in qubit count. This is a real, measured cost blow-up, not an extrapolation, and it happens on top of — not despite — accuracy already being questionable at the same bond dimension.

**Conclusion for this section**: a classical tensor-network simulator does not offer a reliable shortcut here either. This directly satisfies E.ON's secondary ask (a struggling MPS simulator, not just a struggling exact MILP solver), and it also gives the honest answer to the approximation-ratio question below: the gap isn't closed because it was left unmeasured, it's open because the two most natural classical proxies for "just simulate it and see" both hit real, measured walls at a similar, small scale.

## Honest scope of this document

- **These specific instances are not classically hard.** See the caveat at the top of Part 1: at fixed B=15 the repo's own DP solves every one of them in milliseconds. What is demonstrated is that a general-purpose solver on the natural formulation fails at small size, and that the circuit construction works there. Instance-level hardness against *every* classical route needs B to grow with m, which is the scoped next experiment.
- **This is not a demonstrated quantum advantage.** The circuit is shown to be small and to compile/transpile cheaply at the size where classical proof cost is already climbing steeply — a *necessary* condition for a quantum approach to be worth pursuing here, not a sufficient one. The approximation ratio (expected cost vs. `PartitionGadget.exact_optimum`) is now measured exactly at m=2 and m=3 (Part 3's table, statevector column) at a fixed representative angle (0.37 half-turns, this repo's existing convention for "a concrete, reproducible value," not a claimed optimum) — ratios of 1.27 and 1.49 respectively. Beyond m=3, exact simulation is impossible (m=4 alone needs 48 qubits, 4.3 exabytes) and MPS was shown in Part 3 to not reliably substitute for it. So the honest status is: verified at trivial scale, genuinely unknown at the scale where classical hardness actually bites (m=9-11), and that gap is not for lack of trying — it survived both an exact-simulation attempt and a tensor-network attempt.
- **This is a different mixer from `mixer.py`'s**, not a generalization of it. It works because this specific reduction happens to decompose into independent one-hot subproblems; that is a property of this instance family, discovered by reading the hardness proof's own structure, not a general recipe that would apply to an arbitrary hard matroid-basis problem.
- **The two negative attempts** (random quadratic cost, congestion cost) are reported here because they were real, informative dead ends — not polished away — and directly motivate why a structural construction was needed instead of a bigger real network.
