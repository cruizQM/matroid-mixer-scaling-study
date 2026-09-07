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

That is a **partition matroid** (one independent choice per item, out of m options), not the general graphic matroid `mixer.py` handles. The practical consequence: the QAOA circuit needs only `k·m` qubits (one-hot per item), not one qubit per edge in the full gadget graph — and the mixer needs no witness search at all, because swapping one item's bucket choice is valid *unconditionally*, regardless of every other item's state. This is proven exactly, not sampled: `verify_partition_mixer.py`'s `verify_weight_preservation` checks, via the full unitary of the mixer's underlying two-qubit gate (the same RXX+RYY Givens-rotation `mixer.py`'s `_swap_block` already uses), that it preserves total Hamming weight for *every* possible input, and `verify_full_gadget_no_leakage` confirms the same end-to-end on a small full gadget instance.

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

## Honest scope of this document

- **This is not a demonstrated quantum advantage.** The circuit is shown to be small and to compile/transpile cheaply at the size where classical proof cost is already climbing steeply — a *necessary* condition for a quantum approach to be worth pursuing here, not a sufficient one. Solution quality (approximation ratio) is not measured anywhere in this repo: an early, separate, informal check at a trivial 12-qubit scale found the circuit's best sample matched the true optimum but its mean sample sat noticeably above it, with no improvement from more QAOA layers. That check is far too small to say anything about the actual hard instances (m=9-11 above), which are too large to simulate classically at all — the performance question stays genuinely open, not just unmeasured out of laziness.
- **This is a different mixer from `mixer.py`'s**, not a generalization of it. It works because this specific reduction happens to decompose into independent one-hot subproblems; that is a property of this instance family, discovered by reading the hardness proof's own structure, not a general recipe that would apply to an arbitrary hard matroid-basis problem.
- **E.ON's secondary requirement** (showing a classical tensor-network/MPS simulator also struggles) has not been attempted for this instance and remains open work.
- **The two negative attempts** (random quadratic cost, congestion cost) are reported here because they were real, informative dead ends — not polished away — and directly motivate why a structural construction was needed instead of a bigger real network.
