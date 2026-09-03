# The matroid basis-exchange mixer: full construction

This document is the technical reference for *what the mixer actually is*
— the matroid theory it relies on, the exact circuit it compiles to, and
why each design choice is correct — at a level of detail meant to let an
evaluator reconstruct or audit the construction without reading
`scripts/mixer.py` line by line. It complements, and does not repeat,
two other docs:

- `methodology.md` — what is measured (gate count, depth, scope of the
  sweep) and how graphs are generated.
- `docs/circuit-validity.md` — the narrative of the specific bug found
  during development (unconditional swap gates leak) and how it was
  caught. This document assumes that fix and explains the resulting
  construction on its own terms, with the math spelled out.

## 1. The constraint, as a matroid

A **matroid** `M = (E, I)` is a ground set `E` together with a family `I`
of "independent" subsets of `E` satisfying (among other axioms) the
**basis-exchange axiom**: if `B1, B2` are both maximal independent sets
("bases") and `e ∈ B1 \ B2`, there exists `f ∈ B2 \ B1` such that
`B1 - e + f` is also a basis.

The **graphic matroid** of a graph `G = (V, E)` takes independent sets to
be the forests of `G`; its bases are exactly the spanning trees. This is
the constraint this repo targets: **qubit `i` = 1 iff edge `i` of the
feeder graph is closed, and a computational basis state is feasible iff
its edge set is a spanning tree** — i.e., feasible states are exactly the
bases of the graphic matroid. `scripts/graphs.py::enumerate_spanning_trees`
computes this basis set `B ⊆ {0,1}^{n_qubits}` directly (brute force over
`C(n_qubits, k_ties)` edge subsets, cross-checked against Kirchhoff's
Matrix-Tree theorem — see `methodology.md`), rather than relying on any
property of `M` that isn't checked against the actual instance.

## 2. Why an exactly constraint-preserving mixer

Standard QAOA enforces constraints with a penalty term added to the cost
Hamiltonian, trading a tuned penalty weight for search over the *entire*
`2^n` space, most of which is infeasible. Hadfield et al.'s alternative
is a **mixer Hamiltonian `H_M`** that only ever moves probability
amplitude *within* the feasible subspace: for every feasible `x` and
infeasible `y`, `⟨y| U_M(β) |x⟩ = 0` for the mixer unitary
`U_M(β) = exp(-iβH_M)` at every `β`. This is the exact property checked
by `verify_no_leakage` / `verify_term_no_leakage` (§9). Restricting the
mixer's support to the feasible subspace this way means every QAOA
iterate stays feasible by construction, at the cost of a more structured
(and potentially more expensive) mixer circuit — which is precisely the
cost this repo's scaling study measures.

## 3. Move generation: basis exchanges as mixer terms

The basis-exchange axiom (§1) guarantees something stronger than "some
exchange sequence connects any two bases": it guarantees the **exchange
graph** — nodes = bases, an edge between `B1, B2` whenever they differ by
exactly one basis-exchange (`B2 = B1 - e + f`) — is connected using only
these single-element, "weight-2" moves (flip one 0-bit and one 1-bit).
No higher-weight moves are ever structurally necessary for a graphic
matroid. This is why every mixer term in this construction is a
transposition of exactly two qubits `(e, f)` — never a 3-or-more-qubit
permutation — which keeps each term's circuit small regardless of
instance size.

**Candidate generation** considers every pair of edges `(e, f)`,
`e < f`, `C(n_qubits, 2)` pairs total. For each pair, the actual
enumerated basis set `B` (not the matroid axiom in the abstract) is
consulted directly:

```
pairs(e, f) = { (index(T), index(T ^ swap_mask)) :
                T ∈ B, T ^ swap_mask ∈ B, index(T) < index(T ^ swap_mask) }
```

where `swap_mask = (1<<e) | (1<<f)`. A candidate with `pairs(e,f) = ∅` is
skipped (no real exchange exists for this pair in this instance — the
common case, since most edge pairs aren't part of any valid exchange).

**Selection is greedy, Kruskal-style**, via union-find over `B`: process
candidates in `(e, f)` lexicographic order; select a candidate iff it
connects two components of `B` not already connected by previously
selected candidates (`build_matroid_mixer` in `scripts/mixer.py`). Since
every candidate move has equal (weight-2) cost, any selection achieving
full connectivity of `B` is optimal in move-count terms — there is no
weighting problem to solve, unlike constraint classes without a basis-
exchange guarantee. The loop terminates early once `B` is a single
union-find component (`num_components == 1`); on a well-connected
instance this can select far fewer than `C(n_qubits,2)` candidates before
exiting.

**Complexity of this phase**: `O(n_qubits^2)` candidates, each requiring
an `O(|B|)` scan of the basis set — `O(n_qubits^2 |B|)` total before
accounting for witness search (§5, which dominates in practice).

## 4. Why an unconditional swap gate is wrong

The standard "safe-by-construction" building block for a transposition
constraint is the XY-hopping gate, built from `RXX(β) = exp(-iβ X⊗X/2)`
and `RYY(β) = exp(-iβ Y⊗Y/2)` (§6 derives its action exactly). This gate
provably does nothing to `|00⟩` and `|11⟩` — so for constraints where
"both flip or neither" is never itself the concern, it looks safe to
apply unconditionally to every selected `(e, f)`.

It is **not** safe here, because whether `T - e + f` is a valid basis
depends on `T`'s global structure — concretely, whether `f` crosses the
cut created by deleting `e` from the spanning tree `T` — not just on the
values of qubits `e, f`. An unconditional gate fires identically on
*every* basis state with exactly one of `e, f` set (a "trigger" state),
including states where `f` does not cross the right cut and the result is
not a spanning tree at all. §5 in `docs/circuit-validity.md` gives a
concrete failing example found this way. The fix is **witness
conditioning**, defined next.

## 5. Witness sets

For a selected exchange `(e, f)`, define the set of **trigger states**:

```
Trigger(e,f) = { T ∈ B : exactly one of bit e, bit f is set in T }
```

and for each `T ∈ Trigger(e,f)` its **validity**:

```
valid(T) = (T ^ swap_mask) ∈ B
```

— i.e., whether flipping `e` and `f` on `T` lands on another actual basis.
A set of qubits `W ⊆ {0,...,n_qubits-1} \ {e,f}` is a **witness set** for
`(e,f)` iff `valid(T)` is a function of `T` restricted to `W` alone: for
every `T1, T2 ∈ Trigger(e,f)`, `T1|_W = T2|_W ⟹ valid(T1) = valid(T2)`.
`find_witness_set` (`scripts/mixer.py`) searches for a **minimal** such
`W` by trying every subset of the other qubits in increasing size order,
`|W| = 0, 1, 2, ..., MAX_WITNESS_SEARCH_SIZE` (capped at 4), and returns
the first one found together with `valid_patterns` — the set of bit
patterns over `W` under which the exchange is valid. If no witness of
size `≤ 4` exists, the candidate is **dropped** (excluded from the mixer,
never included with an incorrect, unconditioned circuit) and
`dropped_candidates` is incremented.

**This stays small on the synthetic sweep, but does NOT on the whole real
graph** — see `docs/circuit-validity.md` for the full account, including
the fix. On every *synthetic* instance tested (`n_qubits` up to 35), no
witness set larger than 2 was needed, because that generator picks tie
edges as the geometrically *nearest* non-tree pairs, keeping each one's
fundamental cycle (§9) short (2-3 tree edges, measured directly). On the
real IEEE 33-bus topology (37 qubits), whose published tie switches
deliberately connect *distant* parts of the feeder, fundamental cycles
span 6-20 tree edges, and — because those cycles heavily overlap — even a
graph-derived (non-brute-force) witness candidate for the whole graph
averages 22 qubits. Witness size is fundamentally bounded by
fundamental-cycle length, not by sparsity (`k_ties`) alone; a sparse graph
can still have long cycles if its extra edges are long-range, which is
exactly the real-world case here. **The actual fix is not a smarter
witness search but decomposition**: partitioning the graph into zones
(`scripts/zone_decomposition.py`) shrinks each subproblem enough that
witness sizes return to single digits (measured: 0-4 across every tested
configuration, fully leak-free) — see `docs/circuit-validity.md` for the
numbers.

**Complexity, and why this is the actual bottleneck at real-feeder
scale**: for a fixed `(e,f)`, the search examines
`sum_{k=0}^{4} C(n_qubits - 2, k)` candidate witness subsets in the worst
case (a candidate that gets dropped exhausts the entire search), each
requiring a consistency check across `|Trigger(e,f)|` trees. On the
sparse synthetic instances this sweep uses (`n_qubits ≤ 35`, `|B|` in the
tens), this is fast. On the real IEEE 33-bus instance (`n_qubits = 37`,
`|B| = 50,751`), `|Trigger(e,f)|` can run into the thousands, and
`C(35, 4) ≈ 52,360` — a naive Python nested loop over both makes a single
dropped candidate take tens of seconds, dominating total construction
time. `find_witness_set` is implemented with the consistency check
vectorized as numpy bit-packed array operations rather than a per-tree
Python dict, for this reason; the search structure (increasing witness
size, minimal witness returned) is unchanged and was re-verified to
return identical results to the original implementation (same term
counts, same witness sizes, same `verify_no_leakage` outcome) on every
instance in `scripts/verify_correctness.py` before being trusted at real-
feeder scale.

## 6. The two-qubit exchange block, exactly

`_swap_block(β)` composes `RXX(β)` then `RYY(β)` on qubits `(e,f)`. Since
`(X⊗X)(Y⊗Y) = (XY)⊗(XY) = (iZ)⊗(iZ) = -(Z⊗Z) = (Y⊗Y)(X⊗X)`, these two
generators commute, so the composition is exactly
`exp(-iβ(X⊗X + Y⊗Y)/2)`. Evaluating `X⊗X + Y⊗Y` on the computational
basis:

- `(X⊗X + Y⊗Y)|00⟩ = |11⟩ + (i|1⟩)⊗(i|1⟩) = |11⟩ - |11⟩ = 0`, and
  symmetrically `(X⊗X + Y⊗Y)|11⟩ = 0` — the generator **annihilates**
  `|00⟩` and `|11⟩` exactly, so `exp(-iβ(X⊗X+Y⊗Y)/2)` acts as the
  identity there for every `β`. This is the "safe for free" property §4
  refers to.
- `(X⊗X + Y⊗Y)|01⟩ = |10⟩ + (i|1⟩)⊗(-i|0⟩) = |10⟩ + |10⟩ = 2|10⟩`, and
  symmetrically `(X⊗X + Y⊗Y)|10⟩ = 2|01⟩` — on the `{|01⟩, |10⟩}`
  subspace the generator acts as `2σ_x`. So

```
exp(-iβ(X⊗X+Y⊗Y)/2) |{01,10} =  [ cos β    -i sin β ]
                                  [ -i sin β  cos β   ]
```

At generic `β` this is a **partial rotation** — a superposition of
"stayed" and "swapped" — not a deterministic swap (that would require
`β = π/2`). This matters for §9's per-term verification, which checks
confinement to `{|01⟩, |10⟩}`, not collapse onto `|10⟩` alone.

## 7. Witness-conditioned circuit synthesis

`exchange_term_circuit(n_qubits, term, β)` builds one term's circuit:

- If `term.witness_qubits` is empty (witness size 0 — the exchange is
  unconditionally valid whenever triggered), the swap block from §6 is
  applied directly to `(e, f)`. This is the same "safe for free" case as
  the naive construction, used only where it is actually correct.
- Otherwise, for **each cube** in `term.minimized_cubes` — a
  Quine-McCluskey-minimized, exact cover of `term.valid_patterns`
  (`_minimize_patterns`; same truth table, not an approximation, so §9's
  verification needs no special-casing for it), where each cube is a
  length-`k` tuple of `{0, 1, None}` (`None` = don't-care) — the swap
  block is turned into a gate controlled on only the cube's non-`None`
  positions (a possibly-*strict subset* of `witness_qubits`, when some
  bits are don't-cares) and appended, controls first, then `(e, f)`.
  Because a computational basis state has one specific witness-qubit
  value, at most one cube's control condition is ever satisfied for a
  given basis input, so these sequential controlled applications don't
  interact with each other; the circuit is exactly the "apply the swap
  block iff the witness qubits match one of the valid patterns"
  operation, by construction rather than by argument after the fact
  (checked directly in §9). Compiling from `minimized_cubes` rather than
  `valid_patterns` directly changes only how many gates are used and how
  wide each one is, never which basis states the term acts on — see
  "does decomposition keep witness size bounded" in
  `docs/circuit-validity.md` for why this was added (compiling one gate
  per exact pattern, unminimized, was the actual driver of large
  transpiled CX counts on some instances, not witness size itself) and
  the measured improvement (6-8x on the worst cases, ~10-25% on typical
  ones, no case regressed).

`mixer_circuit` composes every selected term's circuit in edge-index
order at a single shared angle `β` (a `qiskit.circuit.Parameter` by
default, or a fixed float for measurement/verification — `β=0.37` is used
throughout this repo's measurements, chosen only to be a generic nonzero
value so no term's rotation accidentally vanishes; see `methodology.md`).

**Gate-count consequence**: a term with witness size `k` compiles to up
to `|minimized_cubes| ≤ |valid_patterns|` controlled two-qubit gates, each
with up to `k` controls but often fewer (don't-care bits). Generic
multi-controlled-unitary synthesis (what `qiskit.transpile` falls back to
for a non-standard controlled gate like this) costs `O(2^k)` CX gates in
the worst case for a `k`-controlled gate; empirically `k ∈ {0,1,2}` on the
synthetic sweep and up to 4 on decomposed real-feeder instances (§5),
which keeps individual gates small — but the *number* of valid patterns a
term needs (unrelated to `k` itself) was found to be the actual driver of
some instances' large transpiled CX counts, which minimization now
reduces directly rather than papering over with a bigger witness-size
budget.

## 8. Correctness verification

Hadfield's validity condition, restated precisely: for the full mixer
circuit's unitary `U`, and every basis (feasible) state `x ∈ B`,
`⟨y|U|x⟩ = 0` for every `y ∉ B`.

**Full-circuit check** (`verify_no_leakage`): builds `U` exactly via
`qiskit.quantum_info.Operator` and checks every column indexed by `x ∈ B`
directly. Exact (not statistical), but requires the full
`2^n_qubits × 2^n_qubits` unitary — infeasible in memory beyond roughly
`n_qubits ≈ 14`, and completely infeasible at the real-feeder scale this
repo also validates (`n_qubits = 37` would need `2^37 × 2^37` complex
amplitudes).

**Per-term local check** (`verify_term_no_leakage`): checks one term in
isolation, on just its own acting qubits (`witness_qubits ∪ {e,f}`,
typically 2-4 qubits total), by building only that term's small local
circuit and its `2^k × 2^k` unitary. This is tractable at any
`n_qubits`, since a term's circuit is provably identity outside its own
acting qubits — the full-register operator is `U_local ⊗ I_rest` up to
qubit relabeling.

**Why checking every term locally is equivalent to checking the full
circuit** (not just a heuristic substitute): a linear operator `U`
satisfies `U(F) ⊆ F` for a subspace `F` (here, the span of the feasible
basis states) iff `U|x⟩ ∈ F` for every basis vector `x` spanning `F` —
membership in a subspace is preserved under linear combination, so
checking basis vectors suffices. If every individual term `U_i` satisfies
`U_i(F) ⊆ F`, then by induction so does the composition
`U = U_T ⋯ U_2 U_1`. `verify_term_no_leakage` checks exactly this
per-term condition, but at the *local* level: it enumerates all `2^k`
local basis states of the term's acting qubits (a superset of the local
patterns that any real basis state in `B` actually projects to) and
confirms each is confined to `{state}` (identity, when not triggered or
the witness pattern isn't valid) or `{state, expected_target}` (the
generic-`β` rotation of §6, when it should fire). Because `valid_patterns`
was derived directly from real basis-set membership (§5) rather than
assumed, "should fire" here agrees exactly with "is a real, verified
exchange" for every basis state that actually occurs — checking the
(larger) set of all local combinations only makes the check more
conservative, not less correct. `scripts/verify_correctness.py`
cross-validates both methods agree on every small instance before the
per-term check is relied on alone at real-feeder scale
(`scripts/run_real_feeder_validation.py`) — this cross-check itself
caught a bug in an earlier version of `verify_term_no_leakage` (see its
docstring in `scripts/mixer.py`): it required amplitude to land
*entirely* on `expected_target`, which is wrong at generic `β` per §6's
derivation (a partial rotation, not a deterministic swap).

## 9. Complexity summary

| Phase | Cost | Dominant at |
|---|---|---|
| Spanning-tree enumeration | `O(C(n_qubits, k_ties) · n_nodes)` (union-find per candidate) | large `k_ties` |
| Candidate move scan | `O(n_qubits^2 · \|B\|)` | large `\|B\|` |
| Witness search (per candidate) | up to `O(sum_{k=0}^{4} C(n_qubits, k) · \|Trigger(e,f)\|)` | real-feeder scale (§5) — the observed bottleneck |
| Circuit transpilation | depends on witness sizes found (§7); `O(2^k)` CX per controlled term | large witness sizes, if any occurred |

On this repo's synthetic sweep (`n_qubits ≤ 35`, `k_ties = 3`, `\|B\|` in
the tens), every phase is fast (well under a second per instance). On the
real IEEE 33-bus instance (`n_qubits = 37`, `\|B\| = 50{,}751`), the
witness search was the measured bottleneck on the whole-graph attempt, at
tens of seconds per candidate that exhausted the search (roughly 4.4
hours total construction time before it was superseded by decomposition
— see `docs/circuit-validity.md`) — `results/real_feeder_results.csv`
has the final per-run numbers.

## 10. Relationship to other mixer constructions

This is a **Hadfield-style constraint-preserving mixer** specialized to
the graphic-matroid case, in the same family as XY-mixers for one-hot /
cardinality constraints — the qualitative difference is that a graphic
matroid's basis-exchange structure is *global* (validity depends on tree
structure, §4) rather than *local* (a one-hot constraint's validity is a
function of the two swapped qubits alone), which is exactly what makes
witness conditioning (§5) necessary here and unnecessary there. This
repo does not implement or compare against a general Fuchs-LX Pauli-term
emission pipeline (that is out of scope — see `README.md`'s "Scope"
section); it is a from-scratch, minimal construction specific to the
graphic-matroid constraint, built to measure this one question.
