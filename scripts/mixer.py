"""Matroid basis-exchange mixer construction and circuit synthesis.

Constraint: the selected edge subset must form a spanning tree of the
feeder graph (a basis of its graphic matroid). This is the natural QAOA
encoding of radial-operation constraints in distribution-feeder
reconfiguration: qubit i = 1 iff edge i is closed, feasible states are
exactly the spanning trees.

## Move generation

The matroid basis-exchange axiom guarantees that for any two bases
(spanning trees) B1, B2 and any e in B1\\B2, there is an f in B2\\B1 with
B1-e+f also a basis -- i.e. the basis-exchange graph (nodes = spanning
trees, edges = single-element exchanges) is *always* connected, using only
weight-2 moves. Candidate exchange pairs are every pair of distinct edges
(e, f); each candidate is verified against the actual enumerated
spanning-tree set (never assumed), and only pairs that connect previously
disconnected components are selected (Kruskal-style, all weight-2
candidates being equal cost).

## Circuit synthesis: why this needs conditioning, unlike the pairwise case

An early version of this assumed the standard XY-hopping gate
(X_eX_f+Y_eY_f)/2 -- which annihilates |00>,|11> on e,f, so it seemed safe
"for free" -- would be enough, by analogy with simpler one-hot/cardinality
constraints. Direct verification (`verify_no_leakage`) caught that this is
WRONG here: whether exchanging e for f preserves the spanning-tree
property depends on whether f crosses the cut created by removing e from
the *current* tree -- a property of the tree's global structure, not just
e and f's own values. An unconditional swap gate fires on every state with
exactly one of e,f active, including trees where the exchange is invalid,
leaking probability outside the feasible subspace.

**The fix, found by direct search rather than assumed**: for every
selected exchange, `find_witness_set` brute-force-searches (bounded, small
cap) for a *minimal set of other qubits* whose values alone determine
validity, then the swap is applied only when the witness qubits match a
pattern under which the exchange is actually valid, verified against the
full enumerated spanning-tree set. On the synthetic sparse-feeder family
this repo's scaling sweep uses (see docs/circuit-validity.md), the
minimal witness never exceeded 2 qubits -- plausibly because a sparse
graph whose extra edges are geometrically short-range has few alternative
ways to reconnect around any given edge, so "does this edge cross the
right cut" reduces to checking a small, fixed number of alternative
edges. **This bound does NOT hold in general** -- on the real, published
IEEE 33-bus feeder, whose tie switches are long-range by design, the
whole-graph witness requirement reaches into the dozens of qubits (see
docs/circuit-validity.md's real-feeder section); `real_feeders.py`,
`zone_decomposition.py` and friends exist specifically because bounded
witness search alone was found insufficient there, and needed
decomposition to fix. If no witness is found within the search cap, that
candidate is dropped (not included with an incorrect circuit) and the
greedy selection tries other candidates.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np
from qiskit import QuantumCircuit
from qiskit.circuit import Parameter
from qiskit.quantum_info import Operator

from graphs import FeederGraph

MAX_WITNESS_SEARCH_SIZE = 4


class UnionFind:
    def __init__(self, n: int):
        self.parent = list(range(n))

    def find(self, x: int) -> int:
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, x: int, y: int) -> bool:
        rx, ry = self.find(x), self.find(y)
        if rx == ry:
            return False
        self.parent[rx] = ry
        return True


@dataclass(frozen=True)
class ExchangeTerm:
    e: int
    f: int
    witness_qubits: Tuple[int, ...]
    valid_patterns: Tuple[Tuple[int, ...], ...]  # each entry: witness_qubits values under which the swap is valid
    minimized_cubes: Tuple[Tuple[Optional[int], ...], ...] = ()
    # each entry: witness_qubits values (0/1) or None (don't-care / uncontrolled) under which the
    # swap is valid -- a Quine-McCluskey-minimized EXACT cover of `valid_patterns` (same truth
    # table, fewer and/or narrower controlled gates), used by exchange_term_circuit instead of
    # valid_patterns directly. See `_minimize_patterns`.

    @property
    def control_count(self) -> int:
        return len(self.witness_qubits)

    @property
    def max_cube_controls(self) -> int:
        """Widest single controlled gate exchange_term_circuit actually
        emits, after minimization -- can be smaller than control_count
        when some valid patterns share a don't-care bit."""
        return max((sum(1 for v in cube if v is not None) for cube in self.minimized_cubes), default=0)

    @property
    def minimized_gate_count(self) -> int:
        return len(self.minimized_cubes)


@dataclass
class MixerConstruction:
    n_qubits: int
    spanning_trees: List[int]
    terms: List[ExchangeTerm]
    fully_connected: bool
    dropped_candidates: int  # candidates that connected something but had no small witness, so were skipped


def fundamental_cycle_structure(graph: FeederGraph) -> Dict[int, Set[int]]:
    """For an arbitrary reference spanning tree `T0` of `graph`, returns
    `cycle_of[i]` = the set of edge indices making up the fundamental
    cycle(s) edge `i` participates in: if `i` is a `T0` non-tree edge, the
    tree-path between its endpoints plus `i` itself; if `i` is a `T0`
    tree edge, the union of every non-tree edge's cycle that happens to
    pass through `i`. `T0`'s specific choice doesn't matter for
    correctness -- it only produces a candidate witness set, which is
    always verified against the real enumerated spanning-tree set before
    being trusted (`find_witness_set` below), never assumed correct from
    graph theory alone.

    An edge with an EMPTY `cycle_of[i]` is a bridge of the whole graph
    `graph` (present in every one of its spanning trees), by a standard
    cycle-space fact: the `r` fundamental cycles of any reference tree
    form a basis (over GF(2)) for the graph's entire cycle space, so a
    tree edge that is a member of none of them individually cannot appear
    in any GF(2) combination of them either -- i.e. it lies on no cycle
    of the graph at all, which is exactly the definition of a bridge. See
    `docs/mixer-construction.md` and `docs/circuit-validity.md`'s
    real-feeder section for why this structural approach was needed: a
    graph's cycle RANK being small (`k_ties`, bounded) does not bound
    individual cycles' LENGTH, and a long fundamental cycle is exactly
    what defeats a small brute-force witness search.
    """
    g = graph.networkx_graph()
    t0 = nx.minimum_spanning_tree(g)
    t0_edges = {tuple(sorted(e)) for e in t0.edges()}
    edge_index = {e: i for i, e in enumerate(graph.edges)}
    non_tree = [i for i, e in enumerate(graph.edges) if e not in t0_edges]

    cycle_of: Dict[int, Set[int]] = {i: set() for i in range(graph.n_edges)}
    for gi in non_tree:
        u, v = graph.edges[gi]
        path_nodes = nx.shortest_path(t0, u, v)
        path_edges = {
            edge_index[tuple(sorted((path_nodes[k], path_nodes[k + 1])))]
            for k in range(len(path_nodes) - 1)
        }
        cycle = path_edges | {gi}
        cycle_of[gi] |= cycle
        for te in path_edges:
            cycle_of[te] |= cycle
    return cycle_of


def _minimize_patterns(
    k: int, valid_patterns: Tuple[Tuple[int, ...], ...]
) -> Tuple[Tuple[Optional[int], ...], ...]:
    """Quine-McCluskey minimization of the Boolean function whose ON-set
    is EXACTLY `valid_patterns` (over `k` witness-qubit variables) and
    whose OFF-set is every other length-`k` bit pattern -- no external
    don't-cares, so the result represents exactly the same truth table,
    just as a cover of (possibly fewer, possibly narrower-controlled)
    cubes instead of one gate per exact pattern.

    Exists because `exchange_term_circuit` originally compiled one
    controlled gate per entry of `valid_patterns` (up to `2^k` of them),
    which was found (`docs/circuit-validity.md`) to be the actual driver
    of large transpiled CX counts on some instances -- not witness size
    itself, which stays small, but the number of valid patterns a term
    happens to need, unrelated to network size. If two valid patterns
    differ in only one bit (e.g. `1011` and `1010`), a controlled gate
    conditioned on `101-` (3 controls, "-" = don't-care/uncontrolled)
    fires on both and on nothing else -- fewer gates, and narrower ones,
    for exactly the same validity function.

    Returns a minimal set of cubes, each a length-`k` tuple with values
    in `{0, 1, None}` (`None` = don't-care). Because this is an exact
    (not approximate) cover of the same ON-set, `verify_term_no_leakage`
    and `verify_all_terms_no_leakage` need no changes -- their
    `should_fire` logic is defined by membership in `valid_patterns`,
    which this doesn't alter.
    """
    if k == 0:
        return ((),) if valid_patterns else ()
    if not valid_patterns:
        return ()

    minterms = [tuple(p) for p in valid_patterns]
    current: Set[Tuple[Optional[int], ...]] = set(minterms)
    prime_implicants: Set[Tuple[Optional[int], ...]] = set()

    while current:
        used: Set[Tuple[Optional[int], ...]] = set()
        next_level: Set[Tuple[Optional[int], ...]] = set()
        current_list = list(current)
        for i in range(len(current_list)):
            for j in range(i + 1, len(current_list)):
                a, b = current_list[i], current_list[j]
                diff = [idx for idx in range(k) if a[idx] != b[idx]]
                if len(diff) == 1:
                    idx = diff[0]
                    if a[idx] is None or b[idx] is None:
                        continue
                    merged = list(a)
                    merged[idx] = None
                    merged_t = tuple(merged)
                    next_level.add(merged_t)
                    used.add(a)
                    used.add(b)
        prime_implicants |= current - used
        current = next_level

    def covers(cube: Tuple[Optional[int], ...], minterm: Tuple[int, ...]) -> bool:
        return all(c is None or c == m for c, m in zip(cube, minterm))

    pi_list = list(prime_implicants)
    coverage = {pi: {m for m in minterms if covers(pi, m)} for pi in pi_list}

    selected: List[Tuple[Optional[int], ...]] = []
    remaining = set(minterms)
    while remaining:
        essential = None
        for m in remaining:
            covering = [pi for pi in pi_list if m in coverage[pi]]
            if len(covering) == 1:
                essential = covering[0]
                break
        if essential is not None:
            chosen = essential
        else:
            chosen = max(
                pi_list,
                key=lambda pi: (len(coverage[pi] & remaining), -sum(1 for x in pi if x is not None)),
            )
        if chosen not in selected:
            selected.append(chosen)
        remaining -= coverage[chosen]

    return tuple(selected)


def _verify_candidate_witness(
    candidate: Tuple[int, ...], trigger: List[int], validity: Dict[int, bool]
) -> Optional[Tuple[Tuple[int, ...], ...]]:
    """Checks whether `candidate` (a specific, already-chosen qubit set,
    as opposed to `find_witness_set`'s search over many candidates) is a
    valid witness set: whether `validity` is in fact a function of `t`
    restricted to `candidate`, for every trigger tree `t`. Returns the
    valid-pattern tuple if so, else None -- never assumed, always checked
    against the real trees passed in."""
    pattern_validity: Dict[Tuple[int, ...], bool] = {}
    for t in trigger:
        key = tuple((t >> q) & 1 for q in candidate)
        v = validity[t]
        if key in pattern_validity and pattern_validity[key] != v:
            return None
        pattern_validity[key] = v
    return tuple(k for k, v in pattern_validity.items() if v)


def find_witness_set(
    trees: List[int],
    tree_set: Set[int],
    n_qubits: int,
    e: int,
    f: int,
    max_size: int = MAX_WITNESS_SEARCH_SIZE,
    cycle_of: Optional[Dict[int, Set[int]]] = None,
    prefer_structural: bool = False,
) -> Tuple[Optional[Tuple[int, ...]], Optional[Tuple[Tuple[int, ...], ...]]]:
    """Two independent ways to find a witness set, tried in an order
    controlled by `prefer_structural`:

    1. **Bounded brute-force search** (the original approach -- see
    `docs/circuit-validity.md`): every subset of the other qubits, in
    increasing size order up to `max_size`, checked for consistency.
    Guaranteed to find the minimal witness when one exists within the
    cap, but its cost grows combinatorially in the cap and in
    `len(trigger)` -- vectorized with numpy (bit-packed array ops instead
    of a Python dict keyed by tuples) since this is the actual bottleneck
    once `trigger` reaches into the thousands of trees (the real IEEE
    33-bus instance) and a candidate without a small witness forces the
    search out to size 3-4 (`C(35,4)~52k` subsets) before giving up.

    2. **Graph-derived structural candidates** (`fundamental_cycle_structure`,
    only used if `cycle_of` is supplied): the union of `e` and `f`'s own
    fundamental cycles (cheap, usually small), then -- guaranteed
    sufficient by the cycle-space argument in `fundamental_cycle_structure`'s
    docstring -- the full union of every fundamental cycle in the graph.
    Each candidate is directly verified (`_verify_candidate_witness`)
    against the real trigger trees, never assumed correct from the graph
    theory alone. This exists because brute force provably cannot succeed
    when the true minimal witness is larger than a graph's fundamental
    cycles justify -- see `docs/circuit-validity.md`'s real-feeder
    section: real (as opposed to this repo's synthetic) tie switches
    connect distant parts of a feeder, giving fundamental cycles of 6-20
    tree edges, far past any brute-force-reachable cap.

    `prefer_structural=False` (the default) tries brute force first --
    this reproduces every previously-verified synthetic-sweep result
    exactly, since brute force already succeeds (and returns the true
    minimal witness) on every synthetic instance tested, so structural
    candidates are never reached there. `prefer_structural=True` tries
    the cheap structural candidates first instead, skipping the expensive
    brute-force search on candidates it would only fail on anyway -- used
    for real-feeder validation, where brute force alone left 91% of
    candidates dropped (see `run_real_feeder_validation.py`).

    Returns (witness_qubits, valid_patterns), or (None, None) if nothing
    was found.
    """
    swap_mask = (1 << e) | (1 << f)
    trigger = [t for t in trees if bin(t & swap_mask).count("1") == 1]
    if not trigger:
        return None, None
    validity: Dict[int, bool] = {t: ((t ^ swap_mask) in tree_set) for t in trigger}

    def try_brute_force() -> Tuple[Optional[Tuple[int, ...]], Optional[Tuple[Tuple[int, ...], ...]]]:
        trigger_arr = np.array(trigger, dtype=np.int64)
        validity_arr = np.array([validity[t] for t in trigger], dtype=np.int64)
        other_qubits = [q for q in range(n_qubits) if q not in (e, f)]
        for size in range(0, max_size + 1):
            for subset in itertools.combinations(other_qubits, size):
                key = np.zeros(len(trigger_arr), dtype=np.int64)
                for i, q in enumerate(subset):
                    key |= ((trigger_arr >> np.int64(q)) & 1) << i
                combined = key * 2 + validity_arr
                uniq_keys, first_idx = np.unique(key, return_index=True)
                if len(np.unique(combined)) == len(uniq_keys):
                    group_validity = validity_arr[first_idx]
                    valid_patterns = tuple(
                        tuple(int((int(k) >> i) & 1) for i in range(size))
                        for k in uniq_keys[group_validity == 1]
                    )
                    return subset, valid_patterns
        return None, None

    def try_structural() -> Tuple[Optional[Tuple[int, ...]], Optional[Tuple[Tuple[int, ...], ...]]]:
        if cycle_of is None:
            return None, None
        local = tuple(sorted((cycle_of.get(e, set()) | cycle_of.get(f, set())) - {e, f}))
        result = _verify_candidate_witness(local, trigger, validity)
        if result is not None:
            return local, result
        full_union = tuple(sorted(set().union(*cycle_of.values()) - {e, f}))
        result = _verify_candidate_witness(full_union, trigger, validity)
        if result is not None:
            return full_union, result
        return None, None

    first, second = (try_structural, try_brute_force) if prefer_structural else (try_brute_force, try_structural)
    witness, patterns = first()
    if witness is not None:
        return witness, patterns
    return second()


def build_matroid_mixer(
    graph: FeederGraph,
    spanning_trees: List[int],
    verbose: bool = False,
    prefer_structural: bool = False,
) -> MixerConstruction:
    """`prefer_structural=True` computes the graph's fundamental-cycle
    structure once up front and tries it before brute-force witness
    search on every candidate -- see `find_witness_set`'s docstring for
    why this matters and when to use it (real-feeder-scale instances,
    where brute force alone leaves most candidates dropped)."""
    import time as _time

    n_qubits = graph.n_edges
    index_of: Dict[int, int] = {mask: i for i, mask in enumerate(spanning_trees)}
    tree_set = set(spanning_trees)
    N = len(spanning_trees)
    cycle_of = fundamental_cycle_structure(graph) if prefer_structural else None

    uf = UnionFind(N)
    num_components = N
    terms: List[ExchangeTerm] = []
    dropped = 0
    total_candidates = n_qubits * (n_qubits - 1) // 2
    examined = 0

    for e in range(n_qubits):
        if num_components == 1:
            break
        for f in range(e + 1, n_qubits):
            if num_components == 1:
                break
            examined += 1
            t0 = _time.perf_counter() if verbose else 0.0
            swap_mask = (1 << e) | (1 << f)
            pairs = []
            for mask in spanning_trees:
                other = mask ^ swap_mask
                if other in index_of and other != mask:
                    a, b = index_of[mask], index_of[other]
                    if a < b:
                        pairs.append((a, b))
            if not pairs:
                if verbose:
                    print(f"    [{examined}/{total_candidates}] ({e},{f}): no valid exchange pairs, skip")
                continue
            connects_something = any(uf.find(a) != uf.find(b) for a, b in pairs)
            if not connects_something:
                if verbose:
                    print(f"    [{examined}/{total_candidates}] ({e},{f}): {len(pairs)} pairs, none new, skip")
                continue

            witness_qubits, valid_patterns = find_witness_set(
                spanning_trees, tree_set, n_qubits, e, f, cycle_of=cycle_of, prefer_structural=prefer_structural
            )
            elapsed = _time.perf_counter() - t0 if verbose else 0.0
            if witness_qubits is None:
                dropped += 1
                if verbose:
                    print(f"    [{examined}/{total_candidates}] ({e},{f}): no witness found within cap, DROPPED ({elapsed:.2f}s)")
                continue

            minimized_cubes = _minimize_patterns(len(witness_qubits), valid_patterns)
            terms.append(
                ExchangeTerm(
                    e=e, f=f, witness_qubits=witness_qubits, valid_patterns=valid_patterns,
                    minimized_cubes=minimized_cubes,
                )
            )
            for a, b in pairs:
                uf.union(a, b)
            if verbose:
                num_components_now = len({uf.find(i) for i in range(N)})
                print(
                    f"    [{examined}/{total_candidates}] ({e},{f}): SELECTED, witness_size={len(witness_qubits)}, "
                    f"{len(pairs)} pairs, components remaining={num_components_now} ({elapsed:.2f}s)"
                )
            num_components = len({uf.find(i) for i in range(N)})

    return MixerConstruction(
        n_qubits=n_qubits,
        spanning_trees=spanning_trees,
        terms=terms,
        fully_connected=(num_components == 1),
        dropped_candidates=dropped,
    )


def _swap_block(beta) -> QuantumCircuit:
    qc = QuantumCircuit(2, name="exchange_swap")
    qc.rxx(beta, 0, 1)
    qc.ryy(beta, 0, 1)
    return qc


def exchange_term_circuit(n_qubits: int, term: ExchangeTerm, beta) -> QuantumCircuit:
    """Builds one term's circuit from `term.minimized_cubes` (a
    Quine-McCluskey-minimized, EXACT cover of `term.valid_patterns` --
    see `_minimize_patterns`), not `valid_patterns` directly: each cube
    may leave some witness qubits as don't-cares, so it becomes a
    controlled gate on only the qubits that cube actually constrains --
    fewer and/or narrower controlled gates than one-gate-per-exact-
    pattern, for the identical validity function."""
    qc = QuantumCircuit(n_qubits, name=f"exchange({term.e},{term.f})")
    block = _swap_block(beta).to_gate()

    if not term.witness_qubits:
        qc.append(block, [term.e, term.f])
        return qc

    for cube in term.minimized_cubes:
        active = [(term.witness_qubits[i], v) for i, v in enumerate(cube) if v is not None]
        if not active:
            qc.append(block, [term.e, term.f])
            continue
        control_qubits = [q for q, _ in active]
        ctrl_state = sum(v << i for i, (_, v) in enumerate(active))
        controlled = block.control(len(control_qubits), ctrl_state=ctrl_state, annotated=False)
        qc.append(controlled, control_qubits + [term.e, term.f])
    return qc


def mixer_circuit(construction: MixerConstruction, beta=None) -> QuantumCircuit:
    if beta is None:
        beta = Parameter("beta")
    qc = QuantumCircuit(construction.n_qubits, name="matroid_basis_exchange_mixer")
    for term in construction.terms:
        qc.compose(exchange_term_circuit(construction.n_qubits, term, beta), inplace=True)
    return qc


def verify_no_leakage(construction: MixerConstruction, beta_value: float = 0.37, tol: float = 1e-9) -> bool:
    """Exact circuit-unitary check: every spanning tree's column must have
    zero amplitude on every non-spanning-tree computational basis state.
    Only feasible for small n_qubits (this is a correctness check on a
    handful of small instances, not part of the scaling measurement) --
    the full circuit's 2^n_qubits x 2^n_qubits unitary is what's built.
    For larger instances (e.g. real feeder topologies), use
    `verify_term_no_leakage` per term instead."""
    n = construction.n_qubits
    tree_set: Set[int] = set(construction.spanning_trees)
    qc = mixer_circuit(construction, beta=beta_value)
    U = Operator(qc).data
    for mask in construction.spanning_trees:
        col = U[:, mask]
        leaked = sum(abs(col[j]) ** 2 for j in range(len(col)) if j not in tree_set)
        if leaked > tol:
            return False
    return True


def verify_term_no_leakage(term: ExchangeTerm, beta_value: float = 0.37, tol: float = 1e-9) -> bool:
    """Exact correctness check for ONE term, independent of total qubit
    count. `exchange_term_circuit` only ever applies gates to `term.e`,
    `term.f`, and `term.witness_qubits` -- by construction it is identity
    on every other qubit -- so the term's full behavior can be verified
    exactly on just those qubits (remapped to a small local circuit),
    regardless of how large the real circuit is. This is what makes exact
    correctness verification tractable at real-feeder scale (e.g. 37
    qubits, where a full 2^37 x 2^37 unitary is not just slow but
    physically impossible to hold in memory): checking every term locally
    is equivalent to checking the full circuit, since the terms act on
    disjoint-in-effect qubit sets from the rest of the register.

    Enumerates every possible local state of (e, f, witness_qubits) and
    checks the gate confines amplitude to the right subspace: (a) if the
    state is a "trigger" (exactly one of e,f set) with a witness pattern
    in term.valid_patterns, amplitude must stay within {state,
    swapped-state} -- note this is a ROTATION at generic beta (a
    superposition of the two), NOT a full deterministic swap; requiring
    all amplitude to land exactly on the swapped state was an earlier bug
    in this function, caught by cross-checking against `verify_no_leakage`
    on small instances before trusting it at real-feeder scale where that
    full check is infeasible. (b) otherwise, it must act as identity (zero
    amplitude anywhere else).
    """
    local_qubits = list(term.witness_qubits) + [term.e, term.f]
    k = len(local_qubits)
    remap = {q: i for i, q in enumerate(local_qubits)}
    local_term = ExchangeTerm(
        e=remap[term.e],
        f=remap[term.f],
        witness_qubits=tuple(remap[q] for q in term.witness_qubits),
        valid_patterns=term.valid_patterns,
        minimized_cubes=term.minimized_cubes,
    )
    qc = exchange_term_circuit(k, local_term, beta_value)
    U = Operator(qc).data

    valid_set = set(local_term.valid_patterns)
    for state in range(2**k):
        bits = [(state >> i) & 1 for i in range(k)]
        e_val, f_val = bits[local_term.e], bits[local_term.f]
        witness_val = tuple(bits[w] for w in local_term.witness_qubits)
        should_fire = (e_val + f_val == 1) and (witness_val in valid_set)
        col = U[:, state]
        if should_fire:
            expected_target = state ^ (1 << local_term.e) ^ (1 << local_term.f)
            allowed = {state, expected_target}
        else:
            allowed = {state}
        leaked = sum(abs(col[j]) ** 2 for j in range(len(col)) if j not in allowed)
        if leaked > tol:
            return False
    return True


def verify_all_terms_no_leakage(construction: MixerConstruction, beta_value: float = 0.37) -> bool:
    return all(verify_term_no_leakage(t, beta_value) for t in construction.terms)
