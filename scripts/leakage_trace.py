"""Exact diagnostic tooling for `truncated_mixer.py`'s bounded-witness
construction: WHY does a term's classically-measured leakage rate
(majority-vote misclassification against a uniform sample of trigger
states, see `truncated_mixer.py`) sometimes not manifest as observed
infeasibility, and sometimes does?

## The mechanism, established here empirically and exactly (not sampled)

A term's declared leakage rate is a property of the ABSTRACT validity
function, measured uniformly over a sample. Whether it actually costs
probability mass in a REAL run depends on whether the SPECIFIC quantum
trajectory (one starting tree, evolving through this exact sequence of
terms) ever accumulates non-negligible amplitude on one of that term's
"danger states" -- a trigger state where the declared rule says fire, but
the true exchange is invalid. `trace_danger_mass` computes this EXACTLY
(full-statevector simulation, tractable only up to a few hundred basis
states), not by sampling.

This is a genuine, starting-tree-DEPENDENT property (confirmed directly on
this repo's own graph family below), not something that "never manifests"
in general. Since it is checkable exactly and cheaply for small instances
(one statevector simulation per candidate starting tree), a starting tree
can be VERIFIED safe before use, rather than hoped to be --
`is_safe_starting_tree` does exactly that check. `trace_danger_mass_sparse`
extends the same exact computation past small instances (real-scale
networks, hundreds of qubits) by tracking only the small subset of basis
states that ever gain amplitude, instead of the full `2**n_qubits`-
dimensional statevector -- see its own docstring for why that subset stays
small, and why that's what makes real-scale tracing possible at all.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from qiskit import QuantumCircuit
from qiskit.quantum_info import Statevector

from graphs import FeederGraph
from mixer import ExchangeTerm, MixerConstruction, exchange_term_circuit


def _is_genuine_tree(graph: FeederGraph, mask: int) -> bool:
    """Union-find tree check, same lightweight approach as
    `graphs.enumerate_spanning_trees` (avoids constructing a networkx
    Graph per candidate, which dominates runtime once this is called
    thousands of times per trace, as `trace_danger_mass_sparse` does at
    real scale)."""
    n = graph.n_nodes
    edges = graph.edges
    bits = [i for i in range(graph.n_edges) if (mask >> i) & 1]
    if len(bits) != n - 1:
        return False
    parent = list(range(n))

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in bits:
        u, v = edges[i]
        ru, rv = find(u), find(v)
        if ru == rv:
            return False
        parent[ru] = rv
    return True


@dataclass
class TermDangerInfo:
    term_index: int
    e: int
    f: int
    danger_mass: float  # probability mass, right before this term's gate, on states where its declared rule misfires
    feasible_mass_after: float  # exact total probability mass on genuine spanning trees, right after this term's gate


def trace_danger_mass(
    graph: FeederGraph, construction: MixerConstruction, start_tree: int, beta: float
) -> List[TermDangerInfo]:
    """Exact (statevector, not sampled) trace of one starting tree through
    `construction`'s full term sequence at a fixed `beta`. Small-instance
    only -- cost is O(2**n_qubits) per term."""
    n = construction.n_qubits
    qc = QuantumCircuit(n)
    for i in range(n):
        if (start_tree >> i) & 1:
            qc.x(i)
    sv = Statevector(qc)

    trace: List[TermDangerInfo] = []
    for k, term in enumerate(construction.terms):
        swap_mask = (1 << term.e) | (1 << term.f)
        valid_set = set(term.valid_patterns)
        witness = term.witness_qubits

        danger_mass = 0.0
        for x in range(2**n):
            p = abs(sv.data[x]) ** 2
            if p < 1e-12 or bin(x & swap_mask).count("1") != 1:
                continue
            witness_pattern = tuple((x >> q) & 1 for q in witness)
            if witness_pattern in valid_set and not _is_genuine_tree(graph, x ^ swap_mask):
                danger_mass += p

        term_qc = exchange_term_circuit(n, term, beta)
        sv = sv.evolve(term_qc)
        feasible_mass = sum(abs(sv.data[x]) ** 2 for x in range(2**n) if _is_genuine_tree(graph, x))
        trace.append(TermDangerInfo(term_index=k, e=term.e, f=term.f, danger_mass=danger_mass, feasible_mass_after=feasible_mass))

    return trace


def trace_danger_mass_sparse(
    graph: FeederGraph,
    construction: MixerConstruction,
    start_tree: int,
    beta: float,
    max_states: int = 200_000,
) -> List[TermDangerInfo]:
    """Same exact computation as `trace_danger_mass`, but WITHOUT
    representing the full `2**n_qubits`-dimensional statevector -- tracks
    only basis states with non-negligible amplitude, in a dict. Starting
    from ONE basis state and applying only LOCAL (2-qubit,
    witness-conditioned) gates, the number of states that ever gain
    amplitude is bounded by how many terms actually fire on already-active
    states, not by the total Hilbert space size -- this is what makes
    real-scale (n_qubits in the hundreds) tracing possible where full
    `Statevector` construction is not.

    The 2-level gate action used per firing term -- exactly
    `[[cos(beta), -i*sin(beta)], [-i*sin(beta), cos(beta)]]` on the
    trigger pair `(x, x^swap_mask)` -- is the exact restriction of
    `mixer._swap_block`'s (Rxx;Ryy) 4x4 unitary to the single-excitation
    subspace, confirmed directly against Qiskit's own `Operator` on that
    block at several angles (`verify_leakage_trace.py`), not assumed from
    the gate's name. That script additionally cross-checks this function's
    OUTPUT against `trace_danger_mass`'s (full statevector, independently
    exact) on the same small instance before this is trusted at any scale
    where that cross-check is impossible.

    Raises `RuntimeError` if the active state count ever exceeds
    `max_states` -- better to fail loudly than silently produce a
    truncated, wrong trace."""
    cos_b, sin_b = math.cos(beta), math.sin(beta)
    state: Dict[int, complex] = {start_tree: complex(1.0)}

    trace: List[TermDangerInfo] = []
    for k, term in enumerate(construction.terms):
        swap_mask = (1 << term.e) | (1 << term.f)
        valid_set = set(term.valid_patterns)
        witness = term.witness_qubits

        danger_mass = 0.0
        for x, amp in state.items():
            p = abs(amp) ** 2
            if p < 1e-15 or bin(x & swap_mask).count("1") != 1:
                continue
            witness_pattern = tuple((x >> q) & 1 for q in witness)
            if witness_pattern in valid_set and not _is_genuine_tree(graph, x ^ swap_mask):
                danger_mass += p

        new_state: Dict[int, complex] = {}
        handled: set = set()
        for x, amp in state.items():
            if x in handled:
                continue
            fires = False
            if bin(x & swap_mask).count("1") == 1:
                witness_pattern = tuple((x >> q) & 1 for q in witness)
                fires = witness_pattern in valid_set
            if not fires:
                new_state[x] = new_state.get(x, 0j) + amp
                handled.add(x)
                continue
            partner = x ^ swap_mask
            partner_amp = state.get(partner, 0j)
            new_state[x] = new_state.get(x, 0j) + cos_b * amp - 1j * sin_b * partner_amp
            new_state[partner] = new_state.get(partner, 0j) - 1j * sin_b * amp + cos_b * partner_amp
            handled.add(x)
            handled.add(partner)

        state = {m: a for m, a in new_state.items() if abs(a) ** 2 > 1e-14}
        if len(state) > max_states:
            raise RuntimeError(
                f"active state count exceeded max_states={max_states} at term {k} -- "
                f"the reachable subspace from this starting tree is larger than expected"
            )
        feasible_mass = sum(abs(a) ** 2 for m, a in state.items() if _is_genuine_tree(graph, m))
        trace.append(
            TermDangerInfo(term_index=k, e=term.e, f=term.f, danger_mass=danger_mass, feasible_mass_after=feasible_mass)
        )

    return trace


def final_feasible_mass(
    graph: FeederGraph, construction: MixerConstruction, start_tree: int, beta: float, sparse: bool = True
) -> float:
    """`sparse=True` (default) uses `trace_danger_mass_sparse`, validated
    exactly equivalent to the dense method (`verify_leakage_trace.py`) but
    without its O(2**n_qubits) ceiling -- pass `sparse=False` only when
    cross-checking against the dense method itself; real instances
    (n_qubits beyond ~25) will simply fail with a `ValueError` from
    `Statevector` otherwise."""
    tracer = trace_danger_mass_sparse if sparse else trace_danger_mass
    trace = tracer(graph, construction, start_tree, beta)
    return trace[-1].feasible_mass_after if trace else float(_is_genuine_tree(graph, start_tree))


def is_safe_starting_tree(
    graph: FeederGraph, construction: MixerConstruction, start_tree: int, beta: float, tol: float = 1e-6,
    sparse: bool = True,
) -> bool:
    """Checks BEFORE using a starting tree whether its trajectory through
    `construction` ever touches a danger state, instead of hoping.
    `sparse=True` (default): see `final_feasible_mass`'s docstring."""
    tracer = trace_danger_mass_sparse if sparse else trace_danger_mass
    trace = tracer(graph, construction, start_tree, beta)
    return all(info.danger_mass < tol for info in trace)


def survey_starting_trees(
    graph: FeederGraph, construction: MixerConstruction, beta: float, sparse: bool = True
) -> List[Tuple[int, float]]:
    """Exact final feasible mass for EVERY tree in `construction.spanning_trees`
    as a starting point. `sparse=True` (default): see `final_feasible_mass`'s
    docstring -- required past small instances."""
    return [
        (start_tree, final_feasible_mass(graph, construction, start_tree, beta, sparse=sparse))
        for start_tree in construction.spanning_trees
    ]
