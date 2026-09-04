"""Bounded-witness (approximate) matroid mixer: trades EXACT constraint
preservation for a witness capped at a fixed small size, accepting some
leakage to infeasible states, in exchange for a construction that doesn't
simply DROP a candidate (e,f) pair the way `mixer.build_matroid_mixer`
does when `find_witness_set` finds nothing within its cap.

## Why this exists

`mixer.py`'s own docstring already states the bound its witness search
finds "does NOT hold in general" -- the real IEEE 33-bus feeder needed
zone decomposition specifically because the whole-graph witness
requirement reached into the dozens of qubits. `tie_density_sweep.py`
(this branch) finds a THIRD failure axis, distinct from both: even on the
repo's own nearest-tie synthetic family (short fundamental cycles by
construction), simply raising tie DENSITY -- not switching to the
long-range tie generator at all -- makes brute-force witness search fail
on a growing fraction of candidates. Zone decomposition is one fix
(shrink each subproblem so its witnesses stay small); this module is a
different one, for cases where decomposition either isn't natural for the
network topology or isn't wanted: accept a BOUNDED witness and the leakage
that comes with it, provided the mixer still starts from and stays
mostly on feasible states -- which only has to beat penalty-QAOA's
"uniform superposition, no structural awareness at all" baseline, not be
perfect.

## What a truncated term actually is

For each candidate (e, f) exchange pair, `mixer.find_witness_set` (brute
force only, small size cap) is tried FIRST, unchanged -- if an exact
witness exists within that cap, it's used, with zero leakage, exactly as
in the exact construction. Only when that fails does this module fall
back to an APPROXIMATE witness: random-restart search over witness qubits
(from the graph-derived cycle-union candidate pool, NOT the whole
register -- random qubits from the full register were confirmed
empirically to carry almost no predictive signal) that minimizes
MAJORITY-VOTE misclassification against the sample of trees actually
seen, capped at `max_witness_size`.

**Deliberately reuses `mixer.ExchangeTerm`/`MixerConstruction` and every
downstream consumer (`mixer_circuit`, `exchange_term_circuit`,
`verify_term_no_leakage`) completely unchanged.** A term doesn't "know"
its witness is approximate -- it's just an `ExchangeTerm` whose
`valid_patterns` happen to be a majority-vote rule instead of an
exactly-verified one. This means:
- Circuit-level correctness (does the compiled circuit implement the
  DECLARED valid_patterns exactly) is checked EXACTLY THE SAME way as
  always, via `verify_term_no_leakage` -- and still passes 100%, since
  nothing about circuit synthesis changed.
- What's NEW and must be tracked separately is whether the DECLARED rule
  itself is semantically correct against the true feasible set -- it
  isn't, by design, for approximate terms. `leakage_rate` on each
  `TruncatedTermInfo` records exactly how often, measured against the
  sample used to build it. `leakage_trace.py` (this branch) goes further,
  tracing whether that declared-rule imperfection actually costs
  probability mass on a REAL trajectory through the full term sequence,
  which is a different (and in practice smaller, but nonzero) question
  than the per-term rate alone.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence, Tuple

import numpy as np

from graphs import FeederGraph
from mixer import (
    ExchangeTerm,
    MixerConstruction,
    UnionFind,
    _minimize_patterns,
    find_witness_set,
    fundamental_cycle_structure,
)


@dataclass
class TruncatedTermInfo:
    e: int
    f: int
    is_exact: bool
    leakage_rate: float  # 0.0 for exact terms
    candidate_pool_size: int


def _majority_rule(
    trigger: List[int], validity: Dict[int, bool], witness: Tuple[int, ...]
) -> Tuple[Tuple[Tuple[int, ...], ...], float]:
    """Majority-vote valid_patterns for `witness`, plus the resulting
    misclassification (leakage) rate against `trigger`/`validity`."""
    pattern_votes: Dict[Tuple[int, ...], List[bool]] = {}
    for t in trigger:
        key = tuple((t >> q) & 1 for q in witness)
        pattern_votes.setdefault(key, []).append(validity[t])
    rule = {k: (sum(v) / len(v)) > 0.5 for k, v in pattern_votes.items()}
    wrong = 0
    for t in trigger:
        key = tuple((t >> q) & 1 for q in witness)
        if rule[key] != validity[t]:
            wrong += 1
    leakage_rate = wrong / len(trigger) if trigger else 0.0
    valid_patterns = tuple(k for k, v in rule.items() if v)
    return valid_patterns, leakage_rate


def _search_truncated_witness(
    trigger: List[int],
    validity: Dict[int, bool],
    candidate_pool: Sequence[int],
    max_size: int,
    rng: np.random.Generator,
    attempts_per_size: int = 60,
) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, ...], ...], float]:
    """Random-restart search: try `attempts_per_size` random subsets at
    each size up to `max_size`, keep whichever minimizes misclassification.

    GREEDY forward selection (add whichever single qubit helps most, then
    the next, ...) was tried first and consistently got stuck at witness
    size 0 -- no single candidate qubit showed ANY
    individual marginal benefit, even though random SUBSETS of size >= 4
    clearly did reduce misclassification. That is the signature of an
    interaction effect (validity depends on a COMBINATION of qubits, not
    any one of them additively) -- greedy provably cannot discover that.
    Random-restart search has no such blind spot, at the cost of not being
    guaranteed optimal either -- a heuristic, stated as one."""
    best_witness: Tuple[int, ...] = ()
    _, best_leak = _majority_rule(trigger, validity, best_witness)

    for size in range(1, max_size + 1):
        if size > len(candidate_pool):
            break
        for _ in range(attempts_per_size):
            # int(q), not numpy's own int64 from rng.choice: Qiskit's
            # `.control(ctrl_state=...)` rejects numpy int types outright.
            witness = tuple(sorted(int(q) for q in rng.choice(candidate_pool, size=size, replace=False)))
            _, leak = _majority_rule(trigger, validity, witness)
            if leak < best_leak - 1e-12:
                best_leak = leak
                best_witness = witness

    valid_patterns, leak = _majority_rule(trigger, validity, best_witness)
    return best_witness, valid_patterns, leak


@dataclass
class TruncatedMixerConstruction:
    construction: MixerConstruction
    term_info: List[TruncatedTermInfo]  # same order as construction.terms

    @property
    def n_exact_terms(self) -> int:
        return sum(1 for info in self.term_info if info.is_exact)

    @property
    def n_approximate_terms(self) -> int:
        return sum(1 for info in self.term_info if not info.is_exact)

    @property
    def mean_leakage_rate(self) -> float:
        approx = [info.leakage_rate for info in self.term_info if not info.is_exact]
        return float(np.mean(approx)) if approx else 0.0

    @property
    def max_leakage_rate(self) -> float:
        approx = [info.leakage_rate for info in self.term_info if not info.is_exact]
        return float(max(approx)) if approx else 0.0


def build_truncated_witness_mixer(
    graph: FeederGraph,
    trees: List[int],
    max_witness_size: int,
    exact_search_max_size: int = 2,
    search_attempts_per_size: int = 60,
    seed: int = 0,
) -> TruncatedMixerConstruction:
    """Same candidate-pair selection loop as `mixer.build_matroid_mixer`
    (union-find over which (e,f) pairs connect new components), but every
    candidate that connects something gets a term -- via an exact witness
    if `find_witness_set` (brute force only, NOT structural, see below)
    finds one within `exact_search_max_size`, else via the random-restart
    truncated approximate witness above, capped at `max_witness_size` --
    instead of being dropped when no exact witness exists.
    `dropped_candidates` on the resulting `MixerConstruction` is therefore
    always 0 by construction.

    `exact_search_max_size` defaults to 2, not `mixer.py`'s own default of
    4: brute-force search at size 4 gets expensive fast on the dense-tie
    instances this module targets (most candidates there have no small
    exact witness at all, so size-4 search fails, expensively, almost
    every time -- see `docs/bounded-witness-mixer.md`). Size 2 costs much
    less and still catches the genuinely cheap exact cases."""
    n_qubits = graph.n_edges
    index_of: Dict[int, int] = {mask: i for i, mask in enumerate(trees)}
    tree_set = set(trees)
    N = len(trees)
    cycle_of = fundamental_cycle_structure(graph)
    full_cycle_pool = tuple(sorted(set().union(*cycle_of.values()))) if cycle_of else ()
    rng = np.random.default_rng(seed)

    uf = UnionFind(N)
    num_components = N
    terms: List[ExchangeTerm] = []
    term_info: List[TruncatedTermInfo] = []

    for e in range(n_qubits):
        if num_components == 1:
            break
        for f in range(e + 1, n_qubits):
            if num_components == 1:
                break
            swap_mask = (1 << e) | (1 << f)
            pairs = []
            for mask in trees:
                other = mask ^ swap_mask
                if other in index_of and other != mask:
                    a, b = index_of[mask], index_of[other]
                    if a < b:
                        pairs.append((a, b))
            if not pairs:
                continue
            connects_something = any(uf.find(a) != uf.find(b) for a, b in pairs)
            if not connects_something:
                continue

            # prefer_structural left at its default (False) AND cycle_of=None
            # here deliberately: find_witness_set's structural path returns
            # the FULL cycle_of[e]|cycle_of[f] union regardless of `max_size`
            # (it is not capped at all) -- passing cycle_of would make every
            # candidate report an "exact" witness of whatever size that union
            # happens to be, defeating the size cap this function exists to
            # respect. Brute-force-only correctly returns None when no small
            # exact witness exists, falling through to the truncated search.
            witness_qubits, valid_patterns = find_witness_set(
                trees, tree_set, n_qubits, e, f, max_size=exact_search_max_size,
                cycle_of=None, prefer_structural=False,
            )

            if witness_qubits is not None:
                minimized_cubes = _minimize_patterns(len(witness_qubits), valid_patterns)
                terms.append(ExchangeTerm(e=e, f=f, witness_qubits=witness_qubits, valid_patterns=valid_patterns, minimized_cubes=minimized_cubes))
                term_info.append(TruncatedTermInfo(e=e, f=f, is_exact=True, leakage_rate=0.0, candidate_pool_size=len(witness_qubits)))
            else:
                trigger = [t for t in trees if bin(t & swap_mask).count("1") == 1]
                validity = {t: ((t ^ swap_mask) in tree_set) for t in trigger}
                # Local union (cycle_of[e]|cycle_of[f]) alone can be too
                # small to carry any signal at all -- fall back to the
                # full-graph cycle union as the search pool, same escalation
                # find_witness_set's own structural path uses.
                local_pool = tuple(sorted((cycle_of.get(e, set()) | cycle_of.get(f, set())) - {e, f}))
                candidate_pool = local_pool if len(local_pool) >= max_witness_size else tuple(
                    q for q in full_cycle_pool if q not in (e, f)
                )
                if not candidate_pool:
                    candidate_pool = tuple(q for q in range(n_qubits) if q not in (e, f))
                approx_witness, approx_patterns, leak = _search_truncated_witness(
                    trigger, validity, candidate_pool, max_witness_size, rng,
                    attempts_per_size=search_attempts_per_size,
                )
                minimized_cubes = _minimize_patterns(len(approx_witness), approx_patterns)
                terms.append(ExchangeTerm(e=e, f=f, witness_qubits=approx_witness, valid_patterns=approx_patterns, minimized_cubes=minimized_cubes))
                term_info.append(TruncatedTermInfo(e=e, f=f, is_exact=False, leakage_rate=leak, candidate_pool_size=len(candidate_pool)))

            for a, b in pairs:
                uf.union(a, b)
            num_components = len({uf.find(i) for i in range(N)})

    construction = MixerConstruction(
        n_qubits=n_qubits, spanning_trees=trees, terms=terms,
        fully_connected=(num_components == 1), dropped_candidates=0,
    )
    return TruncatedMixerConstruction(construction=construction, term_info=term_info)
