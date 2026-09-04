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
    cost_alpha: float = 0.0,
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
    guaranteed optimal either -- a heuristic, stated as one.

    `cost_alpha=0.0` (default) reproduces the original behavior exactly:
    keep whichever witness minimizes leakage alone, regardless of size.
    Found (`docs/bounded-witness-mixer.md`'s circuit-cost investigation)
    to be a real problem: since adding a witness qubit can only weakly
    reduce leakage, never increase it, a leakage-only objective has no
    reason not to walk to `max_size` every time, even when a narrower
    witness would give nearly the same leakage at a fraction of the
    circuit cost (a generic k-controlled 2-qubit gate's cost grows
    steeply in k). `cost_alpha > 0` scores each candidate as
    `leak + cost_alpha * 2**size` instead of `leak` alone -- `2**size` as
    a cheap, size-only proxy for the controlled-gate cost that width
    actually drives, not a transpiled measurement (that would make the
    search itself expensive). `best_leak`/the returned `leak` are always
    the ACTUAL leakage of the selected witness, never the combined score
    -- the score only decides which candidate wins, it isn't reported as
    if it were the leakage rate."""
    best_witness: Tuple[int, ...] = ()
    _, best_leak = _majority_rule(trigger, validity, best_witness)
    best_score = best_leak  # size 0 costs nothing extra: 2**0 * cost_alpha would still be a flat offset, omitted since it's constant across all size-0 candidates (there's only one)

    for size in range(1, max_size + 1):
        if size > len(candidate_pool):
            break
        for _ in range(attempts_per_size):
            # int(q), not numpy's own int64 from rng.choice: Qiskit's
            # `.control(ctrl_state=...)` rejects numpy int types outright.
            witness = tuple(sorted(int(q) for q in rng.choice(candidate_pool, size=size, replace=False)))
            _, leak = _majority_rule(trigger, validity, witness)
            score = leak + cost_alpha * (2 ** size)
            if score < best_score - 1e-12:
                best_score = score
                best_leak = leak
                best_witness = witness

    valid_patterns, leak = _majority_rule(trigger, validity, best_witness)
    return best_witness, valid_patterns, leak


def _search_truncated_witness_adaptive(
    trigger: List[int],
    validity: Dict[int, bool],
    candidate_pool: Sequence[int],
    max_size: int,
    rng: np.random.Generator,
    attempts_per_size: int = 60,
    gain_price: float = 0.01,
) -> Tuple[Tuple[int, ...], Tuple[Tuple[int, ...], ...], float]:
    """Replaces `_search_truncated_witness`'s FIXED `cost_alpha` penalty
    with a PER-TERM one: still picks the size minimizing
    `leak(size) + alpha_term * 2**size`, but `alpha_term` is computed
    from THIS term's own achievable leak range (see `gain_price` below)
    instead of being one hand-picked number shared by every term in the
    construction.

    Why this, instead of picking a better fixed `cost_alpha`: the
    `cost_alpha` sweep (docs/bounded-witness-mixer.md's follow-up)
    found a single global coefficient doesn't fit terms of different
    difficulty -- on easy conditions most of the cost was pure slack
    (pushing alpha higher cut cost 10-100x with NO safety loss), while on
    hard conditions there was little slack to cut at all, and pushing
    alpha too hard on the SMALLEST instances reintroduced the same
    seed-fragility the cap=0/1 check already found. A fixed penalty can't
    tell those cases apart because it doesn't look at how much a given
    term is actually still improving.

    **Not stepwise-greedy stopping, and not a fixed-percentage-of-gain
    elbow either -- both were tried first.** Stepwise stopping (accept a
    size only if it beats the PREVIOUS size by enough) reproduces exactly
    the blind spot `_search_truncated_witness`'s own docstring already
    warns about for greedy forward selection: an interaction effect where
    sizes 1-2 individually buy little but size 3-4 TOGETHER buy a lot
    looks, to a stepwise rule, like "stop, nothing's helping." A
    fixed-percentage-of-the-full-curve elbow (accept the smallest size
    that captures e.g. 85% of the total leak reduction from size 0 to
    `max_size`) fixed that blind spot but didn't reproduce the slack the
    `cost_alpha` sweep found empirically -- on calibration instances it
    kept choosing witnesses nearly as wide as `cost_alpha=0.0` (no cost
    awareness at all), because when a term's leak curve is gradual rather
    than sharply elbowed, "capture most of the achievable gain" demands
    going most of the way to `max_size` regardless of how little that
    gain is worth in absolute terms.

    **A second failure mode found and fixed here, not just in the search
    objective's width term**: pushing `gain_price` (or `exact_search_max_size`
    down, routing more candidates through this search) hard enough, MOST
    terms independently landed on witness=() with `valid_patterns=()` --
    "never fire" -- because a term's (e,f) pair is often invalid for a
    strict MAJORITY of the trigger sample, so "never fire" already has
    low LEAK (misclassifying only the minority of states where it should
    have fired) at literally ZERO circuit cost, an unbeatable combination
    under the `leak + alpha*cost` score alone. On one real instance, 146
    of 151 terms collapsed this way -- a circuit that looks perfect
    (near-zero CX, unsafe_rate=0.0) for the wrong reason: it barely does
    anything, not because it preserves feasibility well.
    `mixer.build_matroid_mixer`'s own union-find only ever selects a
    candidate (e,f) pair as a term because it's confirmed to
    `connects_something` in the tree sample -- so a term resolving to
    "never fire" always defeats the reason it was selected at all, not
    just a leaky-but-acceptable compromise. Fixed by tracking, separately
    from the plain best-leak-at-each-size, the best leak found where
    `valid_patterns` is NON-empty ("active") -- and restricting the final
    size choice below to active options whenever at least one size offers
    one, falling back to an inert witness only if truly no size (0
    through `max_size`) can produce an active rule at all.

    **What this does, given that fix**: search every size from 0 to
    `max_size` INDEPENDENTLY first (fresh random-restart search per size
    -- no stepwise dependency, so the interaction-effect blind spot
    doesn't apply), tracking the best ACTIVE witness at each size, then
    pick the active size minimizing `leak(size) + alpha_term * 2**size`
    -- the SAME additive scoring `_search_truncated_witness` uses, but
    with a PER-TERM `alpha_term` instead of one fixed global value:
    `alpha_term = gain_price / max(total_achievable_gain, eps)`, where
    `total_achievable_gain = leak(0) - best_leak_found_at_any_size` is
    how much THIS term actually has to gain from growing its witness at
    all. A term with a lot of achievable gain (leak drops sharply with
    size -- the genuinely hard cases) gets a SMALL effective alpha, so
    width stays affordable; a term with little to gain (leak barely moves
    regardless of size -- the slack `cost_alpha` found on easy instances)
    gets a LARGE effective alpha, so it settles for a cheap, narrow
    witness instead of paying for width that isn't buying anything. This
    is what makes the per-term adaptation automatic instead of requiring
    a hand-picked global coefficient per condition."""
    sizes = list(range(0, min(max_size, len(candidate_pool)) + 1))
    best_leak_by_size: Dict[int, float] = {}
    best_witness_by_size: Dict[int, Tuple[int, ...]] = {}
    active_leak_by_size: Dict[int, float] = {}
    active_witness_by_size: Dict[int, Tuple[int, ...]] = {}

    patterns0, leak0 = _majority_rule(trigger, validity, ())
    best_leak_by_size[0] = leak0
    best_witness_by_size[0] = ()
    if patterns0:
        active_leak_by_size[0] = leak0
        active_witness_by_size[0] = ()

    for size in sizes[1:]:
        size_best_witness: Tuple[int, ...] = ()
        size_best_leak = leak0
        size_active_best_witness: Tuple[int, ...] = None
        size_active_best_leak = float("inf")
        for _ in range(attempts_per_size):
            witness = tuple(sorted(int(q) for q in rng.choice(candidate_pool, size=size, replace=False)))
            patterns, leak = _majority_rule(trigger, validity, witness)
            if leak < size_best_leak - 1e-12:
                size_best_leak = leak
                size_best_witness = witness
            if patterns and leak < size_active_best_leak - 1e-12:
                size_active_best_leak = leak
                size_active_best_witness = witness
        best_leak_by_size[size] = size_best_leak
        best_witness_by_size[size] = size_best_witness
        if size_active_best_witness is not None:
            active_leak_by_size[size] = size_active_best_leak
            active_witness_by_size[size] = size_active_best_witness

    best_overall_leak = min(best_leak_by_size.values())
    total_achievable_gain = leak0 - best_overall_leak
    alpha_term = gain_price / max(total_achievable_gain, 1e-6)

    if active_leak_by_size:
        # Prefer an active (non-inert) rule at every size that has one --
        # this is the fix: an inert "never fire" option is never allowed
        # to win purely on cost, only accepted below if NOTHING active
        # was found at any size (a rare, genuinely all-invalid case).
        best_leak_by_size, best_witness_by_size = active_leak_by_size, active_witness_by_size

    chosen_size = min(
        best_leak_by_size, key=lambda s: best_leak_by_size[s] + alpha_term * (2 ** s)
    )
    best_witness = best_witness_by_size[chosen_size]

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
    cost_alpha: float = 0.01,
    adaptive: bool = False,
    gain_price: float = 0.01,
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
    less and still catches the genuinely cheap exact cases.

    `adaptive=False` (default) uses the fixed-`cost_alpha` penalty.
    `adaptive=True` tries `_search_truncated_witness_adaptive` instead --
    a PER-TERM effective alpha (`gain_price` divided by how much THIS
    term actually has left to gain from a wider witness) that looked, on
    an initial sweep, like it dominated a fixed global `cost_alpha`. It
    doesn't, once fully validated: re-measured after the never-fire fix
    below (which the initial sweep predates), `adaptive=True` costs MORE
    than `cost_alpha=0.01` on most of the escalating-realism ladder
    (2-7x, e.g. `long_log`/`long_linear` at n_nodes=30-100), not less --
    the per-term formula turned out to be MORE prone to the never-fire
    shortcut than one shared coefficient, not less, because its
    `alpha_term` inflates whenever a term's OWN achievable gain looks
    small, which the never-fire option itself makes look smaller than it
    is. Kept available for comparison and further investigation, not
    because it's recommended -- `docs/bounded-witness-mixer.md` documents
    this as a real, instructive dead end, not a working improvement."""
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
                if adaptive:
                    approx_witness, approx_patterns, leak = _search_truncated_witness_adaptive(
                        trigger, validity, candidate_pool, max_witness_size, rng,
                        attempts_per_size=search_attempts_per_size,
                        gain_price=gain_price,
                    )
                else:
                    approx_witness, approx_patterns, leak = _search_truncated_witness(
                        trigger, validity, candidate_pool, max_witness_size, rng,
                        attempts_per_size=search_attempts_per_size, cost_alpha=cost_alpha,
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
