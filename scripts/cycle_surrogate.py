"""A problem-informed, LOCAL cost surrogate for radial reconfiguration:
lossless real-power flow, linearized along fundamental cycles.

## The obstacle this gets around

The physical objective for reconfiguration -- line loss `sum_e r_e f_e^2`,
or a congestion surrogate `sum_e (f_e / cap_e)^2` -- is a function of the
flow `f_e`, and on a radial network `f_e` is the total demand DOWNSTREAM
of `e` in whatever tree the switches currently form. "Downstream in the
current tree" is a global property of the configuration, not a function
of a few nearby switch bits. There is no low-order polynomial in the
switch variables that expresses it, which is why question 1's mixer has
never had a phase separator to pair with (README, "What is and isn't
shown"). The classical MILP for this objective gets around it with
auxiliary flow variables and conservation constraints; in a circuit that
means extra qubits and a new feasibility constraint the mixer does not
preserve.

## The surrogate

A real feeder sits at a base tree `T0` (backbone closed, ties open).
Each tie `t` has a fundamental cycle `C_t` = the `T0` path between its
endpoints, plus `t`. A single exchange on that cycle -- close `t`, open
some `e` in `C_t` -- reroutes the subtree cut off by `e` through `t`.

The fact everything rests on, CHECKED NUMERICALLY in
`verify_cycle_surrogate.py` rather than assumed: **a single exchange
changes the flow only on the edges of its own cycle.** Every edge off
`C_t` keeps its downstream set; every edge above the cycle's junction
with the root carries the same total. So `delta_f_e(t, e_open)` is a
fixed number per (tie, opened edge, affected edge), computed once from
`T0`, real loads and real topology.

Parameterize a configuration by one choice per tie -- which edge of
`C_t` is open (choosing `t` itself means "tie stays open", the base
state for that cycle). Then

    f_e(choices) ~= f_e^0 + sum_t delta_f_e(t, choice_t)

Flows are SIGNED against a fixed reference orientation per edge (base
tree's parent->child direction for tree edges, low->high bus for ties),
because a rerouted segment genuinely reverses direction and magnitudes
do not superpose. Squaring happens in the objective, after superposition.

Exact for edge-disjoint cycles and for any single exchange; first-order
for overlapping cycles, whose error `verify_cycle_surrogate.py` measures
on every double exchange rather than estimating.

## What this buys

The objective becomes QUADRATIC in one-hot choice variables: Z and ZZ
terms across registers, coupled only where two cycles share an edge.
That is the hard-instance Hamiltonian's shape (`partition_mixer.py`),
reappearing -- there the coupling was "same bucket", here it is "shared
cycle edge". Same local structure, same per-register mixer.

Coverage caveat, stated up front: one-open-edge-per-fundamental-cycle
covers ALL spanning trees only when cycles are edge-disjoint. With
overlap it is a subset (some trees open a shared edge and use the freed
choice elsewhere), and some joint choices -- both cycles open their
shared edge -- disconnect the network. `verify_cycle_surrogate.py`
reports the overlap structure so the size of that caveat is a measured
number, not a guess.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import networkx as nx
import pandapower.networks as pn

from graphs import FeederGraph


@dataclass(frozen=True)
class PhysicalFeeder:
    graph: FeederGraph
    root: int
    tie_idx: tuple[int, ...]        # edge indices that are ties (open in the base configuration)
    tree_idx: tuple[int, ...]       # edge indices of the base tree T0
    load_mw: tuple[float, ...]      # per node
    r_ohm: tuple[float, ...]        # per edge
    orient: tuple[tuple[int, int], ...]  # per edge: (from, to) reference direction for signed flow


def load_ieee33_physical(root: int = 0) -> PhysicalFeeder:
    """Same edge ordering as `real_feeders.load_ieee33` (sorted, ties mixed
    in), plus the physical data that loader deliberately leaves out. Ties
    are recovered from pandapower's `in_service` flag, not inferred."""
    net = pn.case33bw()
    n_nodes = len(net.bus)

    rows: dict[tuple[int, int], dict] = {}
    for _, row in net.line.iterrows():
        e = tuple(sorted((int(row["from_bus"]), int(row["to_bus"]))))
        rows[e] = {"tie": not bool(row["in_service"]), "r": float(row["r_ohm_per_km"]) * float(row["length_km"])}

    edges = tuple(sorted(rows))
    tie_idx = tuple(i for i, e in enumerate(edges) if rows[e]["tie"])
    tree_idx = tuple(i for i, e in enumerate(edges) if not rows[e]["tie"])
    r_ohm = tuple(rows[e]["r"] for e in edges)

    load = [0.0] * n_nodes
    for _, row in net.load.iterrows():
        load[int(row["bus"])] += float(row["p_mw"])

    graph = FeederGraph(n_nodes=n_nodes, k_ties=len(tie_idx), edges=edges)

    # Reference orientation: base-tree parent->child for tree edges (so
    # base flows are all >= 0), low->high bus for ties.
    t0 = nx.Graph()
    t0.add_nodes_from(range(n_nodes))
    t0.add_edges_from(edges[i] for i in tree_idx)
    parent = dict(nx.bfs_predecessors(t0, root))
    orient = []
    for i, (u, v) in enumerate(edges):
        if i in set(tree_idx):
            orient.append((v, u) if parent.get(u) == v else (u, v))
        else:
            orient.append((u, v))
    return PhysicalFeeder(graph, root, tie_idx, tree_idx, tuple(load), r_ohm, tuple(orient))


def signed_flows(feeder: PhysicalFeeder, tree_edges: set[int]) -> dict[int, float]:
    """Lossless radial flow on the tree `tree_edges` (edge indices): the
    signed real power crossing each tree edge in its reference direction
    = (+/-) total demand of the side away from the root. Edges not in the
    tree carry 0. Post-order accumulation from the root; no iteration,
    because a radial network's flow is determined by its topology."""
    g = feeder.graph
    adj: dict[int, list[tuple[int, int]]] = {n: [] for n in range(g.n_nodes)}
    for i in tree_edges:
        u, v = g.edges[i]
        adj[u].append((v, i))
        adj[v].append((u, i))

    parent_edge: dict[int, int | None] = {feeder.root: None}
    order = [feeder.root]
    stack = [feeder.root]
    while stack:
        u = stack.pop()
        for v, ei in adj[u]:
            if v not in parent_edge:
                parent_edge[v] = ei
                order.append(v)
                stack.append(v)
    if len(order) != g.n_nodes:
        raise ValueError("tree_edges do not span the network")

    subtree = list(feeder.load_mw)
    flows = {i: 0.0 for i in range(g.n_edges)}
    for node in reversed(order):
        ei = parent_edge[node]
        if ei is None:
            continue
        u, v = g.edges[ei]
        par = v if u == node else u
        # flow direction is par -> node; positive iff that matches the reference orientation
        sign = 1.0 if feeder.orient[ei] == (par, node) else -1.0
        flows[ei] = sign * subtree[node]
        subtree[par] += subtree[node]
    return flows


def fundamental_cycle(feeder: PhysicalFeeder, tie: int) -> tuple[int, ...]:
    """Edge indices of tie `tie`'s fundamental cycle w.r.t. the base tree:
    the T0 path between its endpoints, plus the tie itself."""
    g = feeder.graph
    t0 = nx.Graph()
    t0.add_nodes_from(range(g.n_nodes))
    index_of = {g.edges[i]: i for i in feeder.tree_idx}
    t0.add_edges_from(index_of)
    u, v = g.edges[tie]
    path = nx.shortest_path(t0, u, v)
    cycle = [index_of[tuple(sorted((path[i], path[i + 1])))] for i in range(len(path) - 1)]
    return tuple(cycle) + (tie,)


def exchange_tree(feeder: PhysicalFeeder, choices: dict[int, int]) -> set[int]:
    """Tree edge set for a set of choices {tie: opened edge}. Choosing the
    tie itself leaves that cycle in its base state. Not guaranteed to be
    a spanning tree when cycles overlap -- callers check via
    `signed_flows`, which raises if it does not span."""
    tree = set(feeder.tree_idx)
    for tie, opened in choices.items():
        if opened == tie:
            continue
        tree.add(tie)
        tree.discard(opened)
    return tree


@dataclass
class Surrogate:
    feeder: PhysicalFeeder
    cycles: dict[int, tuple[int, ...]]                 # tie -> cycle edges (tie included, last)
    base_flow: dict[int, float]
    delta: dict[int, dict[int, dict[int, float]]]      # tie -> opened edge -> {edge: delta flow}

    def flows(self, choices: dict[int, int]) -> dict[int, float]:
        f = dict(self.base_flow)
        for tie, opened in choices.items():
            if opened == tie:
                continue
            for e, d in self.delta[tie][opened].items():
                f[e] += d
        return f

    def loss(self, flows: dict[int, float]) -> float:
        return sum(self.feeder.r_ohm[e] * f * f for e, f in flows.items())

    def overlapping_pairs(self) -> list[tuple[int, int, tuple[int, ...]]]:
        out = []
        for t1, t2 in combinations(sorted(self.cycles), 2):
            shared = tuple(sorted(set(self.cycles[t1]) & set(self.cycles[t2])))
            if shared:
                out.append((t1, t2, shared))
        return out


class SecondOrderLoss:
    """Loss expanded to second order in "which exchanges happened", with
    every single and double exchange tabulated EXACTLY:

        L(choices) ~= L0 + sum_t a_t(o_t) + sum_{t1<t2} b_{t1 t2}(o_1, o_2)

    a_t(o)      = L(single exchange t->o) - L0
    b(o1, o2)   = L(double) - L0 - a_t1(o1) - a_t2(o2)     (0 for disjoint cycles)

    Why loss and not flow: `verify_cycle_surrogate.py` found the
    first-order FLOW surrogate off by up to 346% of true loss on IEEE33's
    overlapping-cycle double exchanges (8 of 10 tie pairs overlap there).
    Expanding the loss instead makes every double exchange exact by
    construction and moves the approximation to 3+ simultaneous
    exchanges -- which `verify_cycle_surrogate_full.py` measures over the
    ENTIRE configuration space, since 5 ties make it enumerable.

    Still a pairwise objective in one-hot variables: a_t(o) is a Z term
    on qubit (t, o); b(o1, o2) is a ZZ term between (t1, o1) and (t2, o2).
    Joint choices that disconnect the network get `infeasible_penalty` as
    their b -- a penalty term, unavoidable because one-open-edge-per-cycle
    admits those combinations when cycles overlap; how many such pairs
    exist is reported, not hidden."""

    def __init__(self, s: "Surrogate", infeasible_penalty: float):
        self.s = s
        self.L0 = s.loss(s.base_flow)
        self.a: dict[int, dict[int, float]] = {}
        self.b: dict[tuple[int, int], dict[tuple[int, int], float]] = {}
        self.n_infeasible_pairs = 0
        for t in s.cycles:
            self.a[t] = {t: 0.0}
            for o in s.delta[t]:
                self.a[t][o] = s.loss(signed_flows(s.feeder, exchange_tree(s.feeder, {t: o}))) - self.L0
        for t1, t2 in combinations(sorted(s.cycles), 2):
            tab: dict[tuple[int, int], float] = {}
            for o1 in s.delta[t1]:
                for o2 in s.delta[t2]:
                    try:
                        L = s.loss(signed_flows(s.feeder, exchange_tree(s.feeder, {t1: o1, t2: o2})))
                        tab[(o1, o2)] = L - self.L0 - self.a[t1][o1] - self.a[t2][o2]
                    except ValueError:
                        tab[(o1, o2)] = infeasible_penalty
                        self.n_infeasible_pairs += 1
            self.b[(t1, t2)] = tab

    def loss(self, choices: dict[int, int]) -> float:
        L = self.L0 + sum(self.a[t][o] for t, o in choices.items())
        active = sorted((t, o) for t, o in choices.items() if o != t)
        for (t1, o1), (t2, o2) in combinations(active, 2):
            L += self.b[(t1, t2)].get((o1, o2), 0.0)
        return L

    def n_pairwise_terms(self) -> int:
        return sum(1 for tab in self.b.values() for v in tab.values() if v != 0.0)


def build_surrogate(feeder: PhysicalFeeder) -> Surrogate:
    base = signed_flows(feeder, set(feeder.tree_idx))
    cycles = {t: fundamental_cycle(feeder, t) for t in feeder.tie_idx}
    delta: dict[int, dict[int, dict[int, float]]] = {}
    for t, cyc in cycles.items():
        delta[t] = {}
        for opened in cyc:
            if opened == t:
                continue
            f = signed_flows(feeder, exchange_tree(feeder, {t: opened}))
            delta[t][opened] = {e: f[e] - base[e] for e in range(feeder.graph.n_edges) if f[e] != base[e]}
    return Surrogate(feeder, cycles, base, delta)
