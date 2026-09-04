"""Uniform random spanning-tree generation via Wilson's algorithm
(loop-erased random walk), used to build mixers past the scale where
`graphs.enumerate_spanning_trees` (brute force over C(n_edges, k_ties)
subsets) is tractable.

## Why this exists

`enumerate_spanning_trees`'s own docstring notes its cost is
"polynomial in n_nodes for the fixed, small k_ties this study uses" --
true for the sparse nearest-tie family, but the TRUE number of spanning
trees a graph has (via Kirchhoff's theorem, the Laplacian cofactor -- no
enumeration needed to COUNT them) grows explosively once tie density
increases (see `tie_density_sweep.py`, and `docs/bounded-witness-mixer.md`
for the actual counts observed here). No enumeration algorithm, however
efficient, changes this -- it is a property of the OUTPUT size, not the
algorithm.

Wilson's algorithm generates an EXACTLY uniformly random spanning tree in
expected polynomial time, independent of how many spanning trees the graph
has (Propp & Wilson, 1996), implemented here in its standard
"next-pointer" form (loop erasure happens implicitly via pointer
overwriting). Not trusted blindly: `verify_random_trees.py` checks every
generated sample is a genuine spanning tree (networkx `is_tree`,
independent of this module's own bitmask bookkeeping) and that, on a
small enumerable instance, the empirical sampling distribution is close to
uniform (cross-checked against `graphs.enumerate_spanning_trees`) before
this is trusted at any scale where that cross-check is impossible.

## What this does NOT provide

A random sample of feasible trees fed into `mixer.build_matroid_mixer`
(which only ever consumes whatever tree list it's given) gives
STATISTICAL, not exhaustive, confidence that the discovered exchange terms
are semantically correct across the true (possibly astronomically large)
feasible set -- weaker than this repo's usual exhaustive
`enumerate_spanning_trees` cross-check. What stays exact regardless of
scale: once a term's witness/valid_patterns are discovered (by whatever
method), `mixer.verify_term_no_leakage` confirms the COMPILED CIRCUIT
correctly implements exactly that declared witness condition -- a
circuit-correctness guarantee, not a claim that the declared condition
itself is semantically complete against the true feasible set. See
`docs/bounded-witness-mixer.md` for how this distinction plays out for the
approximate mixer in `truncated_mixer.py`.
"""

from __future__ import annotations

from typing import List, Set

import numpy as np

from graphs import FeederGraph


def random_spanning_tree(graph: FeederGraph, rng: np.random.Generator) -> int:
    """One uniformly random spanning tree via Wilson's algorithm, as an
    edge-index bitmask (same edge-i-is-qubit-i convention as
    `graphs.enumerate_spanning_trees`)."""
    n = graph.n_nodes
    adj: List[List[tuple]] = [[] for _ in range(n)]
    for i, (u, v) in enumerate(graph.edges):
        adj[u].append((v, i))
        adj[v].append((u, i))

    root = 0
    in_tree = np.zeros(n, dtype=bool)
    in_tree[root] = True
    next_node = np.full(n, -1, dtype=np.int64)
    next_edge = np.full(n, -1, dtype=np.int64)

    for start in range(n):
        u = start
        while not in_tree[u]:
            neighbors = adj[u]
            choice = neighbors[rng.integers(len(neighbors))]
            next_node[u], next_edge[u] = choice
            u = next_node[u]
        u = start
        while not in_tree[u]:
            in_tree[u] = True
            u = next_node[u]

    mask = 0
    for v in range(n):
        if v != root:
            mask |= 1 << int(next_edge[v])
    return mask


def sample_spanning_trees(graph: FeederGraph, n_samples: int, seed: int) -> List[int]:
    """`n_samples` i.i.d. uniform spanning trees, deduplicated (duplicates
    are expected, and increasingly likely as `n_samples` approaches the
    graph's true, possibly-huge, tree count -- the unique set is what a
    downstream witness search actually needs)."""
    rng = np.random.default_rng(seed)
    seen: Set[int] = set()
    for _ in range(n_samples):
        seen.add(random_spanning_tree(graph, rng))
    return list(seen)


def random_walk_exchange_sample(graph: FeederGraph, n_steps: int, seed: int, start_mask: int = None) -> List[int]:
    """A large i.i.d. random sample turns out NOT to be what
    `mixer.build_matroid_mixer` needs once the true tree count is large:
    two INDEPENDENT uniform random trees are, at astronomical tree counts,
    almost never a single edge-swap apart -- the exchange graph is
    astronomically sparse relative to the tree space, so union-find over
    an i.i.d. sample essentially never finds an edge to select as a term.
    What `build_matroid_mixer` actually needs is a sample rich in
    ALREADY-ADJACENT pairs.

    This walks the exchange graph directly instead: from a starting tree,
    repeatedly pick a random currently-active tree edge `e`, compute the
    cut its removal creates (BFS on the tree minus `e`), pick a random
    candidate `f` crossing that cut (a graph-theoretically NECESSARY
    condition for T-e+f to be a tree -- see `mixer.py`'s basis-exchange
    docstring), and record every tree visited. By construction, every
    visited tree is reachable from the start via genuine single-edge
    exchanges, so the resulting sample is dense in exchange-adjacent
    pairs. The tradeoff, stated plainly: this explores (a possibly partial
    view of) the connected component reachable from `start_mask`, not a
    uniform sample of the whole feasible set -- `build_matroid_mixer`'s
    resulting `fully_connected` flag reports connectivity honestly against
    whatever tree list it's actually given, so this is visible, not
    hidden, downstream."""
    rng = np.random.default_rng(seed)
    n = graph.n_nodes
    m = graph.n_edges
    edges = graph.edges

    if start_mask is None:
        start_mask = random_spanning_tree(graph, rng)

    current = start_mask
    visited: Set[int] = {current}

    for _ in range(n_steps):
        tree_edges = [i for i in range(m) if (current >> i) & 1]
        if not tree_edges:
            break
        e = int(tree_edges[rng.integers(len(tree_edges))])

        adj: List[List[int]] = [[] for _ in range(n)]
        for i in tree_edges:
            if i == e:
                continue
            u, v = edges[i]
            adj[u].append(v)
            adj[v].append(u)
        eu, ev = edges[e]
        comp = {eu}
        stack = [eu]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in comp:
                    comp.add(y)
                    stack.append(y)

        candidates = [
            i for i in range(m)
            if i != e and not ((current >> i) & 1) and (edges[i][0] in comp) != (edges[i][1] in comp)
        ]
        if not candidates:
            continue
        f = int(candidates[rng.integers(len(candidates))])
        current = current ^ (1 << e) ^ (1 << f)
        visited.add(current)

    return list(visited)
