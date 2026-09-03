"""Synthetic sparse/near-planar graph generation and spanning-tree
enumeration.

Graph model: place `n_nodes` points at random in the unit square, connect
them with a Euclidean minimum spanning tree (a natural, defensible model of
a distribution feeder's primary backbone -- minimising conductor length is
a real design objective), then add a small, FIXED number `k` of "tie"
edges connecting geographically nearby node pairs not already in the tree
(modelling the normally-open tie switches real feeders use for
reconfiguration/redundancy). This keeps the graph sparse and near-planar by
construction: |E| = n_nodes - 1 + k, with k held constant across a size
sweep.

This is a synthetic model, not derived from real feeder topology data --
stated plainly in methodology.md and here. It is chosen because it is
reproducible, controllable (k fixed lets |E| grow linearly and
predictably with n_nodes), and structurally similar to real feeders
(sparse, spatial, mostly-tree-plus-a-few-ties), not because it has been
validated against real distribution network data.

A second generator, `generate_feeder_graph_long_range_ties`, uses the
same MST backbone but picks tie edges to *maximize* tree-path length
instead of minimizing it -- reproducing the failure mode found on real
topology (see docs/circuit-validity.md) at controllable, arbitrary size.
Everything above about `generate_feeder_graph` (nearest-tie) is unchanged
by its existence; the two are used for different purposes throughout this
repo, not interchangeably.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from typing import List, Tuple

import networkx as nx
import numpy as np


@dataclass
class FeederGraph:
    n_nodes: int
    k_ties: int
    edges: Tuple[Tuple[int, int], ...]  # canonical (u<v) order, fixed indexing = qubit index

    @property
    def n_edges(self) -> int:
        return len(self.edges)

    def networkx_graph(self) -> nx.Graph:
        g = nx.Graph()
        g.add_nodes_from(range(self.n_nodes))
        g.add_edges_from(self.edges)
        return g


def generate_feeder_graph(n_nodes: int, k_ties: int, seed: int) -> FeederGraph:
    if n_nodes < 2:
        raise ValueError("need at least 2 nodes")
    rng = np.random.default_rng(seed)
    points = rng.random((n_nodes, 2))

    complete = nx.Graph()
    complete.add_nodes_from(range(n_nodes))
    for u, v in itertools.combinations(range(n_nodes), 2):
        dist = float(np.linalg.norm(points[u] - points[v]))
        complete.add_edge(u, v, weight=dist)

    mst = nx.minimum_spanning_tree(complete, weight="weight")
    tree_edges = {tuple(sorted(e)) for e in mst.edges()}

    non_tree_by_distance = sorted(
        (e for e in complete.edges() if tuple(sorted(e)) not in tree_edges),
        key=lambda e: complete.edges[e]["weight"],
    )
    tie_edges = [tuple(sorted(e)) for e in non_tree_by_distance[:k_ties]]

    all_edges = sorted(tree_edges | set(tie_edges))
    return FeederGraph(n_nodes=n_nodes, k_ties=len(tie_edges), edges=tuple(all_edges))


def generate_feeder_graph_long_range_ties(n_nodes: int, k_ties: int, seed: int) -> FeederGraph:
    """Like `generate_feeder_graph`, but chooses tie edges to MAXIMIZE
    tree-path length between endpoints (not minimize Euclidean distance).

    `generate_feeder_graph`'s nearest-tie choice keeps fundamental cycles
    short (2-3 tree edges, measured in `docs/circuit-validity.md`), which
    is exactly why the whole-graph witness-conditioned mixer stays cheap
    on it -- but real published tie switches deliberately link *distant*
    parts of a feeder for reconfiguration redundancy, and that long-range
    bias is what was found (on the one real topology available, IEEE
    33-bus) to break the whole-graph mixer construction (see
    `docs/circuit-validity.md`'s real-feeder section). This generator
    reproduces that long-range-tie failure mode at any controllable
    network size, so the zone-decomposition fix's scaling behavior can be
    tested against a family of instances rather than the one fixed real
    size available -- see `run_decomposition_scaling_study.py`.
    """
    if n_nodes < 2:
        raise ValueError("need at least 2 nodes")
    rng = np.random.default_rng(seed)
    points = rng.random((n_nodes, 2))

    complete = nx.Graph()
    complete.add_nodes_from(range(n_nodes))
    for u, v in itertools.combinations(range(n_nodes), 2):
        dist = float(np.linalg.norm(points[u] - points[v]))
        complete.add_edge(u, v, weight=dist)

    mst = nx.minimum_spanning_tree(complete, weight="weight")
    tree_edges = {tuple(sorted(e)) for e in mst.edges()}

    tree_dist = dict(nx.all_pairs_shortest_path_length(mst))
    non_tree_by_tree_path_len = sorted(
        (e for e in complete.edges() if tuple(sorted(e)) not in tree_edges),
        key=lambda e: tree_dist[e[0]][e[1]],
        reverse=True,
    )
    tie_edges = [tuple(sorted(e)) for e in non_tree_by_tree_path_len[:k_ties]]

    all_edges = sorted(tree_edges | set(tie_edges))
    return FeederGraph(n_nodes=n_nodes, k_ties=len(tie_edges), edges=tuple(all_edges))


def enumerate_spanning_trees(graph: FeederGraph) -> List[int]:
    """Returns every spanning tree of `graph` as a bitmask over edge
    indices (bit i set iff edges[i] is in the tree).

    Brute force over all size-(n_nodes-1) edge subsets, which is exactly
    C(n_edges, k_ties) since n_edges = n_nodes-1+k_ties -- polynomial in
    n_nodes for the fixed, small k_ties this study uses, which is the
    whole point of restricting to sparse graphs (see methodology.md).
    """
    n = graph.n_nodes
    m = graph.n_edges
    target_size = n - 1
    edges = graph.edges
    trees = []
    for combo in itertools.combinations(range(m), target_size):
        # Lightweight union-find tree check -- avoids constructing a
        # networkx Graph per candidate, which dominates runtime once
        # C(n_edges, k_ties) reaches the hundreds of thousands (e.g. real
        # feeder topologies with k_ties=5). Cross-checked against the
        # original nx.is_tree-based version on every synthetic instance
        # this repo tests (see tests via verify_correctness.py) before
        # relying on it here.
        parent = list(range(n))

        def find(x: int) -> int:
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        is_tree = True
        for i in combo:
            u, v = edges[i]
            ru, rv = find(u), find(v)
            if ru == rv:
                is_tree = False
                break
            parent[ru] = rv
        if is_tree:
            mask = 0
            for i in combo:
                mask |= 1 << i
            trees.append(mask)
    return trees
