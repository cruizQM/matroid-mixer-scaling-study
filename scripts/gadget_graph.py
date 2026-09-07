"""The FULL gadget graph (root + buckets + items + leaves), for the
classical (CP-SAT) hardness check only -- NOT used by the QAOA circuit
(partition_mixer.py), which only needs the smaller one-hot qubit register
(see partition_gadget.py's module docstring for why the rest is forced
and can be dropped)."""

from __future__ import annotations

from dataclasses import dataclass

from partition_gadget import PartitionGadget


@dataclass(frozen=True)
class GadgetGraph:
    n_nodes: int
    edges: tuple[tuple[int, int], ...]
    root: int
    bipartite_edge_of: dict[int, tuple[int, int]]  # edge index -> (u_node, v_node), free layer only


def build_gadget_graph(gadget: PartitionGadget) -> GadgetGraph:
    m, k = gadget.m, gadget.k
    root = 0
    u_nodes = list(range(1, 1 + m))
    v_nodes = list(range(1 + m, 1 + m + k))
    next_node = 1 + m + k

    edges: list[tuple[int, int]] = [(root, u) for u in u_nodes]
    bipartite_edge_of: dict[int, tuple[int, int]] = {}
    for v in v_nodes:
        for u in u_nodes:
            bipartite_edge_of[len(edges)] = (u, v)
            edges.append((u, v))

    for vi, v in enumerate(v_nodes):
        for _ in range(gadget.items[vi] - 1):
            leaf = next_node
            next_node += 1
            edges.append((v, leaf))

    canon_edges = tuple((min(u, v), max(u, v)) for u, v in edges)
    return GadgetGraph(n_nodes=next_node, edges=canon_edges, root=root, bipartite_edge_of=bipartite_edge_of)
