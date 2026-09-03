"""Diagnostic (not part of the measurement pipeline): investigates *why*
witness search fails to stay small on the real IEEE 33-bus instance
(573/630 candidates dropped -- see results/real_feeder_results.csv),
unlike every synthetic instance in results/scaling_results.csv (0 dropped
candidates everywhere).

Hypothesis: `find_witness_set`'s search for a small witness set is really
searching for a short description of "which side of the cut induced by
removing tree edge e does edge f fall on". For a fixed spanning tree T,
T-e+f is valid iff f lies on the fundamental cycle formed by adding f to
T (equivalently, iff e lies on the tree-path between f's endpoints in T).
The number of independent cycles is bounded by k_ties (5 for IEEE 33-bus,
3 for the synthetic sweep) -- but a cycle's LENGTH (how many tree edges
it passes through) is not bounded by k_ties at all. `graphs.py`'s
synthetic model picks tie edges as the geometrically NEAREST non-tree
pairs, which plausibly keeps fundamental cycles short (few tree edges
between geometrically close endpoints). Real published tie switches exist
specifically to link distant parts of a feeder for reconfiguration
redundancy -- the opposite bias. If real tie edges have much longer
fundamental cycles, that directly predicts the witness-search failures
(a cycle spanning L tree edges means, in the worst case, a witness needs
to distinguish among those L edges' configurations).

This script measures fundamental-cycle length (tree-path length between
a tie edge's endpoints, in edges) for both the real IEEE 33-bus tree and
several synthetic instances at a comparable scale, and reports the
distributions side by side.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import networkx as nx
import pandapower.networks as pn

from graphs import generate_feeder_graph


def real_ieee33_cycle_lengths():
    net = pn.case33bw()
    n_nodes = len(net.bus)
    tree_edges = []
    tie_edges = []
    for _, row in net.line.iterrows():
        e = tuple(sorted((int(row["from_bus"]), int(row["to_bus"]))))
        (tree_edges if row["in_service"] else tie_edges).append(e)

    g = nx.Graph()
    g.add_nodes_from(range(n_nodes))
    g.add_edges_from(tree_edges)
    assert nx.is_tree(g), "in-service lines must form a spanning tree"

    lengths = []
    for u, v in tie_edges:
        path = nx.shortest_path(g, u, v)
        lengths.append(len(path) - 1)  # edges on the path = fundamental cycle length - 1 (plus the tie edge itself)
    return tie_edges, lengths, n_nodes


def synthetic_cycle_lengths(n_nodes: int, k_ties: int, seed: int):
    fg = generate_feeder_graph(n_nodes, k_ties, seed)
    all_edges = set(fg.edges)

    complete = nx.Graph()
    complete.add_nodes_from(range(n_nodes))
    import numpy as np
    import itertools

    rng = np.random.default_rng(seed)
    points = rng.random((n_nodes, 2))
    for u, v in itertools.combinations(range(n_nodes), 2):
        complete.add_edge(u, v, weight=float(np.linalg.norm(points[u] - points[v])))
    mst = nx.minimum_spanning_tree(complete, weight="weight")
    tree_edges = {tuple(sorted(e)) for e in mst.edges()}
    tie_edges = [e for e in fg.edges if e not in tree_edges]

    g = nx.Graph()
    g.add_nodes_from(range(n_nodes))
    g.add_edges_from(tree_edges)
    assert nx.is_tree(g)

    lengths = []
    for u, v in tie_edges:
        path = nx.shortest_path(g, u, v)
        lengths.append(len(path) - 1)
    return tie_edges, lengths


def main() -> None:
    tie_edges, lengths, n_nodes = real_ieee33_cycle_lengths()
    print(f"IEEE 33-bus (real, n_nodes={n_nodes}, k_ties={len(tie_edges)}):")
    for (u, v), L in zip(tie_edges, lengths):
        print(f"  tie edge ({u},{v}): fundamental cycle spans {L} tree edges")
    print(f"  max={max(lengths)}  mean={sum(lengths)/len(lengths):.1f}  min={min(lengths)}")

    print()
    print("Synthetic (Euclidean-MST + nearest ties), n_nodes=33, k_ties=5, seeds 0-4:")
    all_synth_lengths = []
    for seed in range(5):
        tie_edges_s, lengths_s = synthetic_cycle_lengths(33, 5, seed)
        all_synth_lengths.extend(lengths_s)
        print(f"  seed={seed}: lengths={lengths_s}  max={max(lengths_s)}")
    print(
        f"  overall: max={max(all_synth_lengths)}  "
        f"mean={sum(all_synth_lengths)/len(all_synth_lengths):.1f}  min={min(all_synth_lengths)}"
    )


if __name__ == "__main__":
    main()
