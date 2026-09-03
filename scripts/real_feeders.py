"""Loads real, published distribution-feeder topologies to validate the
synthetic-graph scaling study (`results/scaling_results.csv`) against
actual data, rather than only the geometric-MST proxy in `graphs.py`.

## IEEE 33-bus (Baran & Wu, 1989)

Loaded via `pandapower.networks.case33bw` -- the standard IEEE distribution
test feeder, and the network `README.md`'s data source #1 (challenge
statement, IEEE Test Networks link) itself points to. **Not** assumed
correct from the citation: `load_ieee33` checks programmatically that the
33 in-service lines form an actual spanning tree (32 edges, no cycles, one
component) before using them, and raises if that assumption is violated.
The 5 out-of-service lines are the network's own published tie switches --
exactly this study's "k_ties" structure, with real data replacing the
synthetic generator's Euclidean-MST-plus-nearest-ties construction.

## CATS (California Test System) -- investigated, found not applicable

The challenge statement's data source #2. Cloned and inspected directly
(not assumed from the README alone): `CaliforniaTestSystem.m` (MATPOWER
format) has on the order of 8,800 buses and 10,800 branches -- MORE
branches than buses, i.e. the network contains cycles/loops, confirmed by
the repository's own description as built from "California's actual
transmission lines, substations, and power plants." This is a
transmission-scale, meshed network, not a distribution feeder. Real
transmission networks are deliberately operated meshed for reliability,
not radially -- the graphic-matroid spanning-tree constraint this study is
about does not describe how CATS is actually operated, so running this
study's mixer construction on it would not be a meaningful match. Not
used, rather than force-fit for the sake of having a second real dataset.
"""

from __future__ import annotations

import networkx as nx
import pandapower.networks as pn

from graphs import FeederGraph


def load_ieee33() -> FeederGraph:
    net = pn.case33bw()
    n_nodes = len(net.bus)

    tree_edges = []
    tie_edges = []
    for _, row in net.line.iterrows():
        e = tuple(sorted((int(row["from_bus"]), int(row["to_bus"]))))
        (tree_edges if row["in_service"] else tie_edges).append(e)

    g_tree = nx.Graph()
    g_tree.add_nodes_from(range(n_nodes))
    g_tree.add_edges_from(tree_edges)
    if not nx.is_tree(g_tree):
        raise ValueError(
            "case33bw's in-service lines do not form a spanning tree -- "
            "the assumption this loader relies on is violated for this "
            "pandapower version's bundled data; do not proceed silently."
        )

    all_edges = tuple(sorted(set(tree_edges) | set(tie_edges)))
    return FeederGraph(n_nodes=n_nodes, k_ties=len(tie_edges), edges=all_edges)
