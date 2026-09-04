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


def load_cigre_mv() -> FeederGraph:
    """CIGRE medium-voltage benchmark distribution network (Task Force
    C6.04.02) -- loaded via `pandapower.networks.create_cigre_network_mv`.
    A standardized, widely-used realistic benchmark topology (not an
    as-built utility feeder like case33bw, but not a random synthetic
    graph either) -- a second, independently-sourced, differently-sized
    (15 buses vs. case33bw's 33) real-topology check.

    Same discipline as `load_ieee33`: ties are the lines with an open
    switch or `in_service=False` (checked directly, not assumed), plus
    this network's 2 transformers are added as always-tree edges (a
    transformer connecting two voltage levels is not a switchable tie,
    and excluding it here would leave buses on the other side spuriously
    disconnected) -- and the resulting tree-edge set is checked to
    actually form a spanning tree before use, exactly like `load_ieee33`
    raises rather than proceeding silently if that check fails."""
    net = pn.create_cigre_network_mv()
    n_nodes = len(net.bus)

    open_line_switch = set()
    if len(net.switch):
        for _, row in net.switch[net.switch.et == "l"].iterrows():
            if not row["closed"]:
                open_line_switch.add(row["element"])

    tree_edges = []
    tie_edges = []
    for lidx, line in net.line.iterrows():
        e = tuple(sorted((int(line["from_bus"]), int(line["to_bus"]))))
        is_tie = (not line.get("in_service", True)) or (lidx in open_line_switch)
        (tie_edges if is_tie else tree_edges).append(e)
    for _, t in net.trafo.iterrows():
        tree_edges.append(tuple(sorted((int(t.hv_bus), int(t.lv_bus)))))

    g_tree = nx.Graph()
    g_tree.add_nodes_from(range(n_nodes))
    g_tree.add_edges_from(tree_edges)
    if not nx.is_tree(g_tree):
        raise ValueError(
            "CIGRE MV's in-service lines + transformers do not form a spanning tree -- "
            "the assumption this loader relies on is violated for this "
            "pandapower version's bundled data; do not proceed silently."
        )

    all_edges = tuple(sorted(set(tree_edges) | set(tie_edges)))
    return FeederGraph(n_nodes=n_nodes, k_ties=len(tie_edges), edges=all_edges)
