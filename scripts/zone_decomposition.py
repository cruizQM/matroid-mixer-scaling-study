"""Tests a zone-decomposition approach for the matroid basis-exchange
mixer: partition the graph into zones via a min tie-line-cut criterion,
solve a matroid-basis mixer independently per zone (small qubit count
each), plus one small "assembly" mixer on the contracted graph of
boundary/tie-line edges between zones. A standard graphic-matroid
contraction/deletion property guarantees the union of per-zone spanning
trees plus the assembly spanning tree is a spanning tree of the whole
graph -- no reformulation of the mixer construction is needed at either
scale, only the edge set changes.

This exists to test a specific hypothesis raised after the real IEEE
33-bus validation found the full-graph witness-conditioned mixer
impractical (see docs/circuit-validity.md's real-feeder section): if
zones are chosen so that each zone's *internal* tie edges have short
fundamental cycles (unlike the whole graph's long-range ties), does each
zone's witness search stay small (target: < 10), even though the whole
graph's did not?

**Honest caveat**: the IEEE 33-bus test feeder is a single-substation
network (`pandapower`'s bus data has no zone/substation field -- checked
directly, all buses report zone=1.0), unlike a multi-substation
service-area context, where "zones" would naturally align with substation
boundaries. "Zones" here are therefore a graph min-cut partition of one
feeder into sub-regions, standing in for a substation-based partition --
a legitimate test of whether the *technique* (bound witness size by
shrinking each subproblem) works, but not a literal instance of "one zone
per substation" the way a multi-feeder network would provide.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, List, Tuple

sys.path.insert(0, str(Path(__file__).resolve().parent))

import networkx as nx

from graphs import FeederGraph, enumerate_spanning_trees
from mixer import build_matroid_mixer
from real_feeders import load_ieee33


def _split_connected(g: nx.Graph, nodes: set) -> Tuple[set, set]:
    """Splits `nodes` (which must induce a connected subgraph of `g`) into
    two node sets that ALSO each induce connected subgraphs, guaranteed:
    take any spanning tree of the induced subgraph, then remove one tree
    edge. Removing a tree edge always yields exactly two tree components,
    and each is connected by construction (it contains its own spanning-
    tree edges) -- unlike a generic min-cut bisection (e.g. Kernighan-Lin),
    which only minimizes cross edges and can accidentally produce a piece
    whose *induced* subgraph (using the full edge set, not just the tree)
    is disconnected if that piece happens to own none of the graph's
    "extra" edges. This matters because the decomposition property this
    script relies on requires each zone to have its own spanning tree to
    hand off to the per-zone mixer -- an internally-disconnected zone
    breaks that precondition outright (confirmed the hard way: an earlier
    version of this script using plain Kernighan-Lin bisection produced
    exactly this failure on 3 of 4 tested configurations).

    Which tree edge to remove is chosen primarily to keep the two pieces
    balanced in size (matching METIS's own balance constraint -- an
    unbalanced split defeats the point of shrinking each zone's qubit
    count, even if it happens to minimize the cut), tie-broken by
    minimizing the induced subgraph's OTHER (non-tree) edges that end up
    crossing the cut -- an approximation of a "min tie-line-cut
    criterion" (real METIS is not a dependency of this repo). An
    earlier version of this function prioritized cut size first, which
    peeled off small leaf stubs instead of bisecting (one 28-of-33-node
    zone out of 3) -- caught by inspecting the resulting zone sizes, not
    assumed balanced.
    """
    sub = g.subgraph(nodes)
    t0 = nx.minimum_spanning_tree(sub)
    best = None
    for u, v in t0.edges():
        t0_minus = t0.copy()
        t0_minus.remove_edge(u, v)
        comps = list(nx.connected_components(t0_minus))
        a, b = comps[0], comps[1]
        balance = abs(len(a) - len(b))
        crossing = sum(1 for x, y in sub.edges() if (x in a) != (y in a))
        score = (balance, crossing)
        if best is None or score < best[0]:
            best = (score, a, b)
    return best[1], best[2]


def partition_zones(graph: FeederGraph, target_zones: int) -> Dict[int, int]:
    """Recursively splits the graph into `target_zones` zones, each
    GUARANTEED to induce a connected subgraph (see `_split_connected`).
    Returns {node: zone_id}."""
    g = graph.networkx_graph()
    pieces: List[set] = [set(g.nodes())]
    while len(pieces) < target_zones:
        pieces.sort(key=len, reverse=True)
        biggest = pieces.pop(0)
        if len(biggest) < 2:
            pieces.append(biggest)
            break
        a, b = _split_connected(g, biggest)
        pieces.append(a)
        pieces.append(b)
    zone_of: Dict[int, int] = {}
    for zid, piece in enumerate(pieces):
        for node in piece:
            zone_of[node] = zid
    return zone_of


def partition_zones_by_size(graph: FeederGraph, target_zone_size: int) -> Dict[int, int]:
    """Like `partition_zones`, but takes a target *zone size* (node count)
    instead of a target zone *count* -- the number of zones scales with
    `graph.n_nodes / target_zone_size`, which is what a scaling sweep
    over network size needs (a fixed zone count would make zones grow
    with the network; a fixed zone size keeps them bounded)."""
    target_zones = max(1, round(graph.n_nodes / target_zone_size))
    return partition_zones(graph, target_zones)


def build_zone_subgraph(graph: FeederGraph, zone_of: Dict[int, int], zone_id: int) -> FeederGraph:
    nodes = sorted(n for n, z in zone_of.items() if z == zone_id)
    remap = {n: i for i, n in enumerate(nodes)}
    edges = tuple(
        sorted((remap[u], remap[v]) for (u, v) in graph.edges if zone_of[u] == zone_id and zone_of[v] == zone_id)
    )
    return FeederGraph(n_nodes=len(nodes), k_ties=0, edges=edges)


def build_assembly_graph(graph: FeederGraph, zone_of: Dict[int, int], n_zones: int) -> Tuple[FeederGraph, List[Tuple[int, int]]]:
    """Contracted graph: one supernode per zone, one edge per
    boundary/tie-line candidate (edges with endpoints in different
    zones). Kept as a multigraph (parallel supernode edges allowed) --
    `enumerate_spanning_trees`/`build_matroid_mixer` only ever reason
    about edge *indices*, never deduplicate by endpoint pair, so this is
    safe without special-casing."""
    boundary_edges = [(u, v) for (u, v) in graph.edges if zone_of[u] != zone_of[v]]
    edges = tuple(sorted((zone_of[u], zone_of[v]) if zone_of[u] < zone_of[v] else (zone_of[v], zone_of[u]) for u, v in boundary_edges))
    return FeederGraph(n_nodes=n_zones, k_ties=0, edges=edges), boundary_edges


def main() -> None:
    graph = load_ieee33()
    print(f"Full graph: n_nodes={graph.n_nodes}  n_qubits={graph.n_edges}")

    for target_zones in (3, 4, 5, 6):
        print(f"\n=== target_zones={target_zones} ===")
        zone_of = partition_zones(graph, target_zones)
        n_zones = len(set(zone_of.values()))
        zone_sizes = {z: sum(1 for v in zone_of.values() if v == z) for z in set(zone_of.values())}
        print(f"actual zones formed: {n_zones}, sizes={zone_sizes}")

        max_witness_overall = 0
        any_disconnected = False
        for zid in sorted(set(zone_of.values())):
            zg = build_zone_subgraph(graph, zone_of, zid)
            if zg.n_edges == 0 or zg.n_nodes < 2:
                print(f"  zone {zid}: n_nodes={zg.n_nodes} n_qubits={zg.n_edges} (trivial, skipped)")
                continue
            trees = enumerate_spanning_trees(zg)
            if not trees:
                print(f"  zone {zid}: n_nodes={zg.n_nodes} n_qubits={zg.n_edges} -- NOT INTERNALLY CONNECTED (no spanning tree; expected, some/all connectivity comes from boundary edges)")
                any_disconnected = True
                continue
            construction = build_matroid_mixer(zg, trees)
            max_w = max((t.control_count for t in construction.terms), default=0)
            max_witness_overall = max(max_witness_overall, max_w)
            print(
                f"  zone {zid}: n_nodes={zg.n_nodes} n_qubits={zg.n_edges} |B|={len(trees)} "
                f"terms={len(construction.terms)} max_witness={max_w} dropped={construction.dropped_candidates} "
                f"connected={construction.fully_connected}"
            )

        assembly_graph, boundary_edges = build_assembly_graph(graph, zone_of, n_zones)
        print(f"assembly graph: n_zones={n_zones} n_boundary_edges={len(boundary_edges)}  boundary_edges(orig)={boundary_edges}")
        if assembly_graph.n_edges > 0 and assembly_graph.n_nodes >= 2:
            a_trees = enumerate_spanning_trees(assembly_graph)
            if a_trees:
                a_construction = build_matroid_mixer(assembly_graph, a_trees)
                a_max_w = max((t.control_count for t in a_construction.terms), default=0)
                max_witness_overall = max(max_witness_overall, a_max_w)
                print(
                    f"assembly: |B|={len(a_trees)} terms={len(a_construction.terms)} max_witness={a_max_w} "
                    f"dropped={a_construction.dropped_candidates} connected={a_construction.fully_connected}"
                )
            else:
                print("assembly: NOT CONNECTED (no spanning tree over zone supernodes) -- partition choice fails at assembly level")

        print(f"--> max witness size across zones + assembly: {max_witness_overall}  (target: < 10)")


if __name__ == "__main__":
    main()
