"""Generates illustrative diagrams for the README/docs -- what the
fundamental-cycle problem actually looks like on a graph, and what zone
decomposition does to it. These are explanatory pictures, not
measurements: no numbers here feed into `results/*.csv` or any claim in
`methodology.md`. Generated (not hand-drawn) so they stay reproducible
and consistent with the rest of this repo's "nothing hand-edited"
convention; regenerate with `python scripts/plot_illustrations.py`.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from graphs import generate_feeder_graph, generate_feeder_graph_long_range_ties
from zone_decomposition import build_assembly_graph, partition_zones_by_size

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"

# One running example graph, reused across every illustration in this
# script (and hence across the README/docs) for visual coherence: same
# node layout throughout, so a reader isn't asked to re-orient between
# figures. `SMALL_N_NODES` is small enough to label every node clearly;
# the decomposition figure needs more scale to be illustrative, so it
# uses a larger instance built from the same seed instead of the exact
# same graph.
SEED = 3
SMALL_N_NODES = 14


def _mst_and_positions(n_nodes: int, seed: int):
    """Reproduces the same MST + node layout `graphs.py`'s generators use
    internally (same seed -> same random points -> same MST), so this
    stays visually consistent with the actual generators rather than
    drawing something unrelated."""
    rng = np.random.default_rng(seed)
    points = rng.random((n_nodes, 2))
    complete = nx.Graph()
    complete.add_nodes_from(range(n_nodes))
    for u, v in itertools.combinations(range(n_nodes), 2):
        complete.add_edge(u, v, weight=float(np.linalg.norm(points[u] - points[v])))
    mst = nx.minimum_spanning_tree(complete, weight="weight")
    tree_edges = {tuple(sorted(e)) for e in mst.edges()}
    pos = {i: tuple(points[i]) for i in range(n_nodes)}
    return mst, tree_edges, pos


def plot_feasibility_example() -> None:
    """The running example graph (`SMALL_N_NODES` nodes, `SEED`, the
    synthetic/nearest-tie model with one candidate tie switch): what
    "feasible" actually means. Left panel closes exactly the backbone
    (a spanning tree -- radial, every node reached, no loop). Right panel
    additionally closes the one tie switch -- one edge too many for a
    14-node tree, which necessarily creates exactly one loop, highlighted.
    This is the constraint the rest of this repo is about, before any
    circuit or scaling detail."""
    graph = generate_feeder_graph(SMALL_N_NODES, 1, SEED)
    mst, tree_edges, pos = _mst_and_positions(SMALL_N_NODES, SEED)
    tie_edge = next(e for e in graph.edges if e not in tree_edges)
    cycle_path = nx.shortest_path(mst, tie_edge[0], tie_edge[1])
    cycle_path_edges = {tuple(sorted((cycle_path[i], cycle_path[i + 1]))) for i in range(len(cycle_path) - 1)}

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))

    ax = axes[0]
    nx.draw_networkx_edges(mst, pos, ax=ax, edge_color="#2E5FA3", width=2.5)
    nx.draw_networkx_edges(nx.Graph([tie_edge]), pos, ax=ax, edge_color="#B0B0B0", width=2.0, style="dashed")
    nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
    nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title("Feasible: tie switch open\n(radial -- every node reached, no loop)", fontsize=11)
    ax.axis("off")

    ax = axes[1]
    nx.draw_networkx_edges(mst, pos, ax=ax, edgelist=[e for e in tree_edges if e not in cycle_path_edges], edge_color="#2E5FA3", width=2.5)
    nx.draw_networkx_edges(mst, pos, ax=ax, edgelist=list(cycle_path_edges), edge_color="#C0392B", width=3.0)
    nx.draw_networkx_edges(nx.Graph([tie_edge]), pos, ax=ax, edge_color="#C0392B", width=3.0)
    nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
    nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title("Infeasible: tie switch also closed\n(one loop -- highlighted in red)", fontsize=11)
    ax.axis("off")

    fig.suptitle(
        "The problem: a feeder's switches must stay a spanning tree of the network graph\n"
        "-- qubit i = 1 iff switch i is closed; feasible states are exactly the spanning trees",
        fontsize=12,
    )
    handles = [
        plt.Line2D([0], [0], color="#2E5FA3", lw=2.5, label="closed switch (in service)"),
        plt.Line2D([0], [0], color="#B0B0B0", lw=2.0, ls="dashed", label="open switch (tie, unused)"),
        plt.Line2D([0], [0], color="#C0392B", lw=3, label="closed switch, on the loop"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = RESULTS_DIR / "illustration_feeder_problem.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


def plot_basis_exchange_move() -> None:
    """The running example graph again: what one mixer move actually
    does. Closing the tie switch alone would create a loop (the
    infeasible panel in `plot_feasibility_example`) -- a basis exchange
    instead closes the tie switch AND opens one tree edge on the same
    loop at once, landing on a different, still-feasible spanning tree.
    This single operation, repeated and conditioned correctly, is the
    entire mixer; the rest of this repo measures how expensive
    conditioning it correctly turns out to be."""
    graph = generate_feeder_graph(SMALL_N_NODES, 1, SEED)
    mst, tree_edges, pos = _mst_and_positions(SMALL_N_NODES, SEED)
    tie_edge = next(e for e in graph.edges if e not in tree_edges)
    cycle_path = nx.shortest_path(mst, tie_edge[0], tie_edge[1])
    cycle_path_edges = [tuple(sorted((cycle_path[i], cycle_path[i + 1]))) for i in range(len(cycle_path) - 1)]
    edge_out = cycle_path_edges[0]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5.5))

    ax = axes[0]
    nx.draw_networkx_edges(mst, pos, ax=ax, edge_color="#2E5FA3", width=2.5)
    nx.draw_networkx_edges(nx.Graph([tie_edge]), pos, ax=ax, edge_color="#B0B0B0", width=2.0, style="dashed")
    nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
    nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title("Before: one feasible tree", fontsize=11)
    ax.axis("off")

    ax = axes[1]
    nx.draw_networkx_edges(mst, pos, ax=ax, edgelist=[e for e in tree_edges if e != edge_out], edge_color="#2E5FA3", width=2.5)
    nx.draw_networkx_edges(mst, pos, ax=ax, edgelist=[edge_out], edge_color="#C0392B", width=2.5, style="dashed")
    nx.draw_networkx_edges(nx.Graph([tie_edge]), pos, ax=ax, edge_color="#2E8B57", width=3.0)
    nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
    nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title("The move: open one, close the tie\n(one basis exchange)", fontsize=11)
    ax.axis("off")

    ax = axes[2]
    new_tree_edges = [e for e in tree_edges if e != edge_out] + [tie_edge]
    nx.draw_networkx_edges(nx.Graph(new_tree_edges), pos, ax=ax, edge_color="#2E5FA3", width=2.5)
    nx.draw_networkx_edges(nx.Graph([edge_out]), pos, ax=ax, edge_color="#B0B0B0", width=2.0, style="dashed")
    nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
    nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
    ax.set_title("After: a different feasible tree", fontsize=11)
    ax.axis("off")

    fig.suptitle(
        "What we're measuring: the mixer is built entirely from moves like this one\n"
        "-- one closed switch for one open switch, never landing outside the feasible set",
        fontsize=12,
    )
    handles = [
        plt.Line2D([0], [0], color="#2E5FA3", lw=2.5, label="closed switch"),
        plt.Line2D([0], [0], color="#B0B0B0", lw=2.0, ls="dashed", label="open switch"),
        plt.Line2D([0], [0], color="#2E8B57", lw=3, label="switch closing"),
        plt.Line2D([0], [0], color="#C0392B", lw=2.5, ls="dashed", label="switch opening"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    out = RESULTS_DIR / "illustration_basis_exchange_move.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


def plot_fundamental_cycle_comparison() -> None:
    """Same running example graph as `plot_feasibility_example` (same
    node positions, same seed) on the left; the same backbone with a
    long-range tie instead of a nearest-neighbor one on the right. Shows
    the tie edge's fundamental cycle -- the tree path it closes --
    highlighted, so the only visible difference is which tie edge got
    picked and how far its cycle reaches."""
    n_nodes, seed = SMALL_N_NODES, SEED
    short_graph = generate_feeder_graph(n_nodes, 1, seed)
    long_graph = generate_feeder_graph_long_range_ties(n_nodes, 1, seed)
    mst, tree_edges, pos = _mst_and_positions(n_nodes, seed)

    fig, axes = plt.subplots(1, 2, figsize=(13, 6))
    for ax, graph, title in (
        (axes[0], short_graph, "synthetic model: nearest-neighbor tie"),
        (axes[1], long_graph, "real-feeder-like model: long-range tie"),
    ):
        tie_edge = next(e for e in graph.edges if e not in tree_edges)
        cycle_path = nx.shortest_path(mst, tie_edge[0], tie_edge[1])
        cycle_path_edges = {tuple(sorted((cycle_path[i], cycle_path[i + 1]))) for i in range(len(cycle_path) - 1)}

        nx.draw_networkx_edges(mst, pos, ax=ax, edge_color="#B0B0B0", width=1.5)
        nx.draw_networkx_edges(
            mst, pos, ax=ax, edgelist=list(cycle_path_edges), edge_color="#D9822B", width=3.0
        )
        nx.draw_networkx_edges(
            nx.Graph([tie_edge]), pos, ax=ax, edge_color="#C0392B", width=3.0, style="dashed"
        )
        nx.draw_networkx_nodes(mst, pos, ax=ax, node_color="#4472C4", node_size=260, edgecolors="white")
        nx.draw_networkx_labels(mst, pos, ax=ax, font_size=8, font_color="white")
        ax.set_title(f"{title}\nfundamental cycle length = {len(cycle_path_edges)} tree edges", fontsize=11)
        ax.axis("off")

    fig.suptitle(
        "The problem: a tie edge's fundamental cycle is exactly the qubits an exchange's\n"
        "validity can depend on -- short here, long there, same backbone",
        fontsize=12,
    )
    handles = [
        plt.Line2D([0], [0], color="#B0B0B0", lw=1.5, label="tree backbone (not on this cycle)"),
        plt.Line2D([0], [0], color="#D9822B", lw=3, label="tree edges on the fundamental cycle"),
        plt.Line2D([0], [0], color="#C0392B", lw=3, ls="dashed", label="the tie edge"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.02))
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    out = RESULTS_DIR / "illustration_fundamental_cycle.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


def plot_the_decomposition() -> None:
    """One graph, partitioned into zones (colored) with the boundary/tie
    edges that cross zone lines highlighted -- those boundary edges are
    exactly what the small assembly problem (inset) is built from. Same
    `SEED` as the other two illustrations (a larger instance from it --
    the running example graph is too small to make partitioning
    illustrative)."""
    n_nodes, seed, target_zone_size = 24, SEED, 8
    graph = generate_feeder_graph_long_range_ties(n_nodes, 5, seed)
    g = graph.networkx_graph()
    # A graph-layout algorithm (not the raw geometric coordinates used for
    # generation) draws this much more cleanly: zones are connectivity
    # clusters, not necessarily spatial ones, and Kamada-Kawai keeps
    # connected nodes close together, minimizing edge crossings.
    pos = nx.kamada_kawai_layout(g)

    zone_of = partition_zones_by_size(graph, target_zone_size)
    n_zones = len(set(zone_of.values()))
    _, boundary_edges = build_assembly_graph(graph, zone_of, n_zones)

    internal_edges = [e for e in graph.edges if e not in boundary_edges]

    cmap = plt.get_cmap("tab10")
    node_colors = [cmap(zone_of[n] % 10) for n in g.nodes()]

    fig, (ax_main, ax_assembly) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [2.2, 1]})

    nx.draw_networkx_edges(g, pos, ax=ax_main, edgelist=internal_edges, edge_color="#B0B0B0", width=1.5)
    nx.draw_networkx_edges(
        g, pos, ax=ax_main, edgelist=boundary_edges, edge_color="#C0392B", width=2.5, style="dashed"
    )
    nx.draw_networkx_nodes(g, pos, ax=ax_main, node_color=node_colors, node_size=240, edgecolors="white")
    nx.draw_networkx_labels(g, pos, ax=ax_main, font_size=7, font_color="white")
    ax_main.set_title(
        f"Full graph partitioned into {n_zones} zones\n"
        f"(one small matroid mixer solved per zone, independently)",
        fontsize=11,
    )
    ax_main.axis("off")
    handles = [
        plt.Line2D([0], [0], color="#B0B0B0", lw=1.5, label="within-zone edge"),
        plt.Line2D([0], [0], color="#C0392B", lw=2.5, ls="dashed", label="boundary/tie edge"),
    ]
    ax_main.legend(handles=handles, loc="lower center", fontsize=8, bbox_to_anchor=(0.5, -0.08), ncol=2)

    assembly_g = nx.Graph()
    assembly_g.add_nodes_from(range(n_zones))
    zone_centroid = {
        z: (
            sum(pos[n][0] for n in g.nodes() if zone_of[n] == z) / sum(1 for n in g.nodes() if zone_of[n] == z),
            sum(pos[n][1] for n in g.nodes() if zone_of[n] == z) / sum(1 for n in g.nodes() if zone_of[n] == z),
        )
        for z in range(n_zones)
    }
    for u, v in boundary_edges:
        assembly_g.add_edge(zone_of[u], zone_of[v])
    assembly_colors = [cmap(z % 10) for z in assembly_g.nodes()]
    nx.draw_networkx_edges(assembly_g, zone_centroid, ax=ax_assembly, edge_color="#C0392B", width=2.5)
    nx.draw_networkx_nodes(
        assembly_g, zone_centroid, ax=ax_assembly, node_color=assembly_colors, node_size=900, edgecolors="white"
    )
    nx.draw_networkx_labels(assembly_g, zone_centroid, ax=ax_assembly, font_size=10, font_color="white")
    ax_assembly.set_title(
        f"Assembly problem: {n_zones} zones contracted to\nsupernodes, {len(boundary_edges)} boundary edges",
        fontsize=11,
    )
    ax_assembly.axis("off")

    fig.suptitle("The decomposition: same construction, applied once per zone and once for the boundary", fontsize=12)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))
    out = RESULTS_DIR / "illustration_decomposition.png"
    fig.savefig(out, dpi=150, bbox_inches="tight")
    print(f"wrote {out}")


if __name__ == "__main__":
    RESULTS_DIR.mkdir(exist_ok=True)
    plot_feasibility_example()
    plot_basis_exchange_move()
    plot_fundamental_cycle_comparison()
    plot_the_decomposition()
