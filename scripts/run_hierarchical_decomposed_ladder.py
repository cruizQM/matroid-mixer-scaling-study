"""Tests whether HIERARCHICAL decomposition -- recursing zone_decomposition
on the assembly graph itself, when the assembly graph is still large,
instead of measuring it directly -- reduces total circuit cost further
than the flat (one-level) decomposition in run_decomposed_cost_aware_ladder.py.

## Why this might matter, found by direct probing first

Naively shrinking target_zone_size does NOT help: at long_log, n_nodes=150,
total cost went from 856 CX (target_zone_size=8) to 5,836 CX
(target_zone_size=2) as zones shrank, because the ASSEMBLY graph (one
supernode per zone, one edge per cross-zone tie) gets DENSER as zones
shrink -- more, smaller zones means more edges end up crossing zone
boundaries. Shrinking zones just moves cost into an increasingly complex
assembly graph instead of removing it.

Hierarchical decomposition addresses this directly: instead of shrinking
the ORIGINAL zones, keep target_zone_size at its sweet spot (8, matching
this repo's existing convention) at the TOP level, and, if the ASSEMBLY
graph itself is still bigger than a threshold, decompose IT the same way
(zones + a smaller assembly graph), recursively, until every subproblem
at every level is small.

A first version of the recursive step used a FIXED schedule (halve
`target_zone_size` at each level) -- this fixed the worst instability
(one seed's assembly graph concentrating almost all its density into a
single sub-zone, costing 16,120 CX alone) but was still somewhat
arbitrary. `_choose_target_zone_size` replaces the fixed schedule with a
density-aware one: scale `target_zone_size` DOWN in proportion to how
much denser the CURRENT graph is than the original (measured directly at
each level, `graph.n_edges / graph.n_nodes`), since the assembly graph is
systematically denser than the graph it came from (collapsing every
cross-zone tie into a boundary edge on far fewer supernodes). Verified to
do measurably better than the fixed-halving schedule on the same
instances, not just to be more principled in the abstract (long_linear,
n_nodes=150: mean cost 4,308.7 -> 3,767.3 CX across the same 3 seeds,
worst-case seed also improved).

Writes results/hierarchical_decomposed_ladder_results.csv and _summary.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from graphs import FeederGraph
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones_by_size
from run_decomposed_cost_aware_ladder import measure_subproblem
from run_cost_aware_scaling_ladder import CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "hierarchical_decomposed_ladder_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "hierarchical_decomposed_ladder_summary.csv"

BASE_TARGET_ZONE_SIZE = 8
RECURSE_THRESHOLD = 12  # if the assembly graph has more nodes than this, decompose it too instead of measuring directly
MAX_DEPTH = 6
MIN_TARGET_ZONE_SIZE = 2


def _density(graph: FeederGraph) -> float:
    return graph.n_edges / graph.n_nodes


def _choose_target_zone_size(graph: FeederGraph, base_density: float) -> int:
    """Scales `target_zone_size` INVERSELY to how dense `graph` is
    relative to the ORIGINAL graph's density (`base_density`, measured
    once at depth 0) -- a fixed halving-per-level schedule (tried first)
    already fixed the worst instability, but this does measurably better
    (mean cost 4,308.7 -> 3,767.3 on the same long_linear/n_nodes=150
    seeds, and the worst-case seed improves too) because it responds to
    the ACTUAL measured density at each level instead of an arbitrary
    per-level schedule. Mechanism: the assembly graph is systematically
    denser than the graph it came from (it collapses every cross-zone
    tie into one boundary edge on far fewer supernodes) -- a graph twice
    as dense needs roughly half the target zone size to keep each
    resulting sub-zone's OWN tie density comparable to what worked at
    the top level, rather than concentrating it into fewer, harder
    zones."""
    d = _density(graph)
    scaled = round(BASE_TARGET_ZONE_SIZE * base_density / d)
    return max(MIN_TARGET_ZONE_SIZE, min(BASE_TARGET_ZONE_SIZE, scaled))


def decompose_total(graph, seed: int, base_density: float = None, depth: int = 0) -> dict:
    if base_density is None:
        base_density = _density(graph)
    target_zone_size = _choose_target_zone_size(graph, base_density)
    zone_of = partition_zones_by_size(graph, target_zone_size)
    n_zones = len(set(zone_of.values()))
    total_cx = 0
    total_depth_metric = 0
    total_terms = 0
    total_active = 0
    max_levels = depth + 1

    for zid in sorted(set(zone_of.values())):
        zg = build_zone_subgraph(graph, zone_of, zid)
        r = measure_subproblem(zg, seed)
        total_cx += r["cx"]
        total_depth_metric += r["depth"]
        total_terms += r["n_terms"]
        total_active += r["active_terms"]

    assembly_graph, _ = build_assembly_graph(graph, zone_of, n_zones)
    if assembly_graph.n_nodes > RECURSE_THRESHOLD and depth < MAX_DEPTH:
        sub = decompose_total(assembly_graph, seed, base_density, depth + 1)
        total_cx += sub["cx"]
        total_depth_metric += sub["depth"]
        total_terms += sub["n_terms"]
        total_active += sub["active_terms"]
        max_levels = max(max_levels, sub["max_levels"])
    else:
        r = measure_subproblem(assembly_graph, seed)
        total_cx += r["cx"]
        total_depth_metric += r["depth"]
        total_terms += r["n_terms"]
        total_active += r["active_terms"]

    return dict(cx=total_cx, depth=total_depth_metric, n_terms=total_terms, active_terms=total_active, max_levels=max_levels)


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)
    result = decompose_total(graph, seed)
    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits_whole_graph": graph.n_edges,
        "max_hierarchy_levels": result["max_levels"],
        "total_terms": result["n_terms"], "total_active_terms": result["active_terms"],
        "total_cx": result["cx"], "total_depth": result["depth"],
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in ("long_log", "long_linear"):
        for n_nodes in N_NODES_RANGE:
            for seed in range(SEEDS_PER_SIZE):
                row = run_one(condition, n_nodes, seed)
                rows.append(row)
                print(f"  {row}", flush=True)

    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    summary_rows = []
    for condition in ("long_log", "long_linear"):
        for n_nodes in N_NODES_RANGE:
            group = [r for r in rows if r["condition"] == condition and r["n_nodes"] == n_nodes]
            cx_vals = [r["total_cx"] for r in group]
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes,
                "max_levels_mean": round(statistics.mean(r["max_hierarchy_levels"] for r in group), 1),
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                "cx_min": min(cx_vals), "cx_max": max(cx_vals),
            })
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")
    for row in summary_rows:
        print(f"  {row}")


if __name__ == "__main__":
    main()
