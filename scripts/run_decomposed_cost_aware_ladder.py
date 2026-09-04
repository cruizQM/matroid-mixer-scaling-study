"""Tests the fix identified as most promising after the "50-500 CX"
investigation: the cost-aware truncated mixer's floors (2,300-27,700 CX
across the escalating-realism ladder, `run_cost_aware_scaling_ladder.py`)
are driven by TERM COUNT (19-206 active terms), not by individual gate
width -- so a fix that reduces circuit cost per gate (e.g. ancilla-based
multi-controlled synthesis) wouldn't touch the actual bottleneck. Zone
decomposition (`zone_decomposition.py`, headline result 2's fix for the
EXACT construction) directly attacks term count instead: most of a
whole-graph construction's terms only exist to connect distant parts of
the network, and headline result 2 already measured this collapsing
witness size from 12-58 (whole-graph) to a flat ~3 (decomposed) -- this
tests whether the same collapse happens for TOTAL TERM COUNT and CX cost
when combined with the cost-aware truncated mixer, not just witness size
for the exact construction.

## Method

Same four conditions/sizes/seeds as the escalating-realism ladder. At
each instance: partition into zones (`zone_decomposition.partition_zones_by_size`,
`target_zone_size=8`, matching `run_decomposition_scaling_study.py`'s own
convention), then build a cost-aware truncated mixer
(`adaptive=False, cost_alpha=0.01`, the validated default -- see
`measure_subproblem`'s docstring for why this is `adaptive=False`, not
the `adaptive=True` an earlier version of this script used) on EACH zone
subgraph and the assembly graph independently, using exact enumeration
for each (zones are small enough that this stays fast, unlike the
whole-graph case). Reports TOTAL CX/depth summed across every
zone+assembly subproblem, and per-subproblem safety
(`leakage_trace.final_feasible_mass`, reported as the mean and worst
across all subproblems -- a per-subproblem check, not a full joint
end-to-end simulation of the combined mixer, the same scope
`run_decomposition_scaling_study.py` itself used for the exact
construction's leak-freeness check).

Writes results/decomposed_cost_aware_ladder_results.csv (per condition,
size, seed) and _summary.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from graphs import enumerate_spanning_trees
from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer
from math import comb
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones_by_size
from run_cost_aware_scaling_ladder import CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE, BETA

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "decomposed_cost_aware_ladder_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "decomposed_cost_aware_ladder_summary.csv"

MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 2
GAIN_PRICE = 0.01
TARGET_ZONE_SIZE = 8
WALK_STEPS = 300
N_STARTING_TREES_PER_SUBPROBLEM = 30
UNSAFE_TOL = 1e-6


MAX_EXACT_ENUMERATION_COMBINATIONS = 200_000


def measure_subproblem(graph, seed: int) -> dict:
    """Exact enumeration (`graphs.enumerate_spanning_trees`) whenever
    it's cheap -- most zones, being small by construction
    (`TARGET_ZONE_SIZE=8`), qualify, and exact is more faithful than a
    sample for a subproblem this small. Falls back to a walked-exchange-graph
    sample ONLY when exact enumeration's own cost
    (`C(n_edges, n_nodes-1)`) would be too large -- the ASSEMBLY graph is
    the reason this matters: one supernode per zone with every cross-zone
    tie collapsed into a boundary edge can end up genuinely dense even at
    small n_nodes (one real instance: 19 zone-supernodes, 31 boundary
    edges -- brute-force enumeration there alone took 81s, the same
    density trap `tie_density_sweep.py` already found, just reached via
    zone contraction instead of raw tie count). Using a sample
    UNCONDITIONALLY (tried first) gave a meaningfully different, WORSE
    result on the same instance (19,592 CX vs. 4,510 exact) -- the
    walked sample is a partial, connected-component-limited view, less
    informative for the majority-vote search than the true tree set when
    the true tree set is actually affordable to compute, so exact is
    preferred whenever it is."""
    if graph.n_nodes < 2:
        return dict(n_terms=0, active_terms=0, cx=0, depth=0, masses=[1.0])
    target_size = graph.n_nodes - 1
    if graph.n_edges >= target_size and comb(graph.n_edges, target_size) <= MAX_EXACT_ENUMERATION_COMBINATIONS:
        trees = enumerate_spanning_trees(graph)
    else:
        trees = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    if not trees:
        return dict(n_terms=0, active_terms=0, cx=0, depth=0, masses=[1.0])
    truncated = build_truncated_witness_mixer(
        graph, trees, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
        adaptive=False, cost_alpha=GAIN_PRICE,
        # NOT adaptive=True: this script originally used the adaptive search
        # (matching what looked, at the time, like the better default).
        # docs/scaling-ladder-and-decomposition.md later found adaptive does
        # NOT reliably beat fixed cost_alpha once properly re-validated --
        # confirmed to matter here too, not just at whole-graph scale: one
        # specific zone (of 19, long_linear/n_nodes=150/seed=0) cost 3,550 CX
        # under adaptive vs. 292 CX under fixed cost_alpha=0.01, a 12x gap,
        # and was the dominant contributor to that instance's total zone
        # cost. Fixed to use the validated default; every result in this
        # script (and run_best_of_both_ladder.py, run_hierarchical_decomposed_ladder.py,
        # which both build on it) was re-run after this fix, not left stale.
    )
    construction = truncated.construction
    active_terms = sum(1 for t in construction.terms if t.minimized_gate_count > 0)
    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx, depth = 0, 0

    rng = np.random.default_rng(seed)
    n_trees_to_check = min(N_STARTING_TREES_PER_SUBPROBLEM, len(trees))
    masses = []
    for _ in range(n_trees_to_check):
        start_tree = random_spanning_tree(graph, rng)
        masses.append(final_feasible_mass(graph, construction, start_tree, BETA, sparse=True))
    return dict(n_terms=len(construction.terms), active_terms=active_terms, cx=cx, depth=depth, masses=masses)


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)

    zone_of = partition_zones_by_size(graph, TARGET_ZONE_SIZE)
    n_zones = len(set(zone_of.values()))

    total = dict(n_terms=0, active_terms=0, cx=0, depth=0)
    all_masses = []
    for zid in sorted(set(zone_of.values())):
        zg = build_zone_subgraph(graph, zone_of, zid)
        result = measure_subproblem(zg, seed)
        for key in ("n_terms", "active_terms", "cx", "depth"):
            total[key] += result[key]
        all_masses.extend(result["masses"])

    assembly_graph, _ = build_assembly_graph(graph, zone_of, n_zones)
    result = measure_subproblem(assembly_graph, seed)
    for key in ("n_terms", "active_terms", "cx", "depth"):
        total[key] += result[key]
    all_masses.extend(result["masses"])

    unsafe = sum(1 for m in all_masses if m < 1.0 - UNSAFE_TOL)
    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits_whole_graph": graph.n_edges, "n_zones": n_zones,
        "total_n_terms": total["n_terms"], "total_active_terms": total["active_terms"],
        "total_cx": total["cx"], "total_depth": total["depth"],
        "n_masses_checked": len(all_masses),
        "unsafe_rate": round(unsafe / len(all_masses), 4) if all_masses else 0.0,
        "mean_feasible_mass": round(float(np.mean(all_masses)), 4) if all_masses else 1.0,
        "worst_feasible_mass": round(float(min(all_masses)), 4) if all_masses else 1.0,
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in CONDITIONS:
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
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            group = [r for r in rows if r["condition"] == condition and r["n_nodes"] == n_nodes]
            cx_vals = [r["total_cx"] for r in group]
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes, "k_ties": group[0]["k_ties"],
                "n_qubits_whole_graph": group[0]["n_qubits_whole_graph"], "n_zones": group[0]["n_zones"],
                "total_active_terms_mean": round(statistics.mean(r["total_active_terms"] for r in group), 1),
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                "cx_min": min(cx_vals), "cx_max": max(cx_vals),
                "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
                "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
                "worst_feasible_mass_min": round(min(r["worst_feasible_mass"] for r in group), 4),
            })
    with open(SUMMARY_CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        writer.writeheader()
        writer.writerows(summary_rows)
    print(f"wrote {SUMMARY_CSV_PATH}")

    print("\nSummary:")
    for row in summary_rows:
        print(f"  {row}")


if __name__ == "__main__":
    main()
