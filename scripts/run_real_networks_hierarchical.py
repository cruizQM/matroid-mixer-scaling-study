"""Runs this repo's actual, current-best construction -- hierarchical
(density-aware) decomposition + the cost-aware bounded-witness mixer,
`run_hierarchical_decomposed_ladder.py` -- directly on real topology,
not just synthetic proxies calibrated to real anchor points.

## The two real networks used, and why not a third

`real_feeders.load_ieee33` (case33bw, 33 buses, 5 ties -- already this
repo's headline result 2 validation instance) and the newly added
`real_feeders.load_cigre_mv` (CIGRE MV benchmark, 15 buses, 3 ties -- a
second, independently-sourced, differently-sized real-topology check).

A third candidate, `pandapower.networks.mv_oberrhein` (179 buses, a real
German MV network), was investigated and found NOT directly usable: it
is structurally two separate substations' worth of sub-feeders (2
ext_grids), and its "tie" lines are not simple independent redundant
loops the way case33bw's are -- removing all 6 of them splits the
177-tree-edge remainder into 69 components, not 1, meaning several ties
are jointly load-bearing for connectivity in a way this repo's
single-tree-plus-independent-ties model doesn't represent. Its raw
tie-count and tie-span statistics were still usable (see
docs/scaling-ladder-and-decomposition.md's real-topology section) since
those don't require treating it as one clean feeder -- but running the
actual circuit CONSTRUCTION on it would need real data-cleaning work
(e.g. splitting it into its two actual service areas) beyond this
check's scope, the same call this repo already made for the CATS
transmission dataset (real_feeders.py's own docstring).

## Method

For each real network: the EXACT construction (`build_matroid_mixer`,
full enumeration -- both networks are small enough, and fully
deterministic -- no seed dependence, so measured once), AND the
hierarchical, density-aware decomposition + cost-aware bounded-witness
mixer (this branch's current-best construction). Both measured with the
same transpilation/safety methodology as the rest of this investigation.

## Why `decomposed` and `cost_capped` are averaged over multiple seeds
(and `exact` isn't)

An earlier version of this script measured `decomposed` and
`cost_capped` at `seed=0` only, matching neither this repo's own
synthetic-ladder convention (always >=3 seeds) nor, it turned out, real
seed-independence. Checked directly: even though CIGRE MV's zones are
small enough for fully deterministic exact spanning-tree enumeration,
`build_truncated_witness_mixer`'s witness search for approximate terms
still uses randomized restarts internally, and different seeds converge
to different-quality witnesses for the same term -- CIGRE MV's flat
`decomposed` result ranged from 3,668 to 9,178 CX across 5 tested seeds,
a 2.5x spread, with `seed=0` (the only one previously measured) landing
on the cheapest end by chance, not because it's representative.
`cost_capped` turned out to be robust on both networks regardless (its
`cost_alpha`-sweep already searches for the cheapest option, which
apparently collapses the space of viable witnesses enough that restart
order stops mattering) -- but that robustness was never verified either,
just assumed from a single seed. Now both are measured across
`SEEDS_PER_NETWORK` seeds and reported as mean +/- range, exactly like
the synthetic ladder always has been.

Writes results/real_networks_hierarchical_results.csv (per seed) and
results/real_networks_hierarchical_summary.csv (mean/std/min/max per
network + method).
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
from mixer import build_matroid_mixer, mixer_circuit, verify_all_terms_no_leakage
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree
from real_feeders import load_cigre_mv, load_ieee33
from run_decomposed_cost_aware_ladder import measure_subproblem
from run_cost_capped_decomposition import decompose_to_threshold, CX_THRESHOLD
from zone_decomposition import build_assembly_graph, build_zone_subgraph, partition_zones_by_size

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "real_networks_hierarchical_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "real_networks_hierarchical_summary.csv"
BETA = 0.37
N_STARTING_TREES = 200
UNSAFE_TOL = 1e-6
SEEDS_PER_NETWORK = 5


def measure_exact(graph) -> dict:
    trees = enumerate_spanning_trees(graph)
    construction = build_matroid_mixer(graph, trees)
    leak_free = verify_all_terms_no_leakage(construction) if construction.terms else True
    qc = mixer_circuit(construction, beta=BETA)
    tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
    op_counts = tqc.count_ops()

    return {
        "method": "exact_whole_graph",
        "seed": None,
        "n_terms": len(construction.terms),
        "max_witness_size": max((t.control_count for t in construction.terms), default=0),
        "dropped_candidates": construction.dropped_candidates,
        "fully_connected": construction.fully_connected,
        "leak_free_verified": leak_free,
        "cx_count": op_counts.get("cx", 0),
        "depth": tqc.depth(),
        "unsafe_rate": 0.0,
        "mean_feasible_mass": 1.0,
    }


def measure_decomposed(graph, seed: int = 0) -> dict:
    """Flat (one-level) zone decomposition, using the SAME
    `measure_subproblem` the escalating ladder validated (per-subproblem
    exact-or-sampled tree enumeration, `cost_alpha=0.01`, fixed not
    adaptive -- see that function's own docstring for why). Both real
    networks here (15 and 33 buses) are small enough that a single level
    of decomposition (target_zone_size=8, this repo's own convention)
    already gives an assembly graph too small to need the recursive step
    -- checked directly below, not assumed, via `n_zones` (2 for CIGRE
    MV, 4-5 for IEEE33) and the resulting assembly graph's own node
    count, both well under `RECURSE_THRESHOLD=12`."""
    zone_of = partition_zones_by_size(graph, 8)
    n_zones = len(set(zone_of.values()))
    total_cx = 0
    total_depth = 0
    total_terms = 0
    total_active = 0
    all_masses = []
    for zid in sorted(set(zone_of.values())):
        zg = build_zone_subgraph(graph, zone_of, zid)
        r = measure_subproblem(zg, seed)
        total_cx += r["cx"]
        total_depth += r["depth"]
        total_terms += r["n_terms"]
        total_active += r["active_terms"]
        all_masses.extend(r["masses"])

    assembly_graph, _ = build_assembly_graph(graph, zone_of, n_zones)
    r = measure_subproblem(assembly_graph, seed)
    total_cx += r["cx"]
    total_depth += r["depth"]
    total_terms += r["n_terms"]
    total_active += r["active_terms"]
    all_masses.extend(r["masses"])

    unsafe = sum(1 for m in all_masses if m < 1.0 - UNSAFE_TOL)
    return {
        "method": "decomposed",
        "seed": seed,
        "n_terms": total_terms,
        "max_witness_size": None,
        "dropped_candidates": 0,
        "fully_connected": None,
        "leak_free_verified": None,
        "cx_count": total_cx,
        "depth": total_depth,
        "n_zones": n_zones,
        "assembly_n_nodes": assembly_graph.n_nodes,
        "unsafe_rate": round(unsafe / len(all_masses), 4) if all_masses else 0.0,
        "mean_feasible_mass": round(float(np.mean(all_masses)), 4) if all_masses else 1.0,
    }


def measure_cost_capped(graph, seed: int = 0) -> dict:
    """Wraps `run_cost_capped_decomposition.decompose_to_threshold`
    (developed and validated on the synthetic ladder, where every seed at
    every size met `CX_THRESHOLD=500`) -- this is the first time it's run
    on real topology rather than a synthetic proxy."""
    result = decompose_to_threshold(graph, seed)
    unsafe = sum(1 for m in result["masses"] if m < 1.0 - 1e-6)
    return {
        "method": "cost_capped",
        "seed": seed,
        "n_terms": None,
        "max_witness_size": None,
        "dropped_candidates": 0,
        "fully_connected": None,
        "leak_free_verified": None,
        "cx_count": result["cx"],
        "depth": result["depth"],
        "n_zones": None,
        "assembly_n_nodes": None,
        "unsafe_rate": round(unsafe / len(result["masses"]), 4) if result["masses"] else 0.0,
        "mean_feasible_mass": round(float(np.mean(result["masses"])), 4) if result["masses"] else 1.0,
        "n_leaves": result["n_leaves"],
        "max_leaf_cx": result["max_leaf_cx"],
        "met_threshold": result["all_leaves_met_threshold"],
    }


def run(name: str, graph, run_exact: bool = True) -> list:
    print(f"=== {name}: n_nodes={graph.n_nodes} n_qubits={graph.n_edges} k_ties={graph.k_ties} ===", flush=True)
    rows = []
    if run_exact:
        exact = measure_exact(graph)
        exact["network"] = name
        rows.append(exact)
        print(f"  exact: {exact}", flush=True)
    else:
        # IEEE33's whole-graph exact construction is ALREADY measured in this
        # repo (results/real_feeder_results.csv, from run_real_feeder_validation.py)
        # and known to be slow (brute-force-first witness search failing
        # expensively on many of its long-range-tie candidates before
        # dropping them) -- re-running it here would just re-derive an
        # already-established, already-cited result at real cost, not add
        # new information. That existing measurement: 24 terms, 573 of 597
        # candidates dropped, fully_connected=False, cx_count=96 (cheap only
        # because it's mostly incomplete, not because it's a good mixer).
        # Fully deterministic (exact enumeration, no randomized search) --
        # no seed dependence to average over.
        print("  exact: SKIPPED -- already measured in results/real_feeder_results.csv "
              "(24 terms, 573/597 candidates dropped, fully_connected=False)", flush=True)

    for seed in range(SEEDS_PER_NETWORK):
        decomp = measure_decomposed(graph, seed=seed)
        decomp["network"] = name
        rows.append(decomp)
        print(f"  decomposed seed={seed}: {decomp}", flush=True)

    for seed in range(SEEDS_PER_NETWORK):
        capped = measure_cost_capped(graph, seed=seed)
        capped["network"] = name
        rows.append(capped)
        print(f"  cost_capped seed={seed}: {capped}", flush=True)

    return rows


def summarize(rows: list) -> list:
    """One row per (network, method): mean/std/min/max CX and mean
    feasible mass across SEEDS_PER_NETWORK seeds -- exact_whole_graph is
    deterministic (single measurement, no seed loop), reported as-is with
    std=0."""
    summary = []
    seen = set()
    for r in rows:
        key = (r["network"], r["method"])
        if key in seen:
            continue
        seen.add(key)
        group = [x for x in rows if x["network"] == r["network"] and x["method"] == r["method"]]
        cx_vals = [g["cx_count"] for g in group]
        mass_vals = [g["mean_feasible_mass"] for g in group]
        unsafe_vals = [g["unsafe_rate"] for g in group]
        summary.append({
            "network": r["network"], "method": r["method"], "n_seeds": len(group),
            "cx_mean": round(statistics.mean(cx_vals), 1),
            "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
            "cx_min": min(cx_vals), "cx_max": max(cx_vals),
            "mean_feasible_mass_mean": round(statistics.mean(mass_vals), 4),
            "unsafe_rate_mean": round(statistics.mean(unsafe_vals), 4),
        })
    return summary


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    rows.extend(run("CIGRE_MV", load_cigre_mv(), run_exact=True))
    rows.extend(run("IEEE33", load_ieee33(), run_exact=False))

    fieldnames = ["network", "method", "seed", "n_terms", "max_witness_size", "dropped_candidates",
                  "fully_connected", "leak_free_verified", "cx_count", "depth",
                  "n_zones", "assembly_n_nodes", "unsafe_rate", "mean_feasible_mass",
                  "n_leaves", "max_leaf_cx", "met_threshold"]
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    summary_rows = summarize(rows)
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
