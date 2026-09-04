"""Extends `truncated_witness_cap_sweep.py`'s finding (lowering
`max_witness_size` cuts circuit cost far more than it costs safety) to
the two places it hadn't been checked yet:

1. The graph family the actual README/doc claims are built on --
   `generate_feeder_graph_long_range_ties`, not the density-failure
   family the first sweep used. Different failure mechanism (long-range
   ties vs. density), so the tradeoff could plausibly look different.
2. Real scale, not just `n_nodes=8` -- `n_nodes` up to 150, matching
   `run_bounded_witness_safety_survey.py`'s own range, using the same
   walked-sample build method (full enumeration is intractable at this
   scale).

Same measurement as before: transpiled CX/depth (circuit cost) and
`leakage_trace.final_feasible_mass` on a sample of real starting trees
(actual safety), at each (n_nodes, max_witness_size) point.

Writes results/truncated_witness_cap_sweep_longrange_results.csv.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from graphs import generate_feeder_graph_long_range_ties
from leakage_trace import final_feasible_mass
from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
BETA = 0.37

N_NODES_RANGE = [10, 60, 150]
SEEDS = [0, 1, 2]
MAX_WITNESS_SIZE_RANGE = [2, 3, 4, 5, 6]
K_TIES = 5
WALK_STEPS = 400
EXACT_SEARCH_MAX_SIZE = 2
N_STARTING_TREES = 100
UNSAFE_TOL = 1e-6


def run_one(n_nodes: int, seed: int, max_witness_size: int) -> dict:
    graph = generate_feeder_graph_long_range_ties(n_nodes, K_TIES, seed)
    walked = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    truncated = build_truncated_witness_mixer(
        graph, walked, max_witness_size=max_witness_size,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
    )
    construction = truncated.construction

    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx_count, depth = 0, 0

    rng = np.random.default_rng(seed * 1000 + max_witness_size)
    unsafe = 0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1

    return {
        "n_nodes": n_nodes, "seed": seed, "max_witness_size": max_witness_size,
        "n_qubits": graph.n_edges, "n_walked_trees": len(walked),
        "n_terms": len(construction.terms),
        "n_approximate_terms": truncated.n_approximate_terms,
        "mean_term_leakage_rate": round(truncated.mean_leakage_rate, 4),
        "fully_connected": construction.fully_connected,
        "cx_count": cx_count, "depth": depth,
        "n_starting_trees_surveyed": N_STARTING_TREES,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "worst_feasible_mass": round(float(min(masses)), 4),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for n_nodes in N_NODES_RANGE:
        for seed in SEEDS:
            for max_witness_size in MAX_WITNESS_SIZE_RANGE:
                row = run_one(n_nodes, seed, max_witness_size)
                rows.append(row)
                print(f"  {row}")

    with open(RESULTS_DIR / "truncated_witness_cap_sweep_longrange_results.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"wrote {RESULTS_DIR / 'truncated_witness_cap_sweep_longrange_results.csv'}")

    print("\nMean over seeds, per (n_nodes, max_witness_size):")
    for n_nodes in N_NODES_RANGE:
        for max_witness_size in MAX_WITNESS_SIZE_RANGE:
            group = [r for r in rows if r["n_nodes"] == n_nodes and r["max_witness_size"] == max_witness_size]
            mean_cx = sum(r["cx_count"] for r in group) / len(group)
            mean_unsafe = sum(r["unsafe_rate"] for r in group) / len(group)
            mean_mass = sum(r["mean_feasible_mass"] for r in group) / len(group)
            all_connected = all(r["fully_connected"] for r in group)
            print(f"  n_nodes={n_nodes} cap={max_witness_size}: mean_cx={mean_cx:.0f} "
                  f"mean_unsafe_rate={mean_unsafe:.3f} mean_feasible_mass={mean_mass:.4f} "
                  f"all_connected={all_connected}")


if __name__ == "__main__":
    main()
