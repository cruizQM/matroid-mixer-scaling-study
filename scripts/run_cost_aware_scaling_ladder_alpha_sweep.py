"""Follow-up to run_cost_aware_scaling_ladder.py: how does pushing
`cost_alpha` HIGHER than the current default (0.01) affect the four-
condition escalating-realism ladder? The earlier calibration sweep
(truncated_mixer_search_refinement.py, density family only) found
cost_alpha=0.05 occasionally misfired on hard instances (one k_ties=10
seed dropped to mean feasible mass 0.90) -- this checks whether that
risk shows up on the actual ladder (all four conditions, real scale),
not just the small density-only calibration instances, and whether
pushing cost even harder buys further cost reduction or just erodes
safety for little additional gain.

Same method, same conditions/sizes/seeds as run_cost_aware_scaling_ladder.py
-- only cost_alpha varies, swept at {0.03, 0.1} in addition to the
existing alpha=0.01 baseline already in
results/cost_aware_scaling_ladder_results.csv (that data is not
re-generated here; this script's summary references it directly for the
0.01 row of each comparison).

Writes results/cost_aware_scaling_ladder_alpha_sweep_results.csv (per
condition, size, seed, alpha) and
results/cost_aware_scaling_ladder_alpha_sweep_summary.csv (per
condition, size, alpha: mean/std).
"""

from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer
from run_cost_aware_scaling_ladder import (
    CONDITIONS, N_NODES_RANGE, SEEDS_PER_SIZE, WALK_STEPS, MAX_WITNESS_SIZE,
    EXACT_SEARCH_MAX_SIZE, N_STARTING_TREES, UNSAFE_TOL, BETA,
)

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_alpha_sweep_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_alpha_sweep_summary.csv"

ALPHA_VALUES = [0.03, 0.1]


def run_one(condition: str, n_nodes: int, seed: int, cost_alpha: float) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)

    walked = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    truncated = build_truncated_witness_mixer(
        graph, walked, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed, cost_alpha=cost_alpha,
    )
    construction = truncated.construction
    if construction.terms:
        qc = mixer_circuit(construction, beta=BETA)
        tqc = transpile(qc, basis_gates=TRANSPILE_BASIS, optimization_level=1)
        op_counts = tqc.count_ops()
        cx_count, depth = op_counts.get("cx", 0), tqc.depth()
    else:
        cx_count, depth = 0, 0

    rng = np.random.default_rng(seed * 1000 + 2)
    unsafe = 0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1

    return {
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "cost_alpha": cost_alpha,
        "k_ties": k_ties, "n_qubits": graph.n_edges,
        "n_terms": len(construction.terms), "n_approximate_terms": truncated.n_approximate_terms,
        "fully_connected": construction.fully_connected,
        "cx_count": cx_count, "depth": depth,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "worst_feasible_mass": round(float(min(masses)), 4),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for condition in CONDITIONS:
        for n_nodes in N_NODES_RANGE:
            for cost_alpha in ALPHA_VALUES:
                for seed in range(SEEDS_PER_SIZE):
                    row = run_one(condition, n_nodes, seed, cost_alpha)
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
            for cost_alpha in ALPHA_VALUES:
                group = [r for r in rows if r["condition"] == condition and r["n_nodes"] == n_nodes and r["cost_alpha"] == cost_alpha]
                cx_vals = [r["cx_count"] for r in group]
                summary_rows.append({
                    "condition": condition, "n_nodes": n_nodes, "cost_alpha": cost_alpha,
                    "k_ties": group[0]["k_ties"], "n_qubits": group[0]["n_qubits"],
                    "cx_mean": round(statistics.mean(cx_vals), 1),
                    "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                    "cx_min": min(cx_vals), "cx_max": max(cx_vals),
                    "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
                    "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
                    "worst_feasible_mass_min": round(min(r["worst_feasible_mass"] for r in group), 4),
                    "all_connected": all(r["fully_connected"] for r in group),
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
