"""Reshapes the question `run_scaling_study_log_ties.py` raised (headline
result 1's flat-tie-count assumption is fragile -- cost trends invert
under a mild, realistically-calibrated growth model) into a positive
question about the OTHER construction this repo now has:
`truncated_mixer.py`'s cost-aware (`cost_alpha=0.01`) bounded-witness
mixer, the one actually meant for real-scale, real-topology use.

Does IT stay cheap -- roughly constant, or at worst slowly growing, with
small seed-to-seed variance -- across an escalating ladder of
increasingly realistic assumptions, where the EXACT construction either
gets expensive (headline result 1's own fixed-vs-log finding) or fails
outright (tie_density_sweep.py's density failure, docs/circuit-validity.md's
real-feeder range failure)?

## The ladder, each step strictly harder/more realistic than the last

Both axes are varied independently, both calibrated to the same single
real anchor point this repo has (5 ties at the 33-bus feeder,
`k_ties(33) = 5`) -- not chosen to produce any particular result:

- **tie placement**: short-range (`generate_feeder_graph`, nearest-tie --
  headline result 1's family) vs. long-range
  (`generate_feeder_graph_long_range_ties` -- the family
  docs/circuit-validity.md's real-feeder section and this repo's
  bounded-witness safety survey are built on).
- **tie-count growth**: log (`round(1.43*ln(n))`, mild) vs. linear
  (`round(0.152*n)`, aggressive -- reaches 23 ties at n_nodes=150,
  combined with long-range placement the single hardest condition
  tested anywhere in this repo).

Four conditions: short+log, short+linear, long+log, long+linear.

## Method (identical across all four conditions and to
run_bounded_witness_safety_survey.py, for direct comparability)

At each size, builds a walked-exchange-graph sample of trees
(`random_trees.random_walk_exchange_sample` -- full enumeration is
intractable at these densities, particularly the linear-growth
conditions), then `truncated_mixer.build_truncated_witness_mixer` at its
current default (`cost_alpha=0.01`). Reports transpiled CX/depth (mean
AND spread across seeds -- low variance is part of the claim being
tested, not just a low mean) and, via `leakage_trace.final_feasible_mass`
on 100 independent Wilson's-algorithm starting trees per seed, the same
safety metrics already established (`unsafe_rate`, `mean_feasible_mass`).

Writes results/cost_aware_scaling_ladder_results.csv (per condition, per
size, per seed) and results/cost_aware_scaling_ladder_summary.csv (per
condition, per size: mean/std/min/max).
"""

from __future__ import annotations

import csv
import math
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import numpy as np
from qiskit import transpile

from graphs import generate_feeder_graph, generate_feeder_graph_long_range_ties
from leakage_trace import final_feasible_mass
from measure import TRANSPILE_BASIS
from mixer import mixer_circuit
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_results.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "cost_aware_scaling_ladder_summary.csv"
BETA = 0.37

LOG_C = 1.4300  # k_ties_log(33) == 5, same calibration as run_scaling_study_log_ties.py
LINEAR_C = 0.1515  # k_ties_linear(33) == 5, same single real anchor point

N_NODES_RANGE = [10, 30, 60, 100, 150]  # matches run_bounded_witness_safety_survey.py's own range
SEEDS_PER_SIZE = 3
WALK_STEPS = 400
MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 2
N_STARTING_TREES = 100
UNSAFE_TOL = 1e-6


def k_ties_log(n_nodes: int) -> int:
    return max(1, round(LOG_C * math.log(n_nodes)))


def k_ties_linear(n_nodes: int) -> int:
    return max(1, round(LINEAR_C * n_nodes))


CONDITIONS = {
    "short_log": (generate_feeder_graph, k_ties_log),
    "short_linear": (generate_feeder_graph, k_ties_linear),
    "long_log": (generate_feeder_graph_long_range_ties, k_ties_log),
    "long_linear": (generate_feeder_graph_long_range_ties, k_ties_linear),
}


def run_one(condition: str, n_nodes: int, seed: int) -> dict:
    generator, k_ties_fn = CONDITIONS[condition]
    k_ties = k_ties_fn(n_nodes)
    graph = generator(n_nodes, k_ties, seed)

    t0 = time.perf_counter()
    walked = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=seed * 1000 + 1)
    truncated = build_truncated_witness_mixer(
        graph, walked, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
        adaptive=True,  # explicit: this script measures the ADAPTIVE construction specifically
        # (build_truncated_witness_mixer's own default reverted to adaptive=False/cost_alpha=0.01
        # after validation found adaptive doesn't reliably beat it -- see docs/bounded-witness-mixer.md).
    )
    build_elapsed = time.perf_counter() - t0

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
        "condition": condition, "n_nodes": n_nodes, "seed": seed, "k_ties": k_ties,
        "n_qubits": graph.n_edges, "n_walked_trees": len(walked),
        "n_terms": len(construction.terms), "n_approximate_terms": truncated.n_approximate_terms,
        "fully_connected": construction.fully_connected,
        "cx_count": cx_count, "depth": depth,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "worst_feasible_mass": round(float(min(masses)), 4),
        "build_elapsed_s": round(build_elapsed, 1),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    print("k_ties by condition/size:")
    for n in N_NODES_RANGE:
        print(f"  n_nodes={n}: log={k_ties_log(n)} linear={k_ties_linear(n)}")
    print()

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
            cx_vals = [r["cx_count"] for r in group]
            depth_vals = [r["depth"] for r in group]
            summary_rows.append({
                "condition": condition, "n_nodes": n_nodes, "k_ties": group[0]["k_ties"],
                "n_qubits": group[0]["n_qubits"],
                "cx_mean": round(statistics.mean(cx_vals), 1),
                "cx_std": round(statistics.pstdev(cx_vals), 1) if len(cx_vals) > 1 else 0.0,
                "cx_min": min(cx_vals), "cx_max": max(cx_vals),
                "depth_mean": round(statistics.mean(depth_vals), 1),
                "depth_std": round(statistics.pstdev(depth_vals), 1) if len(depth_vals) > 1 else 0.0,
                "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
                "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
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
