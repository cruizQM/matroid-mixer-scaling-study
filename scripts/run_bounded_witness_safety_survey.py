"""Measures, at real-feeder qubit counts where `graphs.enumerate_spanning_trees`
is no longer tractable, how often `truncated_mixer.py`'s bounded-witness
construction actually loses feasible probability mass on a genuine
starting-tree trajectory -- as opposed to the per-term majority-vote
leakage RATE alone (`TruncatedTermInfo.leakage_rate`), which is a property
of the abstract validity function, not of any specific run.

## Why build on the long-range-tie generator, not the nearest-tie one

`generate_feeder_graph_long_range_ties` is the family this repo already
established (`docs/circuit-validity.md`, `run_decomposition_scaling_study.py`)
breaks the EXACT whole-graph construction at real-feeder scale -- the
natural target for an approximate alternative. (`tie_density_sweep.py`,
also this branch, shows the nearest-tie family can independently be pushed
into the same regime by density alone; either family would do here, this
one is chosen to sit next to the existing decomposition study for direct
comparison.)

## Method

At each `n_nodes`:
1. Build the graph, then a `truncated_mixer.build_truncated_witness_mixer`
   term set from a random-walk-on-the-exchange-graph sample of trees
   (`random_trees.random_walk_exchange_sample`) -- not full enumeration,
   which is exactly what's intractable at this scale, and not an i.i.d.
   sample either (i.i.d. trees are essentially never a single edge-swap
   apart at real scale, see `random_trees.py`'s own docstring).
2. Separately sample `N_STARTING_TREES` INDEPENDENT, exactly-uniform
   starting trees via Wilson's algorithm (`random_trees.random_spanning_tree`)
   -- independent of the walk sample used to build the mixer, so this
   measures generalization to genuinely fresh starting points, not just
   trees the construction was built from.
3. For each starting tree, trace its exact trajectory through the full
   term sequence at fixed beta via `leakage_trace.final_feasible_mass`
   (`sparse=True` -- required at this qubit count; the dense method would
   need a `2**n_qubits`-entry statevector).

Reports the fraction of starting trees with any measurable feasible-mass
loss ("unsafe"), and the worst-case (minimum) feasible mass observed.

Writes results/bounded_witness_safety_survey.csv.
"""

from __future__ import annotations

import csv
import statistics
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from graphs import generate_feeder_graph_long_range_ties
from leakage_trace import final_feasible_mass
from random_trees import random_spanning_tree, random_walk_exchange_sample
from truncated_mixer import build_truncated_witness_mixer

import numpy as np

RESULTS_DIR = Path(__file__).resolve().parent.parent / "results"
CSV_PATH = RESULTS_DIR / "bounded_witness_safety_survey.csv"
SUMMARY_CSV_PATH = RESULTS_DIR / "bounded_witness_safety_survey_summary.csv"

N_NODES_RANGE = [10, 30, 60, 100, 150]
SEEDS_PER_SIZE = 3
K_TIES = 5  # matches the real IEEE 33-bus feeder's actual tie-switch count, and run_decomposition_scaling_study.py's convention
WALK_STEPS = 400
MAX_WITNESS_SIZE = 6
EXACT_SEARCH_MAX_SIZE = 2
N_STARTING_TREES = 300
BETA = 0.37
UNSAFE_TOL = 1e-6


def run_one(n_nodes: int, seed: int) -> dict:
    graph = generate_feeder_graph_long_range_ties(n_nodes, K_TIES, seed)

    t0 = time.perf_counter()
    walk_rng_seed = seed * 1000 + 1
    walked_trees = random_walk_exchange_sample(graph, n_steps=WALK_STEPS, seed=walk_rng_seed)
    truncated = build_truncated_witness_mixer(
        graph, walked_trees, max_witness_size=MAX_WITNESS_SIZE,
        exact_search_max_size=EXACT_SEARCH_MAX_SIZE, seed=seed,
    )
    build_elapsed = time.perf_counter() - t0

    t1 = time.perf_counter()
    rng = np.random.default_rng(seed * 1000 + 2)
    unsafe = 0
    worst_mass = 1.0
    masses = []
    for _ in range(N_STARTING_TREES):
        start_tree = random_spanning_tree(graph, rng)
        mass = final_feasible_mass(graph, truncated.construction, start_tree, BETA, sparse=True)
        masses.append(mass)
        if mass < 1.0 - UNSAFE_TOL:
            unsafe += 1
        worst_mass = min(worst_mass, mass)
    survey_elapsed = time.perf_counter() - t1

    return {
        "n_nodes": n_nodes,
        "seed": seed,
        "n_qubits": graph.n_edges,
        "n_walked_trees": len(walked_trees),
        "n_terms": len(truncated.construction.terms),
        "n_exact_terms": truncated.n_exact_terms,
        "n_approximate_terms": truncated.n_approximate_terms,
        "mean_term_leakage_rate": round(truncated.mean_leakage_rate, 4),
        "max_term_leakage_rate": round(truncated.max_leakage_rate, 4),
        "fully_connected_on_walk_sample": truncated.construction.fully_connected,
        "n_starting_trees_surveyed": N_STARTING_TREES,
        "n_unsafe": unsafe,
        "unsafe_rate": round(unsafe / N_STARTING_TREES, 4),
        "worst_feasible_mass": round(worst_mass, 4),
        "mean_feasible_mass": round(float(np.mean(masses)), 4),
        "build_elapsed_s": round(build_elapsed, 1),
        "survey_elapsed_s": round(survey_elapsed, 1),
    }


def main() -> None:
    RESULTS_DIR.mkdir(exist_ok=True)
    rows = []
    for n_nodes in N_NODES_RANGE:
        for seed in range(SEEDS_PER_SIZE):
            row = run_one(n_nodes, seed)
            rows.append(row)
            print(f"  {row}")

    fieldnames = list(rows[0].keys())
    with open(CSV_PATH, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nwrote {CSV_PATH}")

    summary_rows = []
    for n_nodes in N_NODES_RANGE:
        group = [r for r in rows if r["n_nodes"] == n_nodes]
        summary_rows.append({
            "n_nodes": n_nodes,
            "n_seeds": len(group),
            "n_qubits_mean": round(statistics.mean(r["n_qubits"] for r in group), 1),
            "unsafe_rate_mean": round(statistics.mean(r["unsafe_rate"] for r in group), 4),
            "unsafe_rate_min": round(min(r["unsafe_rate"] for r in group), 4),
            "unsafe_rate_max": round(max(r["unsafe_rate"] for r in group), 4),
            "worst_feasible_mass_min": round(min(r["worst_feasible_mass"] for r in group), 4),
            "mean_feasible_mass_mean": round(statistics.mean(r["mean_feasible_mass"] for r in group), 4),
            "mean_term_leakage_rate_mean": round(statistics.mean(r["mean_term_leakage_rate"] for r in group), 4),
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
